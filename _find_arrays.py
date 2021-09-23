import os
import re
import sys
import subprocess
import pickle
import binascii
import yaml

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from dit.divergences import jensen_shannon_divergence

from multiprocessing import Pool
from functools import partial
from itertools import groupby

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))
SPACER_LENGTH_MIN = 15
SPACER_LENGTH_MAX = 40
SPACER_OVERHANG = 15

def parse_assembly_summary(file):

    with open(file) as f:
        table = pd.read_table(f, sep='\t', index_col=0, header=0, lineterminator='\n')

    return table


def annotate_array(array, table):

    # print('.', file=sys.stderr, end='')

    gbk_id = '_'.join([array['filename'].split('_')[0], array['filename'].split('_')[1]])
    assembly = table.loc[gbk_id]

    array.update({'assembly_info': assembly})

    return array


def parse_CRISPRDetect_reports(directory):

    arrays = []

    for filename in os.listdir(directory):
        try:
            statinfo = os.stat(os.path.join(ROOT_PATH, 'crispr_arrays_p', filename))

            if statinfo.st_size > 0 and filename.endswith('.txt'):

                with open(os.path.join(ROOT_PATH, 'crispr_arrays_p', filename), 'r') as f:
                    for line in f:

                        # ------------ start of an array
                        if 'Predicted by CRISPRDetect' in line:
                            array = {}
                            array['spacers'] = []

                            array['filename'] = filename
                            array['array_n'] = line.split(' ')[1]
                            array['pos'] = line.split(' ')[2]

                        # ------------ description
                        if line.startswith('>'):
                            array['description'] = line.strip()

                        # ------------ array type
                        if line.startswith('# Array family'):
                            array['type'] = line.split(':')[1].strip()

                        # ------------ extract spacer
                        pattern = re.compile("^\d+\t")
                        if pattern.match(line.replace(' ', '')):
                            line_split = line.replace(' ', '').strip().split('\t')

                            if len(line_split) != 5 :
                                spacer_fields = ['position', 'repeat', 'identity', 'spacer', 'repeat_seq', 'spacer_seq', 'indel']

                                spacer = dict(zip(spacer_fields, line_split))

                                if spacer['spacer_seq'] != '|':
                                    array['spacers'].append(spacer)

                        # ------------ direction
                        if 'Final direction:' in line and 'HIGH' in line:
                            array['direction_prob'] = 'high'
                        elif 'Final direction:' in line and 'HIGH' not in line:
                            array['direction_prob'] = 'low'

                        if 'Final direction:' in line:
                            array['direction'] = line.strip().split(':')[1].replace(' ', '')[0]

                        # ------------ array score
                        if '# Questionable array :' in line:
                            line_p = line.strip().split('\t')
                            array['questionable'] = line_p[0].split(':')[1].replace(' ', '')
                            array['score'] = line_p[1].split(':')[1].replace(' ', '')

                        # ------------ end of an array
                        if line.startswith('//'):
                            arrays.append(array)
        except:
            print('Nop', filename)
    return arrays


def run_blast(array, database):

    print('.', file=sys.stderr, end='')

    mismatches = sys.argv[2]

    reported_aln = 10
    genome_index_path = database

    spacers_filtered = []

    temp_spacer_file = 'tmp/' + array['filename']
    blast_results_file = 'blast_results/' + array['filename'] + '_' + array['array_n'] + '.xml'

    if len(array['spacers']) > 0:

        with open(temp_spacer_file, 'wb') as f:

            for i, spacer in enumerate(array['spacers']):
                if len(spacer['spacer_seq']) >= SPACER_LENGTH_MIN and len(spacer['spacer_seq']) <= SPACER_LENGTH_MAX:
                    spacers_filtered.append(spacer['spacer_seq'])

                    f.write('>{}_{}_{}\n'.format(array['filename'], array['array_n'], i))
                    f.write('{}\n'.format(spacer['spacer_seq']))


        blastx_cline = NcbiblastxCommandline(cmd='blastn',
                                            db=database,
                                            outfmt=5,
                                            word_size=17,
                                            out=blast_results_file,
                                            query=temp_spacer_file,
                                            evalue=0.1)

        # print(blastx_cline)
        stdout, stderr = blastx_cline()

    return array


def parse_blast_xml(array):

    array['targets'] = []

    blast_results_file = 'blast_results/' + array['filename'] + '_' + array['array_n'] + '.xml'

    if os.path.exists(blast_results_file) and os.path.getsize(blast_results_file) > 0:
        with open(blast_results_file, 'r') as r:

            blast_records = NCBIXML.parse(r)

            targets = []

            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:

                        mismatches = int(record.query_length) - int(hsp.positives)

                        m_matches = hsp.match.count('|')
                        m_length = len(hsp.match)

                        m_diff = m_length - m_matches


                        if record.gapped == 0 and mismatches < 6:

                            # print('{}\t{}\t{}\t{}\t{}'.format(int(record.query_length), m_matches, int(hsp.positives), hsp.match, mismatches))

                            target_name = alignment.hit_def.split(' ')[0]

                            if hsp.frame[1] == -1:
                                b_start = hsp.sbjct_end
                                b_end = hsp.sbjct_start
                                b_strand = '-'
                            else:
                                b_start = hsp.sbjct_start
                                b_end = hsp.sbjct_end
                                b_strand = '+'

                            target = {
                                'strand': b_strand,
                                'target': target_name,
                                'start': int(b_start),
                                'seq': hsp.sbjct,
                                'mismatches': mismatches
                            }

                            targets.append(target)

            array['targets'] = targets

    return array


def run_bowtie(array, database):

    print('.', file=sys.stderr, end='')

    array['targets'] = []

    mismatches = sys.argv[2]

    reported_aln = 10
    genome_index_path = database

    spacers_filtered = []

    for spacer in array['spacers']:
        if len(spacer['spacer_seq']) >= SPACER_LENGTH_MIN and len(spacer['spacer_seq']) <= SPACER_LENGTH_MAX:
            spacers_filtered.append(spacer['spacer_seq'])

    spacers = ",".join(spacers_filtered)

    # run bowtie alignment
    bowtie_command = 'bowtie -a -v {} -l 23 --suppress 6 {} -c {}'.format(mismatches, genome_index_path, spacers)
    rproc = subprocess.Popen(bowtie_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = rproc.communicate()

    targets = []

    for line in stdout.split(os.linesep):
        if line:
            b_id, b_strand, b_target, b_start, b_seq, b_reps, b_mm = line.split('\t')

            target = {
                'strand': b_strand,
                'target': b_target,
                'start': int(b_start),
                'reps': int(b_reps),
                'seq': b_seq,
                'mismatches': b_mm
            }

            targets.append(target)

    array['targets'] = targets

    return array


def chunk_it(spacer_records, chunk_size):
    return [spacer_records[i:i + chunk_size] for i in range(0, len(spacer_records), chunk_size)]


def extract_flanking_regions(array, database):

    # print('.', file=sys.stderr, end='')
    prime_end=3

    flanks = []
    database = '{}.fna'.format(database)

    for target in array['targets']:

        strand = target['strand']

        if prime_end == 3:
            if strand == '+':
                start = target['start'] + len(target['seq'])
                end = start + 15

            elif strand == '-':
                start = target['start'] - SPACER_OVERHANG
                end = target['start']

        elif prime_end == 5:
            if strand == '+':
                start = target['start'] - SPACER_OVERHANG
                end = target['start']

            elif strand == '-':
                start = target['start'] + len(target['seq'])
                end = start + SPACER_OVERHANG

        target_seq_id = target['target']

        # run samtools faidx to extract flanks
        extract_command = 'samtools faidx {} {}:{}-{}'.format(database, target_seq_id, start, end)
        rproc = subprocess.Popen(extract_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = rproc.communicate()

        # print(stdout, stderr)

        if len(stdout.split('\n')[:2]) == 2:
            flank_id, flank_seq = stdout.split('\n')[:2]

            if strand == '-':
                flanks.append(str(Seq(flank_seq, generic_dna).reverse_complement()))
            else:
                flanks.append(flank_seq)

    array.update({'flanks': flanks})

    return array


def create_logo(name, flanks, folder):

    print('---- weblogo ' + name, file=sys.stderr)

    filename = os.path.join(ROOT_PATH, 'analysis', 'flanks', 'blast', 'flanks_' + folder, name + '_flanks.txt').replace(' ', '_')
    weblogo_save = os.path.join(ROOT_PATH, 'analysis', 'logos', 'blast', 'logos_' + folder, name + '_flanks.pdf')

    if not os.path.exists(os.path.join(ROOT_PATH, 'analysis', 'flanks', 'blast', 'flanks_' + folder)):
        os.makedirs(os.path.join(ROOT_PATH, 'analysis', 'flanks', 'blast', 'flanks_' + folder))

    if not os.path.exists(os.path.join(ROOT_PATH, 'analysis', 'logos', 'blast', 'logos_' + folder)):
        os.makedirs(os.path.join(ROOT_PATH, 'analysis', 'logos', 'blast', 'logos_' + folder))

    weblogo_command = 'weblogo --format PDF -c classic --size large'

    with open(filename, "wb") as f:
        for i, flank in enumerate(flanks):
            if len(flank) > 5:
                f.write(flank + '\n')

    f = open(filename, "r")
    print('----------- ' + filename, file=sys.stderr)

    with open(weblogo_save, 'w') as w:
        p = subprocess.Popen(weblogo_command.split(), stdin=f, stdout=w)
        p.wait()
        w.flush()

    f.close()


def filter_and_group_flanks(arrays, n_threshold, folder_name):

    folder_name = folder_name.split('/')[-1]
    arrays = sorted(arrays, key= lambda x: x['assembly_info'].taxid)
    tax_flanks = {}
    tax_names = {}

    for key, group in groupby(arrays, lambda x: x['assembly_info'].taxid):
        for array in group:

            if key not in tax_flanks.keys():
                tax_flanks[key] = []
                tax_names[key] = ''

            if array['direction_prob'] == 'high' and array['questionable'] != 'YES' and array['type'] != 'NA':
                tax_flanks[key].extend(array['flanks'])
                tax_names[key] = array['assembly_info'].organism_name

    for taxid, flanks in tax_flanks.items():
        if len(set(flanks)) > n_threshold:
            # print(tax_names[taxid], len(set(flanks)))
            create_logo(tax_names[taxid], set(flanks), folder_name)


def calculate_jsd(column, nrows, index):
    """ Calculate Jensen-Shannon divergence """
    P = dit.ScalarDistribution(['A', 'C', 'T', 'G'], [column['A'], column['C'], column['T'], column['G']])
    Q = dit.ScalarDistribution(['A', 'C', 'T', 'G'], [0.25, 0.25, 0.25, 0.25])
    
    return jensen_shannon_divergence([P, Q])


def calculate_jsd_2(column_p, column_q):
    """ Calculate Jensen-Shannon divergence """
    
    P = dit.ScalarDistribution(['A', 'C', 'T', 'G'], [column_p['A'], column_p['C'], column_p['T'], column_p['G']])
    Q = dit.ScalarDistribution(['A', 'C', 'T', 'G'], [column_q['A'], column_q['C'], column_q['T'], column_q['G']])
    
    return np.sqrt(jensen_shannon_divergence([P, Q]))


# SCORING FUNCTIONS FOR BITWISE INFORMATION
def pams_and_confidence():
    for item in species_type_dict:
        # figure out whether to use back or front flanks
        if species_type_dict[item].find("-II-") == -1:
            flank_type = "_front"
            five_prime = True
        else:
            flank_type = "_back"
            five_prime = False
        cur_species = item + flank_type
        metrics = get_metrics(all_pams[cur_species], five_prime)
        final_pams_scores[cur_species] = metrics


def get_metrics(flank_list, five_prime):
    # species-wide variables
    N = len(flank_list)
    ni = len(flank_list[0])
    en = 1.039721/N
    # i specific variables
    fai_matrix = list()
    Ri_list = list()
    for i in range(0,ni):
        freq = {"A":0,"T":0,"C":0,"G":0}
        for flank in flank_list:
            freq[flank[i]] += 1
        # get relative frequencies for each nucleotide
        for item in freq:
            freq[item] = freq[item]/N
        fai_matrix.append(freq)
        Ri_list.append(get_ri(freq,en))
    bby_pam = get_pam(Ri_list,fai_matrix,five_prime)
    return bby_pam


def get_ri(freq_table,en):
    hi_score = 0
    for item in freq_table:
        # to make sure the log doesn't see a zero:
        if freq_table[item] != 0:
            hi_score += freq_table[item]*math.log(freq_table[item],2)
    ri_score = 2+hi_score-en
    return ri_score


def get_pam(Ris,fais,is_five_prime):
    pam_nt_and_pos = dict()
    risdev = list()
    # Find the significant positions and put them into the sig_pos container
    sig_pos = list()
    Ri_mean = np.mean(Ris)
    Ri_stdev = np.std(Ris)
    for i in range(len(Ris)):
        # print(i, Ris[i], Ri_mean, Ri_stdev)
        if Ris[i] > Ri_mean+Ri_stdev/2:
            # Check if it is close to sequence PAM is unlikely past first 5 nucleotides
            if is_five_prime:
                if i < len(Ris)-5:
                    if Ris[i] > Ri_mean + Ri_stdev:
                        sig_pos.append(i)
                        risdev.append(i)
                else:
                    sig_pos.append(i)
                    risdev.append(Ris[i] - Ri_mean)
            else:
                if i > 5:
                    if Ris[i] > Ri_mean + Ri_stdev:
                        sig_pos.append(i)
                        risdev.append(i)
                else:
                    sig_pos.append(i)
                    risdev.append(Ris[i]-Ri_mean)
    # check to see if there were any significant sequences at all:
    if not sig_pos:
        return "No PAM identified."
    for pos in sig_pos:
        for nt in fais[pos]:
            # single consensus nucleotide
            if fais[pos][nt] > 0.5:
                pam_nt_and_pos[pos] = nt
            # possible degenerate code
            elif fais[pos][nt] > 0.25:
                if pos in pam_nt_and_pos:
                    pam_nt_and_pos[pos] += nt
                else:
                    pam_nt_and_pos[pos] = nt

    # find the start and end of the pam depending on the type:
    if is_five_prime:
        pam_start = min(pam_nt_and_pos.keys())
        pam_end = len(Ris)
    else:
        pam_start = 0
        pam_end = max(pam_nt_and_pos.keys())+1

    # Make the PAM sequence by putting N's in the places that are not significant
    PAM = str()
    for i in range(pam_start,pam_end):
        if i in pam_nt_and_pos:
            if len(pam_nt_and_pos[i]) > 1:
                PAM += degenerate_nucleotides(pam_nt_and_pos[i])
            else:
                PAM += pam_nt_and_pos[i]
        else:
            PAM += "N"
    return PAM, risdev, sig_pos

# ----------------------- HELPER FUNCTIONS ------------- #

def degenerate_nucleotides(nts):
    if nts == "AG":
        return "R"
    elif nts == "TC":
        return "Y"
    elif nts == "CG":
        return "S"
    elif nts == "AT":
        return "W"
    elif nts == "TG":
        return "K"
    elif nts == "AC":
        return "M"
    elif nts.find("A") == -1:
        return "B"
    elif nts.find("T") == -1:
        return "V"
    elif nts.find("C") == -1:
        return "D"
    elif nts.find("G") == -1:
        return "H"
    else:
        return "N"


def get_pam_metrics(flank_block):

    pam, metrics, sig_pos = get_metrics(flank_block, False)
    block = pd.Series(flank_block).apply(lambda x: pd.Series(list(x.upper())))
    ncols = len(block.columns)
    nrows = len(block.index)
    f_table = pd.DataFrame(np.zeros((4, ncols)), index=['A', 'C', 'T', 'G'])

    non_sig_pos = list(set(list(block.columns)) - set(sig_pos))

    jsd_table = pd.Series(index=range(0,10), dtype=float)
    jsd_table_old = pd.Series(index=range(0,10), dtype=float)

    sig_block = block[sig_pos]
    non_sig_block = block[non_sig_pos]
    bg_freq = dict(pd.Series(non_sig_block.values.flatten()).value_counts() / len(non_sig_block.values.flatten()))

    for col in block:
        f_table[col] = block[col].value_counts() / nrows
        f_table[col].fillna(0, inplace=True)

        jsd_table[col] = calculate_jsd_2(f_table[col], bg_freq)
        jsd_table_old[col] = calculate_jsd(f_table[col], 1, 1)


    jsd_sig_mean = jsd_table[sig_pos].mean()
    jsd_non_sig_mean = jsd_table[non_sig_pos].mean()

    jsig = np.sqrt(jsd_table_old[sig_pos]).mean()
    jnonsig = np.sqrt(jsd_table_old[non_sig_pos]).mean()

    return pam, metrics, sig_pos, jsd_sig_mean, jsd_non_sig_mean, jsig, jnonsig

# -----------------------------------
# main

viral_databases = [
    '/mnt/data1/blast/viral/viral.genomic',
    # '/mnt/data1/blast/mycobacteriophages/mycobacteriophages',
    # '/mnt/data1/blast/aclame/aclame_db',
    # '/mnt/data1/blast/camera/camera',
]

for vd in viral_databases:

    n_mismatches = sys.argv[2]

    db_name = vd.split('/')[-1]
    pickle_dir = os.path.join(ROOT_PATH, 'analysis', 'pickled_arrays', db_name)
    # pickle_dir = os.path.join(ROOT_PATH, 'analysis', 'pickled_arrays_blast', db_name)

    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir)

    pickle_path = os.path.join(pickle_dir, '{}_arrays_flanks.p'.format(db_name))

    threads = int(sys.argv[1])
    pool = Pool(threads)

    print('\n================= Viral database:{} Mismatches:{} ================='.format(vd, n_mismatches))
    print('\n----------------- Parsing CRISPRDetect reports')
    arrays = parse_CRISPRDetect_reports('/mnt/flash/crispr/PAMs/crispr_arrays_p/')

    # print('\n----------------- Running search')
    # arrays = pool.map(partial(run_blast, database=vd), arrays)
    arrays = pool.map(partial(run_bowtie, database=vd), arrays)

    # print('\n----------------- Parsing BLAST XML')
    # arrays = pool.map(parse_blast_xml, arrays)

    # for a in arrays:
    #     arrays = parse_blast_xml(a)

    print('\n----------------- Extracting flanks')
    arrays = pool.map(partial(extract_flanking_regions, database=vd), arrays)

    print('\n----------------- Annotating')
    assembly_summary = parse_assembly_summary('/mnt/flash/crispr/PAMs/tables/assembly_summary.txt')
    arrays = pool.map(partial(annotate_array, table=assembly_summary), arrays)

    print('\n----------------- Pickling')
    pickle.dump(arrays, open(pickle_path, "wb"))
    print('\n---------------------- saved: {}'.format(pickle_path))


# -----------------------------------------------------------------------------------
# FILTERING
# -----------------------------------------------------------------------------------
ROOT_PATH = '/Users/nace/imperial/crispr_pipeline/_tdata'

for mm in [0, 1, 2, 3]:
    for strand in ['+', '+-']:
        with open(os.path.join(ROOT_PATH, 'analysis', 'pickled_arrays', f'a_arrays_flanks_{mm}mm_{strand}.p'), 'rb') as f:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
        
        pdf = pd.DataFrame(p)
        pdf['strand'] = strand
        pdf['mismatches_count'] = mm
        dt_a.append(pdf)

df_a = pd.concat(dt_a)

asinfo = df_a['assembly_info'].to_dict()
df_a['asinfo'] = df_a.apply(lambda x: x['assembly_info'].to_dict(), axis=1)
df_a['species_taxid'] = df_a.apply(lambda x: asinfo[x.name]['species_taxid'], axis=1)
df_a['taxid'] = df_a.apply(lambda x: asinfo[x.name]['taxid'], axis=1)
df_a['organism_name'] = df_a.apply(lambda x: asinfo[x.name]['organism_name'], axis=1)

# get metrics for collected flanks
pam = {}
predicted_pam, metrics, sig_pos, jsd_sig_mean, jsd_non_sig_mean, jsig, jnonsig = get_pam_metrics(flanks)
pam['pam'] = predicted_pam
pam['sig_pos'] = sig_pos
pam['ris_sig'] = np.mean(metrics)
pam['jsd_sig_mean'] = jsd_sig_mean
pam['jsd_non_sig_mean'] = jsd_non_sig_mean
pam['flanks_count'] = len(flanks)