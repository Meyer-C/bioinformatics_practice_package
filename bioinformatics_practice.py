import math

class FASTA_seq:

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

class FASTQ_seq(FASTA_seq):

    def __init__(self, name, seq, sign, quality):
        super().__init__(name, seq)
        self.sign = sign
        self.quality = quality

def FASTA_open(FASTA_filepath):
    # opens a FASTA file and creates objects for each sequence and name
    all_FASTA = []
    if FASTA_filepath.endswith('.gz'):
        with gzip.open(FASTA_filepath) as new_FASTA_filepath:
            FASTA_readlines(FASTA)

    elif FASTA_filepath.endswith('.txt') or FASTA_filepath.endswith('.fasta'):
        FASTA_readlines(FASTA_filepath)
    with open(FASTA_filepath,'r') as FASTA:
        a = FASTA.readlines()
        # first, find all lines in a FASTA formatted .txt file
        for x in range(len(a)-1):
            line_info = a[x]
            if line_info[0] == '>':
                seq_name = line_info
                seq = ''
                for y in range(len(a)-(x+1)):
                    line = a[x+1+y]
                    if line[0] != '>':
                        seq += line.strip()
                    elif line[0] == '>':
                        break
                    elif x+1+y == len(a):
                        break
                all_FASTA.append(FASTA_seq(seq_name, seq))
    return all_FASTA

def FASTA_readlines(file):
    all_FASTA = []
    a = file.readlines()
    # first, find all lines in a FASTA formatted .txt file
    for x in range(len(a) - 1):
        line_info = a[x]
        if line_info[0] == '>':
            seq_name = line_info
            seq = ''
            for y in range(len(a) - (x + 1)):
                line = a[x + 1 + y]
                if line[0] != '>':
                    seq += line.strip()
                elif line[0] == '>':
                    break
                elif x + 1 + y == len(a):
                    break
            all_FASTA.append(FASTA_seq(seq_name, seq))

    return all_FASTA

def bp_count(seq):
    ACGT = [0,0,0,0]
    for x in range(len(seq)):
        if seq[x] == 'A':
            ACGT[0] += 1
        elif seq[x] == 'C':
            ACGT[1] += 1
        elif seq[x] == 'G':
            ACGT[2] += 1
        elif seq[x] == 'T':
            ACGT[3] += 1
    return ACGT

def t_to_u(seq):
    new_seq = ''
    for x in range(len(seq)):
        if seq[x] == 'T':
            new_seq += 'U'
        else:
            new_seq += seq[x]
    return new_seq

def reverse_complement(seq):
    # this should now be inclusive if there are RNA seqs or DNA seqs
    complement_seq = ''
    if 'U' not in seq:
        complement_base = {'A':'T','T':'A','C':'G','G':'C'}
    else:
        complement_base = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    for x in range(len(seq)):
        for base,complement in complement_base.items():
            if seq[x] == base:
                complement_seq += complement
    return complement_seq [::-1]

def high_GC_from_FASTA(FASTA_filepath):
    with open(FASTA_filepath,'r') as FASTA:
        a = FASTA.readlines()
        all_gc = []
        seq_names = []
        # first, find all lines in a FASTA formatted .txt file
        for x in range(len(a)-1):
            line_info = a[x]
            seq = ''
            if line_info[0] == '>':
                seq_names.append(line_info)
                for y in range(len(a)-(x+1)):
                    line = a[x+1+y]
                    if line[0] != '>':
                        seq += line.strip()
                    elif line[0] == '>':
                        break
                    elif x+1+y == len(a):
                        break
            if len(seq) != 0:
                all_bp = bp_count(seq)
                all_gc.append(((all_bp[1]+all_bp[2])/len(seq))*100)
        max_gc = max(all_gc)
        max_index = all_gc.index(max_gc)
        return seq_names[max_index].strip(), max_gc


def point_mutation_SL(seq1,seq2):
    # count point mutations for 2 sequences of the same length(SL)
    mut_count = 0
    for x in range(len(seq1)):
        if seq1[x] != seq2[x]:
            mut_count += 1
    return mut_count

# finish this one up, SOY CONFUSCIOS
# ---------------------------------------------------------------------------------------------------------------------
def mendel(homo_dom,hetero,homo_res):
    # put in number of homo pos, hetero, and homo neg and return the prob that 2 random selected individuals will return a hetero genotype
    probs_dom_allele = [1,0.5,0]

# ---------------------------------------------------------------------------------------------------------------------

def dna_transcribe(seq):
    transcribed_dna = ''
    for x in range(len(seq)):
        if seq[x] == 'T':
            transcribed_dna += 'U'
        else:
            transcribed_dna += seq[x]
    return transcribed_dna


# just finish the dictionary, takes too long to write it all
def rna_translate(rna_seq):
    aa_from_rna = {'UUU':'F','CUU':'L','AUU':'I','GUU':'V','UUC':'F','CUC':'L','AUC':'I','GUC':'V','UUA':'L','CUA':'L','AUA':'I','GUA':'V','UUG':'L','CUG':'L','AUG':'M','GUG':'V','UCU':'S','CCU':'P','ACU':'T','GCU':'A','UCC':'S','CCC':'P','ACC':'T','GCC':'A','UCA':'S','CCA':'P','ACA':'T','GCA':'A','UCG':'S','CCG':'P','ACG':'T','GCG':'A','UAU':'Y','CAU':'H','AAU':'N','GAU':'D','UAC':'Y','CAC':'H','AAC':'N','GAC':'D','UAA':'Stop','CAA':'Q','AAA':'K','GAA':'E','UAG':'Stop','CAG':'Q','AAG':'K','GAG':'E','UGU':'C','CGU':'R','AGU':'S','GGU':'G','UGC':'C','CGC':'R','AGC':'S','GGC':'G','UGA':'Stop','CGA':'R','AGA':'R','GGA':'G','UGG':'W','CGG':'R','AGG':'R','GGG':'G'}
    aa_seq = ''
    for x in range(int((len(rna_seq)/3))):
        codon = rna_seq[x*3:x*3+3]
        for rna,aa in aa_from_rna.items():
            if rna == codon and aa != 'Stop':
                aa_seq += aa
            elif rna == codon and aa == 'Stop':
                return aa_seq

def motif_find(motif,seq):
    motif_pos = []
    for x in range(len(seq)-len(motif)+1):
        if seq[x:x+len(motif)] == motif:
            motif_pos.append(x+1)
    return motif_pos

def long_com_motif(FASTA_filepath):
    # input a list of seqs and return longest common sub seq
    all_fasta = FASTA_open(FASTA_filepath)
    seqs = []
    for x in range(len(all_fasta)):
        seqs.append(all_fasta[x].seq)
    start_seq = seqs[0]
    motif_len = 1
    long_motifs = []
    for y in range(len(start_seq) - motif_len):
        potential_motif = start_seq[y:y+motif_len]
        for z in range(1,len(seqs)):
            comp_seq = seqs[z]
            counter = 0
            for a in range(len(comp_seq) - motif_len):
                if comp_seq[a:a+motif_len] == potential_motif:
                    counter += 0
            if counter == (len(seqs) - 1):
                long_motifs.append(potential_motif)
    return long_motifs
    # fuccc me, it hard


def compliment_DNA(seq):
    new_seq = ''
    base_dict = {'A':'T','T':'A','G':'C','C':'G'}
    for x in range(len(seq)):
        for base,comp in base_dict.items():
            if seq[x] == base:
                new_seq += comp
    return new_seq

def open_reading_frame(FASTA_filepath):
    fasta = FASTA_open(FASTA_filepath)
    prot_seqs = []
    for x in fasta:
        seq = dna_transcribe(x.seq)
        rev_seq = reverse_complement(seq)
        for z in range(len(rev_seq)-3):
            if rev_seq[z:z+3] == 'AUG':
                prot_seq = rna_translate(rev_seq[z::])
                if prot_seq != 'None':
                    prot_seqs.append(prot_seq)
        for y in range(len(seq)-3):
            codon = seq[y:y+3]
            if codon == 'AUG':
                prot_seq = rna_translate(seq[y::])
                if prot_seq != 'None':
                    prot_seqs.append(prot_seq)

    return [*set(prot_seqs)]

def calc_protein_mass(prot_seq):
    masses = {'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,'M':131.04049,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333}
    tot_mass = 0
    for aa in prot_seq:
        for prot,mass in masses.items():
            if aa == prot:
                tot_mass += mass
    return tot_mass


# okay, so hear me out, I found a new way to write for loops that is much more effiecient
#
# rather than:
#
#   for x in range(len(blank)):
#       return blank[x]
#
# use this:
#
#   for x in blank:
#       return x
#
# this is essentially the same as for key value in dict
# but instead do for whatever it is(can be named anything) in list, return element of list


# Let's find some protein mutations
def protein_seq_mutations(seq1, seq2):

    this_seq_comparison = ''

    if len(seq1) == len(seq2):
        this_seq_comparison +=  seq_comparison(seq1, seq2)

    elif len(seq1) < len(seq2):
        seq_ham_dist = {}
        for x in range(len(seq2) - len(seq1)):
            new_seq2 = seq2[x: x+len(seq1)]
            seq_ham_dist[new_seq2] = hamming_dist(seq1, new_seq2)
        this_seq_comparison = seq_comparison(min(seq_ham_dist, key=seq_ham_dist.get), seq1)

    else:
        seq_ham_dist = {}
        for x in range(len(seq1) - len(seq2)):
            new_seq1 = seq1[x: x + len(seq2)]
            seq_ham_dist[new_seq1] = hamming_dist(seq2, new_seq1)
        this_seq_comparison = seq_comparison(min(seq_ham_dist, key=seq_ham_dist.get), seq2)

    return this_seq_comparison

def hamming_dist(seq1, seq2):
    hamming_distance = 0
    for x, symbol in enumerate(seq1):
        if symbol != seq2[x]: hamming_distance += 1
    return hamming_distance

def seq_comparison(seq1, seq2):
    amino_acid_properties = {'A':'H', 'I':'H', 'L':'H', 'M':'H', 'V':'H', 'F':'A', 'W':'A', 'Y':'A', 'N':'P', 'C':'P',
                             'Q':'P', 'S':'P', 'T':'P', 'D':'+', 'E':'+', 'R':'-', 'H':'-', 'K':'-', 'G':'U', 'P':'U'}
    comparison_1 = ''
    comparison_2 = ''
    for y, amino_acid in enumerate(seq1):
        if amino_acid == seq2[y]:
            comparison_1 += '|'
            comparison_2 += '|'
        else:
            comparison_1 += amino_acid_properties.get(amino_acid)
            comparison_2 += amino_acid_properties.get(seq2[y])
    return str(f'{seq1}\n{comparison_1}\n{comparison_2}\n{seq2}')


def protein_seq_properties(seq):
    amino_acid_properties = {'A': 'H', 'I': 'H', 'L': 'H', 'M': 'H', 'V': 'H', 'F': 'A', 'W': 'A', 'Y': 'A', 'N': 'P',
                             'C': 'P', 'Q': 'P', 'S': 'P', 'T': 'P', 'D': '+', 'E': '+', 'R': '-', 'H': '-', 'K': '-',
                             'G': 'U', 'P': 'U'}
    properties_seq = ''
    for amino_acid in seq:
        properties_seq += amino_acid_properties.get(amino_acid)
    return str(f'{seq}\n{properties_seq}')

seq1 = 'AIYMAFLVGLYC'
seq2 = 'ADYMAKAVGSYC'
print(protein_seq_mutations(seq1, seq2))
print(protein_seq_properties(seq1))
print(protein_seq_properties(seq2))




















