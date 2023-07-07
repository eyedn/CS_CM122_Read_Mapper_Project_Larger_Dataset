###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_functions.py
###############################################################################
import project1b_classes as cl
import datetime

# given a file of reads, create a list of read objects
def create_reads_list(file, subdiv_size, min_div_size, genome, max_errors):
    print(f'{datetime.datetime.now()}: reading and mapping reads')
    mutations_and_errors = {}
    reads_list = []
    with open(file, "r") as f:
        total_reads = sum(1 for line in f if '>' in line)
        f.seek(0)  # reset file pointer to the beginning of the file
        i = 0  # initialize counter variable
        for line in f:
            curr_line = line.strip()
            if '>' in curr_line:
                # names will only have read number and pair number
                read_name = curr_line[6:]
                continue
            # associate this sequence with the current read name
            curr_read = cl.Read(curr_line, read_name, subdiv_size, min_div_size)
            # map read to the genome
            curr_read.map_subdivisions(genome)
            curr_read.map_read(genome.sequence, max_errors)
            # record this read's mutations and errors
            record_all_mutations(curr_read, mutations_and_errors)
            # if pair reads are detected...
            if '/' in read_name and read_name.split('/')[1] == '2': 
                curr_read.assign_pair(reads_list[-1])
                reads_list[-1].assign_pair(curr_read)
            reads_list.append(curr_read)
            i += 1  # increment counter variable
            if i % 10000 == 0:  # print progress message every 100,000 reads
                print(f'{datetime.datetime.now()}: {i} of {total_reads} reads mapped')
            if i in range(20):
                print(curr_read.sequence)
                print(curr_read.division_idx)
    return reads_list, mutations_and_errors

# find the hamming distance between two strings
def calc_hamming_dist(str_1, str_2):
    h_dist = 0
    # min_len = min(len(str_1), len(str_2))
    for i in range(len(str_1)):
        # hamming dist increases when the strings are not the same at pos i
        if str_1[i] != str_2[i]:
            h_dist += 1
    return h_dist

# when mapping reads, add their recorded mutations to a dictionary of all muts
def record_all_mutations(read, all_mutations):
    # conglomerate the occurances of all mutations and errors
    if read.mutations:
        for mut in read.mutations:
            try:
                all_mutations[mut] += 1
            except:
                all_mutations[mut] = 1
    return all_mutations

# once all reads are mapped, find the mutations that are likely real mutations
def concensus_mutations(all_mutations, reads_list, av_read_len, genome_len):
    consensus = []
    coverage_est = round((len(reads_list)*av_read_len)/genome_len)
    print(f'{datetime.datetime.now()}: found {len(all_mutations)} possible mutations')
    print(f'{datetime.datetime.now()}: estimated coverage is {coverage_est}')
    # only consider mutations that occur >= the minimum accepted occurance lvl
    for mut in all_mutations:
        # score mutations if they are substitutions
        if mut[1]== 'S':
            if all_mutations[mut] >= (0.4*coverage_est):
                consensus.append(mut)
        # score mutations if they are insertions
        elif mut[1]== 'I': 
            if all_mutations[mut] >= (0.15*coverage_est):
                consensus.append(mut)
        # score mutations if they are insertions
        else: 
            if all_mutations[mut] >= (0.15*coverage_est):
                consensus.append(mut)
    print(f'{datetime.datetime.now()}: {len(consensus)} mutations confirmed')
    return consensus