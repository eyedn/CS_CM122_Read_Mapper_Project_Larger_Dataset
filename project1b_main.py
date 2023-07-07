###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_main.py
###############################################################################
import project1b_classes as cl
import project1b_functions as fx
from sys import argv
import datetime

def main():
    genome_file = argv[1]
    reads_file = argv[2]
    output_file = argv[3]
    # create genome objcet and dictionary of reads objects
    print(f'{datetime.datetime.now()}: reading in genome')
    my_genome = cl.Genome(genome_file)
    # reading in reads and mapping them to the genome
    my_reads, read_mutations_and_errors = fx.create_reads_list(reads_file, 10, 7, my_genome, 3)
    # finding consensus mutations
    my_mutations = fx.concensus_mutations(read_mutations_and_errors, my_reads, 50, len(my_genome.sequence))
    with open(output_file, "w") as f:
        for mut in my_mutations:
            f.write(mut + '\n')


if __name__ == "__main__":
    main()