# CS CM122 Project 1B: Read Mapper

This is a continuation of Project 1A. For this project, I created a read mapper that mapped read to a reference genome using a Burrows Wheeler Transform. The output is list of all mutation that were contained within the reads, exluding sequencing errors. The side of the reference genome was 1,000,000 bases.

The code for this project was run in the following format
    python3 project1b_main.py [reference genome] [reads]

For the results on Codalab, I ran
    python3 project1b_main.py project1b_1000000-data/project1b_1000000_reference_genome.fasta project1b_1000000-data/project1b_1000000_with_error_paired_reads.fasta