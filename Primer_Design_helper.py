#!/usr/bin/python
import os, primer3

def rev_comp(seq): 
### This function takes a sequence and return a reverse and complement sequence
	comp_dict = { 'A': 'T',
			'T': 'A',
			'G': 'C',
			'C': 'G',
			'a': 't',
			't': 'a',
			'g': 'c',
			'c': 'g',
			} # Creat dictionary for complement A,T,G,C
	rev_seq = seq[::-1] # reverse the input sequence
	rev_comp_seq = '' # Creat empty string for output sequence
	for n in rev_seq: # Creat the string by adding the complement NA one by one
		rev_comp_seq += comp_dict[n]
	return rev_comp_seq

def locate_seq(seq, genomic_seq):
### This function takes the query sequence and return its starting position in the genome
	q_seq_index = genomic_seq.find(seq)
	if q_seq_index == -1:
		q_seq_index = genomic_seq.find(rev_comp(seq))
	return q_seq_index

def primer_design(seq_id, seq_temp, loc, loc2):
### This function utilizes the primer3 package, takes the sequence id, sequence template, and two position coordinates to exclude the desired amplicon region
	primer = (primer3.binding.designPrimers
	{
	'SEQUENCE_ID': 'seq_id'
	'SEQUENCE_TEMPLATE': seq_temp
	'SEQUENCE_EXCLUDE_REGION': [loc1, loc2 - loc1]
	},
	{
	'PRIMER_OPT_SIZE' : 20,
	'PRIMER_MIN_SIZE': 18,
	'PRIMER_MAX_SIZE': 24,
	'PRIMER_OPT_TM': 60.0,
	'PRIMER_MIN_TM': 57.0,
	'PRIMER_MAX_TM': 63.0,
	'PRIMER_MIN_GC': 20.0,
	'PRIMER_MAX_GC': 80.0,
	'PRIMER_MAX_POLY_X': 5,
	'PRIMER_SALT_MONOVALENT': 50.0,
	'PRIMER_DNA_CONC': 50.0,
	'PRIMER_MAX_NS_ACCEPTED': 0,
	'PRIMER_MAX_SELF_ANY': 8,
	'PRIMER_MAX_SELF_END': 3,
	'PRIMER_PAIR_MAX_COMPL_ANY': 12,
	'PRIMER_PAIR_MAX_COMPL_END': 8,
	'PRIMER_PRODUCT_SIZE_RANGE': [300, len(seq_temp)], 
	})) # primer3 generates a dictionary "primer" contains many primer information
	return [primer['PRIMER_LEFT_0_SEQUENCE'], primer['PRIMER_RIGHT_0_SEQUENCE'], primer['PRIMER_LEFT_0_TM'], primer['PRIMER_RIGHT_0_TM']]
	# The function extracts and returns the first pair of primers, and their melting temperatures
	
def amplicon_extractor(gene_name, Chr-num, gene_start, gene_end, gRNA_seq1, gRNA_seq2):
### This function takes the the gene name, chromosome coordinates, and target sgRNA sequences from the input file, and return a list containing the sequence of the amplicon and positions of the gRNA within the amplicon
	file_name = 'Danio_rerio.GRCz10_CHr.' + Chr_num + '.fa'
	file_directory = 'Users/jinliu/documents/Dre_genome'
	open_file_path = file_directory + file_name
	file_handler = open(open_file_path)
	chromosome_file = file_handler.readlines()
	chromosome_seq = chromosome_file[1].strip()
	chromosome_seq = chromosome_seq.upper()
	gene_seq = chromosome_sq[gene_start : gene_end + 1].upper()
	sgRNA_loc1 = locate_seq(gRNA_seq1, gene_seq)
	sgRNA_loc2 = locate_seq(gRNA_seq2, gene_seq)
	if sgRNA_loc1 > 150:
		location_mark_1 = min([sgRNA_loc1, sgRNA_loc2]) - 150
	else:
		location_mark_1 = 0
	if len(gene_seq) - sgRNA_loc2 > 150:
		location_mark_2 = max([sgRNA-loc1, sgRNA_loc2]) + 170
	else:
		location_mark_2 = len(gene_seq) # The amplicon contains 150bp upstream and downstream of the two target gRNA sites. So that the gRNA are away enough from primer binding sites, avoiding misreading from sanger sequencing.
	if location_mark_1 >= 200:
		amplicon_start = location_mark_1 - 200
	else:
		amplicon_start = 0
	amplicon_end = location_mark_2 + 200
	amplicon_seq = gene_seq[amplicon_start:amplicon_end] # The gene fragment for primer design template goes further 200bp upstream and downstream of the amplicon, to give primer3 enough space to pick optimal primers

	file_handler.close()
	return [amplicon_seq, amplicon_start, location_mark_1, location_mark_2, gene_seq]

input_file = open('inputfile.txt')
input_lst = input_file.readlines()
primer_file = open('primers_lst.txt', 'w')
primer_header = 'Gene_Name' + '\tForward_primer' + '\tFPrimer_TM' + '\tReverse_Primer' + '\tRprimer_TM' + '\tFragment_size' + '\n'
primer_file.write(primer_header)

for ln in input_lst:
	line = ln.split()
	Chr_num = str(line[1])
	gene_name = str(line[0])
	gene_start = int(line[2]) - 1
	gene_end = int(line[3])
	gRNA_sq1 = str(line[4])
	if len(line) == 6:
		gRNA_seq2 = str(line[5])
	else: 
		gRNA_seq2 = '0'
	amplicon = amplicon_extractor(gene_name, Chr_num, gene_start, gene_end, gRNA_seq1, gRNA_seq2)
	amplicon_seq = amplicon[0]
	primers = primer_design(gene_name, amplicon[0], amplicon[2] - amplicon[1], amplicon[3] - amplicon[1])
	primer_line = gene_name + '\t' + primers[0] + '\t' + str(primers[2]) + '\t' + str(primers[3])
	primer_file.write(primer_line)
	amplicon_file_name = gene_name + '_amplicon.seq'
	amplicon_file = open(amplicon_file_name, 'w')
	amplicon_start = locate_seq(primers[0], amplicon_seq)
	amplicon_end = locate_seq(primers[1], amplicon_seq)
	gRNA_loc1 = locate_seq(gRNA_seq1, amplicon_seq)
	gRNA_loc2 = locate_seq(gRNA_seq2, amplicon_seq)
	if gRNA_loc1 > gRNA_loc2:
		amplicon_seq_output = amplicon_seq[amplicon_start:gRNA_loc2].lower() + amplicon_seq[gRNA_loc2:gRNA_loc2 + len(gRNA_seq2)] + amplicon_seq[gRNA_loc2 + len(gRNA_seq2):gRNA_loc2].lower() + amplicon_seq[gRNA_loc1:gRNA_loc1 + len(gRNA_seq1)] + amplicon_seq[gRNA_loc1 + len(gRNA_seq1):amplicon_end + len(primers[1]) + 1].lower()
	else:
		amplicon_seq_output = amplicon_seq[amplicon_start:gRNA_loc1].lower() + amplicon_seq[gRNA_loc1:gRNA_loc1 + len(gRNA_seq1)] + amplicon_seq[gRNA_loc1 + len(gRNA_seq1):gRNA_loc2].lower() + amplicon_seq[gRNA_loc2:gRNA_loc2 + len(gRNA_seq2)] + amplicon_seq[gRNA_loc2 + len(gRNA_seq2):amplicon_end + len(primers[1]) + 1].lower() # The sgRNA sites are upper cased, while the rest are lower cased in the output file
	amplicon_file.write(amplicon_seq_output)
	amplicon_file.close()
	amplicon_size = len(amplicon_seq_output)
	primer_output_continue = '\t' + str(amplicon_size) + '\n'
	primer_file.write(primer_output_continue)
	amplicon_cor_1 = location_seq(amplicon_seq[amplicon_start:amplicon_end + len(primers[1]) + 1], amplicon[4])
	amplicon_cor_2 = amplicon_cor_1 + amplicon_size

primer_file.close()
	
