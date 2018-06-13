import string
import numpy as np
import random
import argparse

dna2int = {}
int2dna = {}
dna2int['A'] = 0;
dna2int['C'] = 1;
dna2int['G'] = 2;
dna2int['T'] = 3;
int2dna[0] = 'A';
int2dna[1] = 'C';
int2dna[2] = 'G';
int2dna[3] = 'T';


def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', type=str, default="output.txt")
    parser.add_argument('--sample_file', type=str, default="sample.txt")
    parser.add_argument('--coverage', type=float, default=1.0)
    parser.add_argument('--num_chunks', type=int, default=100)
    parser.add_argument('--SNP_prob', type=float, default=0.0)
    return parser


def randomly_sample_reads(input_file, sample_file,num_sample_reads,SNP_prob):
    
    f_input = open(input_file, "r");
    f_sample = open(sample_file, "w");
    
    input_lines = f_input.readlines();
    #print num_sample_reads
    for i in range(num_sample_reads):
        clean_sample = random.choice(input_lines).rstrip('\n');
	clean_sample_arr = np.array([dna2int[c] for c in clean_sample])
	noise = np.random.choice([0,1,2,3],size=len(clean_sample),p=[1.0-SNP_prob,SNP_prob/3.0,SNP_prob/3.0,SNP_prob/3.0])
	output_sample_arr = np.mod(clean_sample_arr+noise,4)
	output_sample = ''.join([int2dna[j] for j in output_sample_arr])
        f_sample.write("%s\n" %output_sample)

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    num_samples = int(config.coverage*config.num_chunks)
    
    randomly_sample_reads( config.output_file, config.sample_file, num_samples,config.SNP_prob);

if __name__ == '__main__':
    main()
