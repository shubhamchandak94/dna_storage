import numpy as np
import random
import argparse

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', type=str, default="output.txt")
    parser.add_argument('--sample_file', type=str, default="sample.txt")
    parser.add_argument('--coverage', type=float, default=1.0)
    return parser

def randomly_sample_reads(input_file, sample_file,coverage):
    
    f_input = open(input_file, "r");
    f_sample = open(sample_file, "w");
    
    input_lines = f_input.readlines();
    num_input_lines = len(input_lines)
    num_sample_reads = int(coverage*num_input_lines)
    print num_sample_reads
    for i in range(num_sample_reads):
        f_sample.write("%s" %random.choice(input_lines));

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    randomly_sample_reads( config.output_file, config.sample_file, config.coverage);

if __name__ == '__main__':
    main()
