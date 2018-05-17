import numpy as np
import random
import argparse
import json
import subprocess
import sys

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--recon_file', type=str)
    return parser

def decode(input_file, recon_file):
	decode_raptor_script = "./rq --debug decode "    
	f_input = open(input_file, "r");
	data = json.load(f_input)
	f_input.close()
	reads_dict = {}
	# load reads and store according to index
	for read in data['symbols']:
		if read[0] in reads_dict:
			reads_dict[read[0]].append([int(c) for c in read[1]])
		else:
			reads_dict[read[0]] = [[int(c) for c in read[1]]]
	consensus_reads = [] #to store index and consensus bitstring
	count_dict = {} #to store the number of reads for each index (for sorting later)
	# convert to numpy arrays and find consensus
	for k in reads_dict.keys():
		count_dict[k] = len(reads_dict[k])
		read_array = np.array(reads_dict[k], dtype = int)
		majority_list = np.array(np.mean(read_array,axis = 0)>0.5,dtype=int)
		consensus_reads.append([k, ''.join([str(c) for c in majority_list])])
	consensus_reads.sort(key=lambda x: count_dict[x[0]])
	output_data = dict(data)
	output_data['symbols'] = consensus_reads
	f_intermediate = open("tmpfile", "w");
	f_intermediate.write(json.dumps(output_data, sort_keys = 'False', indent=2, separators=(',', ': ')))
	f_intermediate.close()
	ret = subprocess.call([decode_raptor_script+" tmpfile "+ recon_file], shell=True)
	return ret

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    ret =  decode( config.input_file, config.recon_file);
    return ret

if __name__ == '__main__': sys.exit(main())
