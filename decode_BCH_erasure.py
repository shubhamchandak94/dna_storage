import numpy as np
import random
import argparse
import json
import subprocess
import sys
import bchlib
import binascii
BCH_POLYNOMIAL = 285


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
	bch = bchlib.BCH(BCH_POLYNOMIAL, data['BCH_bits'])
	corrected_reads = []
	# put corrected reads into corrected reads and drop non-corrected reads (no consensus)
	for read in data['symbols']:
		read_bytes = binascii.unhexlify(((hex(int(read[1],2)))[2:-1]).zfill(2*(data['symbol_size']+data['BCH_bits'])))
		read_data,read_ecc = read_bytes[:-bch.ecc_bytes], read_bytes[-bch.ecc_bytes:]
		(bitflips, read_data, read_ecc) = bch.decode(read_data, read_ecc)
		if bitflips >= 0: #success
			corrected_str = bin(int(binascii.hexlify(read_data), 16))[2:].zfill(data['symbol_size']*8)
			corrected_reads.append([read[0], corrected_str])
	output_data = dict(data)
	output_data['symbols'] = corrected_reads
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
