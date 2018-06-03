import subprocess
import numpy as np
import os
import csv
import filecmp
import argparse
import distance

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_experiments', type=int, default=20)
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--num_chunks', type=int, default=100)           
    parser.add_argument('--eps', type=float, default=0.0) 
    return parser


def main():
	parser = get_argument_parser()
	config = parser.parse_args()

	pchk_mat = " LDPC-codes/matrices/ldpc_9.pchk "
	gen_mat = " LDPC-codes/matrices/ldpc_9.gen "	   
	alpha = 0.5
	encode_script="./LDPC-codes/encode "
	sample_script="sample_generation_LDPC.py "
	decode_script="./LDPC-codes/decode "
	data_dir = "data"

	num_chunks=config.num_chunks
	num_experiments=config.num_experiments
	input_file = config.input_file.strip()
	eps = config.eps

	if not os.path.exists(data_dir):
		os.makedirs(data_dir)
	log_file="log_LDPC_" + str(num_chunks) +"_chunks_"+str(eps)+"_eps_"+str(alpha)+"_alpha.txt"
	f_log = open(log_file,'w')
	log_str = "alpha = " + str(alpha) + "\n\n"
	print(log_str)
	f_log.write(log_str)
	output_file=" data/encoded-file "
	    
	### Generate the codewords        
	arg_string = pchk_mat + gen_mat + input_file + output_file
	 
	encode_command = encode_script + arg_string
	subprocess.call([encode_command], shell=True) 
	f_enc = open(output_file.strip(),'r')
	enc_str = f_enc.read().replace('\n', '')
	f_enc.close()

	coverage_list = [2.0,2.1,2.2,2.3]
	for coverage in coverage_list:
		log_str = "Coverage = " + str(coverage) + "\n"
		print(log_str)
		f_log.write(log_str)
	    	num_success = 0
	    	for it in range(num_experiments): 
			sample_file=" data/received-file "

			### Generate the sample files
			arg_string  = " --output_file " + output_file
			arg_string += " --sample_file " + sample_file
			arg_string += " --coverage " + str(coverage) 
			arg_string += "  --num_chunks " + str(num_chunks)
			arg_string += "  --eps " + str(eps)

			sample_command = "python " + sample_script + arg_string
			subprocess.call([sample_command] , shell=True)

			### Perform decoding
			recon_file=" data/decoded-file "
			### Generate the codewords        
			arg_string = pchk_mat + sample_file + recon_file +" misc 0.0 prprp 100"
		 
			decode_command = decode_script + arg_string
			subprocess.call([decode_command], shell=True) 
			
			f_recon = open(recon_file.strip(),'r')
			recon_str = f_recon.read().replace('\n', '')
			f_recon.close()
			hamm = distance.hamming(enc_str,recon_str)
			if hamm == 0:
				num_success += 1
			log_str = "Bit errors = " + str(hamm) + " (abs), " + str(1.0*hamm/len(recon_str)) +  " (rel)\n"
			print(log_str)
			f_log.write(log_str)

		success_percentage = int(num_success*100.0/num_experiments)
		log_str = "Coverage = " + str(coverage) + ", Success percentage = " + str(success_percentage) + "%\n\n"
		print(log_str)
		f_log.write(log_str)
	f_log.close()	

if __name__ == '__main__':
    main()
