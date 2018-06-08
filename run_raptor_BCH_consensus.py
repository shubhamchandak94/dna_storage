import subprocess
import numpy as np
import os
import csv
import filecmp
import argparse

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
    
    encode_script="python encode_BCH_raptor.py "
    sample_script="sample_generation_raptor.py"
    decode_script="python decode_BCH_consensus.py "
    data_dir = "data"

    num_chunks=config.num_chunks
    num_experiments=config.num_experiments
    input_file = config.input_file
    eps = config.eps

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    log_file="log_BCH_" + str(num_chunks) +"_chunks_"+str(eps)+"_eps.csv"
    log_data = [];
    log_data.append(["alpha_raptor","BCH_bits","alpha","coverage","success-percentage"]);

    for alpha_raptor in [0.28]:
        for BCH_bits in [4]:
		#calculate effective alpha and num_chunks
		read_length = (31-BCH_bits) - (31-BCH_bits)%4 + BCH_bits
		alpha = (1+alpha_raptor)*float(read_length)/(read_length-BCH_bits) - 1
		num_chunks_eff = int(32.0*num_chunks/read_length)
		# effective number of chunks in original file, to calculate coverage
		
		output_file="data/output.txt"
		    
		### Generate the codewords        
		arg_string = " --input_file " + input_file
		arg_string += " --output_file " + output_file
		arg_string += " --BCH_bits " + str(BCH_bits)
		arg_string += " --alpha_raptor " + str(alpha_raptor)
		 
		encode_command = encode_script + arg_string
		subprocess.call([encode_command], shell=True) 
		assert os.path.isfile(output_file),"The codebook did not get generated"
		
		coverage_list = [2.5]
		for coverage in coverage_list:
		    num_success = 0
		    for iter in range(num_experiments): 
			sample_file="data/sample.txt"

			### Generate the sample files
			arg_string  = " --output_file " + output_file
			arg_string += " --sample_file " + sample_file
			arg_string += " --coverage " + str(coverage) 
			arg_string += "  --num_chunks " + str(num_chunks_eff)
			arg_string += "  --eps " + str(eps)

			sample_command = "python " + sample_script + arg_string
			subprocess.call([sample_command] , shell=True)
			assert os.path.isfile(sample_file),"The sample did not get generated"

			### Perform decoding
			recon_file="data/recon_" + str(alpha) + "_alpha_" + str("%.2f" % round(coverage,2)) + "_coverage.txt"
			### Generate the codewords        
			arg_string  = "--input_file " + sample_file + " "
			arg_string += "--recon_file " + recon_file
		 
			decode_command = decode_script + arg_string
			ret = subprocess.call([decode_command], shell=True) 
			print "The process exited: ", ret
			if (ret == 0):
			    ret = filecmp.cmp(input_file,recon_file);
			    num_success += ret
			    print "filecmp is: ", ret

		    success_percentage = int(num_success*100.0/num_experiments)
		    log_data.append([alpha_raptor,BCH_bits,str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),success_percentage]);
		    print([alpha_raptor,BCH_bits,str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),success_percentage])
         
        writer = csv.writer(open(log_file, 'w'))
        writer.writerows(log_data);

if __name__ == '__main__':
    main()
