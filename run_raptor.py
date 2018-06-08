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
    return parser


def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    
    encode_script="./rq --debug encode "
    sample_script="sample_generation_raptor.py"
    decode_script="./rq --debug decode "
    data_dir = "data"

    num_chunks=config.num_chunks
    data_block=32 #in terms of bytes
    num_experiments=config.num_experiments
    input_file = config.input_file
    
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    log_file="log_" + str(num_chunks) + "_chunks.csv"
    log_data = [];
    log_data.append(["alpha","coverage","success-percentage"]);

    for alpha in [0.5]:
        
        output_file="data/output_" + str(alpha) + "_alpha.txt"
            
        ### Generate the codewords        
        arg_string  = " -s " + str(32)
        arg_string += " -m " + str(1000000) 
        arg_string += " --repair-symbols-rate " + str(alpha)
	arg_string += " " 
        arg_string += input_file + " "
        arg_string += output_file
         
        encode_command = encode_script + arg_string
        subprocess.call([encode_command], shell=True) 
        assert os.path.isfile(output_file),"The codebook did not get generated"
   	optimal_coverage = (1+alpha)*np.log(1+1/alpha)
	coverage_list = np.around(optimal_coverage + np.array([0.1]),2)
        for coverage in coverage_list:
            num_success = 0
            for iter in range(num_experiments): 
                sample_file="data/sample_" + str(alpha) + "_alpha_" + str("%.2f" % round(coverage,2)) + "coverage.txt"

                ### Generate the sample files
                arg_string  = " --output_file " + output_file
                arg_string += " --sample_file " + sample_file
                arg_string += " --coverage " + str(coverage) 
                arg_string += "  --num_chunks " + str(num_chunks)

                sample_command = "python " + sample_script + arg_string
                subprocess.call([sample_command] , shell=True)
                assert os.path.isfile(sample_file),"The sample did not get generated"

                ### Perform decoding
                recon_file="data/recon_" + str(alpha) + "_alpha_" + str("%.2f" % round(coverage,2)) + "_coverage.txt"
                ### Generate the codewords        
                arg_string  = sample_file + " "
                arg_string += recon_file
         
                decode_command = decode_script + arg_string
                ret = subprocess.call([decode_command], shell=True) 
                print "The process exited: ", ret
                if (ret == 0):
                    ret = filecmp.cmp(input_file,recon_file);
                    num_success += ret
                    print "filecmp is: ", ret

            success_percentage = int(num_success*100.0/num_experiments)
            log_data.append([alpha,str("%.2f" % round(coverage,2)),success_percentage]);
         
        writer = csv.writer(open(log_file, 'w'))
        writer.writerows(log_data);

if __name__ == '__main__':
    main()
