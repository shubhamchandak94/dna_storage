import subprocess
import numpy as np
import os
import csv
import filecmp
import argparse

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_experiments', type=int, default=10)
    parser.add_argument('--num_chunks', type=int, default=100)            
    return parser


def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    
    input_script="input_generation.py"
    encode_script="encode.py"
    sample_script="sample_generation.py"
    decode_script="decode.py"
    data_dir = "data"

    num_chunks=config.num_chunks
    data_block=32 #in terms of bytes
    num_experiments=config.num_experiments
    
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    log_file="log_" + str(num_chunks) + "_chunks.csv"
    log_data = [];
    log_data.append(["alpha","coverage","success-percentage"]);

    input_file="data/input.txt"
     
    ### Generate the input files
    arg_string  = "  --input_file " + input_file
    arg_string += "  --num_chunks " + str(num_chunks)

    generation_command = "python " + input_script + arg_string
    subprocess.call([generation_command] , shell=True)
    assert os.path.isfile(input_file),"The data did not get generated"

    for alpha in [0.1,0.2,0.3,0.5,0.8,1.2,1.5,2.0]:
        
        output_file="data/output_" + str(alpha) + "_alpha.txt"
            
        ### Generate the codewords        
        arg_string  = " -f " + input_file
        arg_string += " -m " + str(1000) 
        arg_string += " --gc " + str(0.5) 
        arg_string += " --rs " + str(0) 
        arg_string += " --out " + output_file
        arg_string += " --no_fasta "
        arg_string += " --alpha " + str(alpha)
        arg_string += " -l " + str(data_block) 
         
        encode_command = "python " + encode_script + arg_string
        subprocess.call([encode_command], shell=True) 
        assert os.path.isfile(output_file),"The codebook did not get generated"
    
        for coverage in [0.8,1.0,1.5,2.0,3.0,5.0]:
            num_success = 0
            for iter in range(num_experiments): 
                sample_file="data/sample_" + str(alpha) + "_alpha_" + str(coverage) + "coverage.txt"

                ### Generate the sample files
                arg_string  = " --output_file " + output_file
                arg_string += " --sample_file " + sample_file
                arg_string += " --coverage " + str(coverage) 
                arg_string += "  --num_chunks " + str(num_chunks)

                sample_command = "python " + sample_script + arg_string
                subprocess.call([sample_command] , shell=True)
                assert os.path.isfile(sample_file),"The sample did not get generated"

                ### Perform decoding
                recon_file="data/recon_" + str(alpha) + "_alpha_" + str(coverage) + "_coverage.txt"
                ### Generate the codewords        
                arg_string  = " -f " + sample_file
                arg_string += " -d " + str(4) 
                arg_string += " -n " + str(num_chunks) 
                arg_string += " --rs " + str(0) 
                arg_string += " --gc " + str(0.5) 
                arg_string += " -m " + str(1000) 
                arg_string += " --out " + recon_file
                #arg_string += " --truth " + input_file
         
                decode_command = "python " + decode_script + arg_string
                ret = subprocess.call([decode_command], shell=True) 
                print "The process exited: ", ret
                if (ret == 0):
                    ret = filecmp.cmp(input_file,recon_file);
                    num_success += ret
                    print "filecmp is: ", ret

            success_percentage = int(num_success*100.0/num_experiments)
            log_data.append([alpha,coverage,success_percentage]);
         
        writer = csv.writer(open(log_file, 'w'))
        writer.writerows(log_data);

if __name__ == '__main__':
    main()
