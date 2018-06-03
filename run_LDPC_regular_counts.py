import subprocess
import numpy as np
import os
import csv
import filecmp
import argparse
import scipy.io
import scipy.sparse
from tqdm import tqdm

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_experiments', type=int, default=10)
    parser.add_argument('--num_oligos', type=int, default=100)
    parser.add_argument('--k', type=int, default=6)
    parser.add_argument('--l', type=int, default=3)
    parser.add_argument('--eps', type=float, default=0.005) 
    return parser


def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    
    k = config.k
    l = config.l
    num_oligos = config.num_oligos
    alpha = 1.0*l/(k-l) 
    num_experiments=config.num_experiments
    eps = config.eps
    oligo_length = 256
    block_length = oligo_length*num_oligos
    print("Constructing code ...")	
    (check_node_nbr,var_node_nbr) = generate_code(l,k,block_length)
    print("Code construction successful")
    log_file="log_LDPC_regular_k_"+str(k)+"_l_"+str(l)+"_numoligos_"+str(num_oligos)+"_eps_" +str(eps)+".csv"
    log_data = [];
    log_data.append(["alpha","coverage","block_error","bit_error"]);

    coverage_list = [1.8]
    for coverage in coverage_list:
	    print('coverage',coverage)	
	    num_success = 0
	    bit_errors = 0	
	    for iter in range(num_experiments): 
                # perform experiment
                total_counts = np.zeros(num_oligos)
		zero_counts = np.zeros((num_oligos,oligo_length))
		for i in range(int(coverage*num_oligos/(1+alpha))):
			rid = np.random.randint(num_oligos)
			total_counts[rid] += 1
			zero_counts[rid,:]  += (np.random.rand(oligo_length) > eps)
		print('Unique reads: ',np.sum(total_counts!=0))
		#generate log likelihood vector
		l_y = np.zeros(block_length)
		for i in range(num_oligos):
			for j in range(oligo_length):
				l_y[oligo_length*i+j] = (2*zero_counts[i,j]-total_counts[i])*np.log((1-eps)/eps)
                x_hat = message_passing(l_y,check_node_nbr,var_node_nbr,30) 
                bit_errors += np.sum(x_hat)
		if np.sum(x_hat) == 0:
			print('success')
		    	num_success += 1
		else:
			print('failure')
			print('number of errors', np.sum(x_hat))

	    success_percentage = int(num_success*100.0/num_experiments)
	    log_data.append([str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),100-success_percentage,1.0*bit_errors/block_length]);
	    print([str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),success_percentage,1.0*bit_errors/block_length]);
    writer = csv.writer(open(log_file, 'w'))
    writer.writerows(log_data);

def generate_code(l,k,n):
    m = int(l*n/k)
    check_node_nbr = np.zeros((m,k),dtype = int)
    var_node_nbr = np.zeros((n,l), dtype= int)
    while(True):
        temp_array = np.array([[i]*k for i in range(m)])
        temp_array = temp_array.flatten()
        np.random.shuffle(temp_array)
        duplicate_found = False
        for i in range(n):
            if(np.size(np.unique(temp_array[l*i:l*(i+1)])) != l):#duplicates
                duplicate_found = True
                break
        if(duplicate_found == True):        
            continue
        var_node_nbr = np.reshape(temp_array,(n,l))
        for i in range(m):
            check_node_nbr[i,:] = np.where(var_node_nbr == i)[0]
        break
    return (check_node_nbr,var_node_nbr)    

def message_passing(l_y,check_node_nbr,var_node_nbr,num_iter):
    (m,k) = np.shape(check_node_nbr)
    (n,l) = np.shape(var_node_nbr)
    x_hat = np.zeros(n,dtype = int)
    msg_var_to_check = {}
    msg_check_to_var = {}
    # first message from var to check
    for i in range(n):
        for j in var_node_nbr[i,:]:
            msg_var_to_check[(i,j)] = l_y[i]
    for it in range(num_iter):
        # calculating messages from check to var nodes
        for i in range(m):
            for j in check_node_nbr[i,:]:
                prod = 1.0 
                for t in check_node_nbr[i,:]:
                    if t != j:
                        prod *= np.tanh(msg_var_to_check[(t,i)]/2)
                    msg_check_to_var[(i,j)] = 2*np.arctanh(prod)

        # calculating messages from var to check nodes
        for i in range(n):
            #first calculate sum of all messages
            for j in var_node_nbr[i,:]:
                msg_var_to_check[(i,j)] = l_y[i]
                for t in var_node_nbr[i,:]:
                    if t != j:
                        msg_var_to_check[(i,j)] += msg_check_to_var[(t,i)]

    #final decoding
    for i in range(n):
        sum_message = l_y[i]
        for j in var_node_nbr[i,:]:
            sum_message += msg_check_to_var[(j,i)]    
        if sum_message >= 0:
            x_hat[i] = 0
        else:
            x_hat[i] = 1

    return x_hat

if __name__ == '__main__':
    main()
