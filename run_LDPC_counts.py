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
    parser.add_argument('--eps', type=float, default=0.0) 
    return parser


def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    
    # load LDPC parity check matrix
    d = scipy.io.loadmat('code-ar4ja_n=20480_c=45.mat')	
    H = d['H'] # parity check matrix
    H = scipy.sparse.lil_matrix(H)
    block_length = 20480
    #shuffle non-punctured columns of H - just to get randomness because
    #our channel works on blocks
    H = H.T
    index = np.array(range(block_length))		
    np.random.shuffle(index)
    H[:block_length,:] = H[index,:]
    H = H.T 
    m,n = np.shape(H)
    (r,c) = H.nonzero()
    check_node_nbr = []
    for i in range(m):
	check_node_nbr.append(c[np.nonzero(r == i)])
    var_node_nbr = []
    for i in range(n):
	var_node_nbr.append(r[np.nonzero(c == i)])
    alpha = 0.25
    num_experiments=config.num_experiments
    eps = config.eps
    oligo_length = 256
    num_oligos = block_length/oligo_length
    log_file="log_LDPC_0.8_20480_" +str(eps)+"_eps.csv"
    log_data = [];
    log_data.append(["alpha","coverage","success-percentage"]);

    coverage_list = [2.0,2.5,3.0,4.0,5.0]
    for coverage in coverage_list:
	    print('coverage',coverage)	
	    num_success = 0
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
		l_y = np.zeros(n)
		for i in range(num_oligos):
			for j in range(oligo_length):
				l_y[oligo_length*i+j] = (2*zero_counts[i,j]-total_counts[i])*np.log((1-eps)/eps)
                x_hat = message_passing(l_y,check_node_nbr,var_node_nbr,30) 
		if np.sum(x_hat[:block_length]) == 0:
			print('success')
		    	num_success += 1
		else:
			print('failure')

	    success_percentage = int(num_success*100.0/num_experiments)
	    log_data.append([str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),success_percentage]);
	    print([str("%.2f" % round(alpha,2)),str("%.2f" % round(coverage,2)),success_percentage]);
    writer = csv.writer(open(log_file, 'w'))
    writer.writerows(log_data);

def message_passing(l_y,check_node_nbr,var_node_nbr,num_iter):
    m = len(check_node_nbr)
    n = len(var_node_nbr)
    x_hat = np.zeros(n,dtype = int)
    msg_var_to_check = {}
    msg_check_to_var = {}
    # first message from var to check
    for i in range(n):
        for j in var_node_nbr[i]:
            msg_var_to_check[(i,j)] = l_y[i]
    for it in range(num_iter):
        # calculating messages from check to var nodes
        for i in range(m):
            for j in check_node_nbr[i]:
                prod = 1.0 
                for t in check_node_nbr[i]:
                    if t != j:
                        prod *= np.tanh(msg_var_to_check[(t,i)]/2)
                    msg_check_to_var[(i,j)] = 2*np.arctanh(prod)

        # calculating messages from var to check nodes
        for i in range(n):
            #first calculate sum of all messages
            for j in var_node_nbr[i]:
                msg_var_to_check[(i,j)] = l_y[i]
                for t in var_node_nbr[i]:
                    if t != j:
                        msg_var_to_check[(i,j)] += msg_check_to_var[(t,i)]

    #final decoding
    for i in range(n):
        sum_message = l_y[i]
        for j in var_node_nbr[i]:
            sum_message += msg_check_to_var[(j,i)]    
        if sum_message >= 0:
            x_hat[i] = 0
        else:
            x_hat[i] = 1

    return x_hat

if __name__ == '__main__':
    main()

