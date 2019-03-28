import util_final as util
file_name = 'myfile_160K'
ldpc_code_prefix = 'LDPC-codes/matrices/ldpc_9'
oligo_length = 100
LDPC_alpha = 0.5
bin_index_len = 14
sub_prob = 0.005
ins_prob = 0.0005
del_prob_list = [0.01]
sync_list = ['','AGT']
file_size = 160000
num_trials = 3
frac_random_reads = 0.15
print('frac_random_reads:',frac_random_reads)

for sync in sync_list:
    for BCH_bits in [1,2,3,4]:
        for del_prob in del_prob_list:
            print('BCH_bits:',BCH_bits)
            print('sync:',sync)
            print('del_prob:',del_prob)
            index_len = (bin_index_len+6*BCH_bits)//2
            sync_pos = (oligo_length-index_len)//2+index_len
            min_coverage = util.find_min_coverage(file_name,oligo_length,BCH_bits,LDPC_alpha,ldpc_code_prefix,bin_index_len,file_size,sub_prob,sub_prob,num_trials,ins_prob = ins_prob, del_prob = del_prob, start_coverage = 2.2, sync = sync, sync_pos=sync_pos,frac_random_reads=frac_random_reads)
            print('min_coverage',min_coverage)
