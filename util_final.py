import itertools
import math
import numpy as np
import bchlib
import binascii
import random
import os
import subprocess
import distance
import shutil

#some constants
LDPC_dim = 256000
index_block_len = 16
index_block_len_noRLL = 12
# we operate in the BCH field of 2^6 = 64 bits, so we need 6 bits per error correction
BCH_POLYNOMIAL = 67
BCH_bits_per_error = 6
LDPC_max_iter = 100

# parameters for PRP x -> ax+b mod 2^n
# obtained by running np.random.randint(2**n)
prp_a = {14: 15713, 16: 30627, 18: 111891, 20: 468527, 22: 3683343, 24: 3776737}
prp_b = {14: 10394, 16: 55140, 18: 12297, 20: 575746, 22: 1535258, 24: 10356263}
prp_a_inv = {14: 9889, 16: 29707, 18: 154907, 20: 123087, 22: 277231, 24: 8626977}
# modular inverse of prp_a mod 2^n (obtained by running modinv(prp_a,2**n))
# prp_a must be odd for the inversion to work

# gcd and modular inverse from https://stackoverflow.com/questions/4798654/modular-multiplicative-inverse-function-in-python
# These were used to create a pseudorandom permutation over 0 ... 2^n-1 to avoid issues with sequences of homopolymers
# in the same LDPC block

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

# below function from https://codereview.stackexchange.com/questions/151329/reverse-complement-of-a-dna-string
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

# functions for reading Kalign consensus
def one_hot(row, alphabet):
    _temp = np.zeros((len(row),alphabet))
    _temp[np.arange(len(row)),np.array(row)] = 1.0
    return _temp

def read_in_kalign_fasta(file_name):
    f = open(file_name)
    reads = []
    for i, line in enumerate(f):
        if i%2 == 1:
            reads.append(line.rstrip('\n'))
    f.close()
    return reads

def get_read_counts(reads):
    mapping_dict = {'A': 0, 'C': 1, 'G':2, 'T':3 , '-':4}
    reads_id = [[mapping_dict[c] for c in read] for read in reads]
    reads_one_hot = np.array([one_hot(read_id,5) for read_id in reads_id])
    read_counts = np.sum(reads_one_hot,axis=0)
    return read_counts

def filter_counts(read_counts, payload_length, sync, sync_pos):
    # returns numpy arrays total_count and zero_count of size 2*payload_length each
    # returns None if failed to recover any part of the read
    payload_length_after_sync = payload_length + len(sync)
    int2base = {0:'A',1:'C',2:'G',3:'T'}
    total_counts = np.zeros(2*payload_length)
    zero_counts = np.zeros(2*payload_length)
    pos_in_count_array = 0

    failure_flag = False
    if len(read_counts) < payload_length_after_sync:
        failure_flag = True # failed
    elif len(read_counts) > payload_length_after_sync:
        read_count_w_index = [(i,read_counts[i]) for i in range(len(read_counts))]
        read_count_w_index.sort(key=lambda a: a[1][4]) # sort by increasing number of blanks
        if read_count_w_index[payload_length_after_sync][1][4] == read_count_w_index[payload_length_after_sync-1][1][4]:
            # bad: here we are unable to precisely separate and so we fail
            failure_flag = True
        threshold = read_count_w_index[payload_length_after_sync-1][1][4]
        # we'll keep positions where the number of blanks is leq to the threshold
    else:
        threshold = np.sum(read_counts[0])
        # chosen so that all positions are picked
    if not failure_flag:
        num_done = 0
        for row in read_counts:
            if row[4] <= threshold:
                if sync != '' and (num_done in [sync_pos+i for i in range(len(sync))]):
                    num_done += 1
                    continue # no need to add sync positions to LDPC counts
                num_done += 1
                total_counts[pos_in_count_array] += np.sum(row[:4])
                total_counts[pos_in_count_array+1] += np.sum(row[:4])
                zero_counts[pos_in_count_array] += row[0] + row[1] # A and C have 0 at first bit
                zero_counts[pos_in_count_array+1] += row[0] + row[2] # A and G have 0 at second bit
                pos_in_count_array += 2
        assert pos_in_count_array == 2*payload_length
    else:
        if sync == '':
            return None
        else:
            # try to recover using sync (see if part before or after sync marker is likely to be good)
            consensus = ''.join([int2base[np.argmax(row[:4])] for row in read_counts])
            sync_pos_from_end = payload_length_after_sync-sync_pos
            if len(consensus) >= sync_pos+len(sync) and consensus[sync_pos:sync_pos+len(sync)] == sync:
                pos_in_count_array = 0
                pos_read_counts_start = 0
                pos_read_counts_end = sync_pos
            elif len(consensus) >= sync_pos_from_end and consensus[-sync_pos_from_end:-sync_pos_from_end+len(sync)] == sync:
                pos_in_count_array = 2*sync_pos
                pos_read_counts_start = -sync_pos_from_end+len(sync)
                pos_read_counts_end = -1
            else:
                return None
            for row in read_counts[pos_read_counts_start:pos_read_counts_end]:
                total_counts[pos_in_count_array] += np.sum(row[:4])
                total_counts[pos_in_count_array+1] += np.sum(row[:4])
                zero_counts[pos_in_count_array] += row[0] + row[1] # A and C have 0 at first bit
                zero_counts[pos_in_count_array+1] += row[0] + row[2] # A and G have 0 at second bit
                pos_in_count_array += 2
    return (total_counts, zero_counts)

def bin2dna_2bpb(bin_string):
	'''
	Generate dna according to mapping 00-A,01-C,10-G,11-T.
	Throw error if binary string length is not even.
	'''
	if len(bin_string)%2 != 0:
		raise Exception('binary string length not even')
	d = {'00':'A','01':'C','10':'G','11':'T'}
	return ''.join([d[bin_string[2*i:2*(i+1)]] for i in range(len(bin_string)/2)])


def dna2bin_2bpb(dna_string):
	'''
	Generate binary according to mapping A-00,C-01,G-10,T-11.
	'''
	d = {'A':'00','C':'01','G':'10','T':'11'}
	return ''.join([d[dna_string[i]] for i in range(len(dna_string))])

def binary_string_to_bytes(bin_string):
	'''
	Convert binary string to bytes (needed for BCH).
	Throw error if bin_string length not multiple of 8.
	'''
	if len(bin_string)%8 != 0:
		raise Exception('binary string length not multiple of 8')
	return binascii.unhexlify((((hex(int(bin_string,2)))[2:]).rstrip('L')).zfill(len(bin_string)/4))

def bytes_to_binary_string(byte_string):
	'''
	Convert bytes to binary string (needed for BCH).
	'''
	return bin(int(binascii.hexlify(byte_string), 16))[2:].zfill(len(byte_string)*8)

def add_index(num_oligos,BCH_bits,infile_name,outfile_name):
	'''
	Generate DNA encoding of all indexes from 0 to num_oligos-1, each protected with BCH_bits protection.
	The index DNA strings are appropriately concatenated with corresponding lines from infile_name and the resulting oligos are written to outfile_name, one line per index.
	Throw error if num_oligos > 2**24 = 16777216.
	'''
	if num_oligos > MAX_OLIGOS_NO_RLL:
		raise Exception('Too many oligos')
	block_len = index_block_len_noRLL
	bin_block_len = 2*block_len

	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = BCH_bits*BCH_bits_per_error/2
		num_bases_index = block_len + num_bases_BCH
	index = 0
	with open(infile_name) as infile, open(outfile_name, 'w') as outfile:
		for line in infile:
			index_prp = (prp_a*index+prp_b)%MAX_OLIGOS_NO_RLL
			bin_string = bin(index_prp)[2:].zfill(bin_block_len)
			dna = bin2dna_2bpb(bin_string)
			if BCH_bits != 0:
				bits_ecc = bytes_to_binary_string(bch.encode(binary_string_to_bytes(bin_string)))
				bits_ecc = bits_ecc[:BCH_bits*BCH_bits_per_error]
				dna_ecc = bin2dna_2bpb(bits_ecc)
				outfile.write(dna+dna_ecc+line)
			else:
				outfile.write(dna+line)
			index += 1
			if index == num_oligos:
				break

def remove_index(num_oligos,BCH_bits,infile_name,outfile_data,outfile_index,mode="correct",attempt_indel_cor=False):
	'''
	Decode index from a collection of (noisy) reads in infile_name and write data and index to outfile_data and outfile_index, line by line, skipping positions where index failed to decode.
	Mode can be "correct" or "detect":
	Correct - erase if correction fails
	Detect - don't try to correct, erase if any error detected
	'''
	if num_oligos > 2**24:
		raise Exception('Too many oligos')
	block_len = index_block_len_noRLL
	bin_block_len = 2*block_len
	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = BCH_bits*BCH_bits_per_error/2
		num_bases_index = block_len + num_bases_BCH
	else:
		num_bases_index = block_len
	count = 0
	deletion_corrected = 0
	f_failed = open('failed_indices','w')
	count_failed = 0
	with open(infile_name) as infile, open(outfile_data, 'w') as f_data, open(outfile_index, 'w') as f_index:
		for line in infile:
			dna_data = line[num_bases_index:]
			dna_index = line[:num_bases_index]
			if BCH_bits != 0:
				dna_ecc = dna_index[-num_bases_BCH:]
				bin_ecc = dna2bin_2bpb(dna_ecc)
				bin_ecc = bin_ecc[-BCH_bits*BCH_bits_per_error:]
				extra_len = 8*int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/8))-len(bin_ecc)
				bin_ecc = bin_ecc+'0'*extra_len
				(bitflips,cor_index,cor_ecc) = bch.decode(binary_string_to_bytes(dna2bin_2bpb(dna_index[:-num_bases_BCH])),binary_string_to_bytes(bin_ecc))
				if bitflips >= 0: #success in correction
					if mode == "detect" and bitflips > 0:
						continue
					bin_index = bytes_to_binary_string(cor_index)
					index_prp = int(bin_index,2)
					index = (prp_a_inv*(index_prp-prp_b))%MAX_OLIGOS_NO_RLL
					if index < num_oligos:
						f_data.write(dna_data)
						f_index.write(str(index)+'\n')
						count += 1
					else:
						f_failed.write('@BCH_success_but_not_within_num_oligos_'+str(count_failed)+'\n'+line)
						count_failed+=1
				else:
					if attempt_indel_cor:
					# len <= because it includes new line
					# deletion correction mode
					# try to correct 1 deletion by inserting each base at each possible position
					# and try BCH decoding with 0 errors
						dna_index = line[:num_bases_index-1]
						successful_indices = []
						ins_flag = False # found in insertion
						for pos in range(num_bases_index):
							for new_base in ['A','C','G','T']:
								new_dna_index = dna_index[:pos]+new_base+dna_index[pos:]
								dna_ecc = new_dna_index[-num_bases_BCH:]
								bin_ecc = dna2bin_2bpb(dna_ecc)
								bin_ecc = bin_ecc[-BCH_bits*BCH_bits_per_error:]
								extra_len = 8*int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/8))-len(bin_ecc)
								bin_ecc = bin_ecc+'0'*extra_len
								(bitflips,cor_index,cor_ecc) = bch.decode(binary_string_to_bytes(dna2bin_2bpb(new_dna_index[:-num_bases_BCH])),binary_string_to_bytes(bin_ecc))
								if bitflips == 0:
									successful_indices.append(cor_index)
								if len(successful_indices) > 1:
									break # fail
                                        if len(successful_indices) == 0:
                                            ins_flag = True
                                            # try insertion (i.e., delete positions one by one and see if it satisfies BCH w/o error)
                                            dna_index = line[:num_bases_index+1]
                                            for pos in range(num_bases_index):
                                                new_dna_index = dna_index[:pos]+dna_index[pos+1:]
                                                dna_ecc = new_dna_index[-num_bases_BCH:]
                                                bin_ecc = dna2bin_2bpb(dna_ecc)
                                                bin_ecc = bin_ecc[-BCH_bits*BCH_bits_per_error:]
                                                extra_len = 8*int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/8))-len(bin_ecc)
                                                bin_ecc = bin_ecc+'0'*extra_len
                                                (bitflips,cor_index,cor_ecc) = bch.decode(binary_string_to_bytes(dna2bin_2bpb(new_dna_index[:-num_bases_BCH])),binary_string_to_bytes(bin_ecc))
                                                if bitflips == 0:
                                                    successful_indices.append(cor_index)
                                                if len(successful_indices) > 1:
                                                    break # fail

                                        if len(successful_indices) == 1:
                                            # exactly one successful decoding
                                            cor_index= successful_indices[0]
                                            bin_index = bytes_to_binary_string(cor_index)
                                            index_prp = int(bin_index,2)
                                            index = (prp_a_inv*(index_prp-prp_b))%MAX_OLIGOS_NO_RLL
                                            if index < num_oligos:
                                                if ins_flag:
                                                    dna_data = line[num_bases_index+1:]
                                                else:
                                                    dna_data = line[num_bases_index-1:]
                                                f_data.write(dna_data)
                                                f_index.write(str(index)+'\n')
                                                count += 1
                                                deletion_corrected += 1
                                            else:
                                                f_failed.write('@del_corr_failed_one_found_but_not_within_num_oligos_'+str(count_failed)+'\n'+line)
                                                count_failed += 1
                                        elif len(successful_indices) == 0:
                                            f_failed.write('@del_corr_failed_none_found_'+str(count_failed)+'\n'+line)
                                            count_failed+=1
                                        else:
                                            f_failed.write('@del_corr_failed_more_than_one_found_'+str(count_failed)+'\n'+line)
                                            count_failed+=1
                                    else:
                                            f_failed.write('@BCH_failed_del_cor_not_attempted_'+str(count_failed)+'\n'+line)
                                            count_failed+=1
			else:
				bin_index =  dna2bin_2bpb(dna_index)
				if status == 0:
                                        index_prp = int(bin_index,2)
                                        index = (prp_a_inv*(index_prp-prp_b))%MAX_OLIGOS_NO_RLL
					if index < num_oligos:
						f_data.write(dna_data)
						f_index.write(str(index)+'\n')
						count += 1
        f_failed.close()
	print "Successfully decoded",count,"indices"
        print("Deletion corrected",deletion_corrected)

def rll1_pad(dna_str,padded_len):
	'''
	Pad dna_str to padded_len such that the padding does not introduce any runs of >= 2
	'''
	extra_len = padded_len - len(dna_str)
	if extra_len == 0:
		return dna_str
    if len(dna_str) == 0:

    next_symbol =   {
                    'A':['C','G','T'],
                    'C':['A','G','T'],
                    'G':['A','C','T'],
                    'T':['A','G','C']
                    }

    for i in range(extra_len):
	d = {'A':'C','C':'G','G':'T','T':'A'}
	d1 = {'A':0,'C':1,'G':2,'T':3}
	d2 = {0:'A',1:'C',2:'G',3:'T'}
	delta_coded = np.random.choice([1,2,3],extra_len)
	delta_coded[0] = d1[d[dna_str[-1]]]
	quaternary_coded = np.mod(np.cumsum(delta_coded),4)
	return dna_str+''.join([d2[i] for i in quaternary_coded])

def encode_data(infile,oligo_length,outfile,BCH_bits,LDPC_alpha,LDPC_prefix,sync='',sync_pos=-1):
	'''
	Encode binary data in infile to oligos written to oufile.
	LDPC_prefix.pchk and LDPC_prefix.gen are the LDPC matrices. LDPC_prefix.systematic contains positions of systematic bits.
	int(LDPC_dim*LDPC_alpha) should be the number of parity check bits.
	infile is a bytestream file.
	'''
        oligo_length_before_sync = oligo_length-len(sync)
	f_in = open(infile,'r')
	data = f_in.read()
	bin_data = bytes_to_binary_string(data)

	# calculate various parameters for encoding
	data_len = len(bin_data)
        num_LDPC_blocks = int(math.ceil(1.0*data_len/LDPC_dim))
        parity_bits_per_LDPC_block = int(LDPC_alpha*LDPC_dim)

	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = BCH_bits*BCH_bits_per_error/2
		num_bases_index = index_block_len_noRLL + num_bases_BCH
	else:
		num_bases_index = index_block_len_noRLL

	num_bases_payload = oligo_length_before_sync - num_bases_index
	bits_per_oligo = num_bases_payload*2
	num_oligos_data_per_LDPC_block = int(math.ceil(1.0*LDPC_dim/(bits_per_oligo)))
	num_oligos_parity_per_LDPC_block = int(math.ceil(1.0*parity_bits_per_LDPC_block/bits_per_oligo))
	num_oligos_per_LDPC_block = num_oligos_data_per_LDPC_block+num_oligos_parity_per_LDPC_block
	overall_rate = 1.0*data_len/(num_oligos_per_LDPC_block*num_LDPC_blocks*oligo_length)
	print 'overall rate:' ,overall_rate, 'bpb'
	print 'num oligos:', num_LDPC_blocks*num_oligos_per_LDPC_block
	print 'oligo length:', oligo_length, 'bases'
	print 'bases per oligo for index + index parity:', num_bases_index
	print 'fraction of oligos used for parity check:', 1.0*num_oligos_parity_per_LDPC_block/num_oligos_per_LDPC_block
	print 'number of LDPC blocks:', num_LDPC_blocks
	f_out = open(outfile+'.tmp','w')
	# pad bin_data to multiple of LDPC_dim
	extra_len = LDPC_dim*num_LDPC_blocks - len(bin_data)
	bin_data = bin_data + ''.join(np.random.choice(['0','1'],extra_len))

	# find positions of parity bits in LDPC encoded file
	f_sys = open(LDPC_prefix+".systematic",'r')
	sys_pos = f_sys.readlines()
	f_sys.close()
	sys_pos = np.array([int(i) for i in sys_pos])
	mask = np.zeros(parity_bits_per_LDPC_block+LDPC_dim,dtype=bool)
	mask[sys_pos] = True
	mask = ~mask
	parity_pos = np.nonzero(mask)[0]

        # encode LDPC code
        f_LDPC_input = open(outfile+'.tmp.1','w')
        f_LDPC_input.write(bin_data)
        f_LDPC_input.close()
        subprocess.call(["./LDPC-codes/encode "+LDPC_prefix+".pchk "+LDPC_prefix+".gen "+outfile+'.tmp.1 '+outfile+'.tmp.2 '], shell=True)
        # read parity bits
        f_LDPC_output = open(outfile+'.tmp.2','r')
	for i in range(num_LDPC_blocks):
		bin_data_block = bin_data[i*LDPC_dim:(i+1)*LDPC_dim]
#                f_LDPC_input = open(outfile+'.tmp.1','w')
#                f_LDPC_input.write(bin_data_block)
#                f_LDPC_input.close()
#		# encode LDPC
#		subprocess.call(["./LDPC-codes/encode "+LDPC_prefix+".pchk "+LDPC_prefix+".gen "+outfile+'.tmp.1 '+outfile+'.tmp.2 '], shell=True)
#		# read parity bits
#		f_LDPC_output = open(outfile+'.tmp.2','r')
		encoded_LDPC = f_LDPC_output.readline().rstrip('\n')
#		f_LDPC_output.close()
		parity_bits_block = ''.join([encoded_LDPC[j] for j in parity_pos])
		# pad to multiple of bits_per_oligo
		extra_len = bits_per_oligo*num_oligos_parity_per_LDPC_block - len(parity_bits_block)
		parity_bits_block = parity_bits_block + ''.join(np.random.choice(['0','1'],extra_len))
		# Write dna oligos (without index)
		# First write for data
		for j in range(num_oligos_data_per_LDPC_block-1):
                    f_out.write(bin2dna_2bpb(bin_data_block[j*bits_per_oligo:(j+1)*bits_per_oligo])+'\n')
		# for last oligo, might need to pad
                f_out.write(rll2_pad(bin2dna_2bpb(bin_data_block[(num_oligos_data_per_LDPC_block-1)*bits_per_oligo:]),num_bases_payload)+'\n')
		# Now write parity bits
		for j in range(num_oligos_parity_per_LDPC_block):
			f_out.write(bin2dna_2bpb(parity_bits_block[j*bits_per_oligo:(j+1)*bits_per_oligo])+'\n')
	f_out.close()
        f_LDPC_output.close()
	add_index_noRLL(num_oligos_per_LDPC_block*num_LDPC_blocks,BCH_bits,outfile+'.tmp',outfile+'.tmp.1')
        if sync != '':
            with open(outfile+'.tmp.1') as fin, open(outfile,'w') as fout:
                for line in fin:
                    fout.write(line[:sync_pos]+sync+line[sync_pos:])
            os.remove(outfile+'.tmp.1')
        else:
            os.rename(outfile+'.tmp.1',outfile)
	os.remove(outfile+'.tmp')
	os.remove(outfile+'.tmp.2')
	return

def decode_data(infile,oligo_length,outfile,BCH_bits,LDPC_alpha,LDPC_prefix,file_size, eps, mode = 'correct', MSA=False, sync='',sync_pos=-1):
	'''
	Decoder corresponding to encoder encode_data. Need same parameters as that.
	infile is file containing reads of the same length as the oligo_length.
	Returns status - 0 for success, 1 for failure
	In case of success, resulting decoded file is written to outfile.
	file_size is the size of the original file in bytes.
	mode is correct or detect as in remove_index
	eps is the error to be used in LDPC LLR
	'''
	data_len = 8*file_size
	decoded_data = []

        # NEW
	# calculate various parameters for encoding
        num_LDPC_blocks = int(math.ceil(1.0*data_len/LDPC_dim))
        parity_bits_per_LDPC_block = int(LDPC_alpha*LDPC_dim)

	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = BCH_bits*BCH_bits_per_error/2
		num_bases_index = index_block_len_noRLL + num_bases_BCH
	else:
		num_bases_index = index_block_len_noRLL
        oligo_length_before_sync = oligo_length - len(sync)
	num_bases_payload = oligo_length_before_sync - num_bases_index
	bits_per_oligo = num_bases_payload*2
	num_oligos_data_per_LDPC_block = int(math.ceil(1.0*LDPC_dim/(bits_per_oligo)))
	num_oligos_parity_per_LDPC_block = int(math.ceil(1.0*parity_bits_per_LDPC_block/bits_per_oligo))
	num_oligos_per_LDPC_block = num_oligos_data_per_LDPC_block+num_oligos_parity_per_LDPC_block
	overall_rate = 1.0*data_len/(num_oligos_per_LDPC_block*num_LDPC_blocks*oligo_length)
	print 'overall rate:' ,overall_rate, 'bpb'
	print 'num oligos:', num_LDPC_blocks*num_oligos_per_LDPC_block
	print 'oligo length:', oligo_length, 'bases'
	print 'bases per oligo for index + index parity:', num_bases_index
	print 'fraction of oligos used for parity check:', 1.0*num_oligos_parity_per_LDPC_block/num_oligos_per_LDPC_block
	print 'number of LDPC blocks:', num_LDPC_blocks

	# find positions of parity bits in LDPC encoded file
	f_sys = open(LDPC_prefix+".systematic",'r')
	sys_pos = f_sys.readlines()
	f_sys.close()
	sys_pos = np.array([int(i) for i in sys_pos])
	mask = np.zeros(parity_bits_per_LDPC_block+LDPC_dim,dtype=bool)
	mask[sys_pos] = True
	mask = ~mask
	parity_pos = np.nonzero(mask)[0]
	tmp_index_file = infile+'.tmp.index'
	tmp_data_file = infile+'.tmp.data'
	# first decode index
    remove_index(num_LDPC_blocks*num_oligos_per_LDPC_block,BCH_bits,infile,tmp_data_file,tmp_index_file,mode,attempt_indel_cor=True)
	# Now store counts in an array
	total_counts = np.zeros((num_LDPC_blocks,LDPC_dim+parity_bits_per_LDPC_block))
	zero_counts = np.zeros((num_LDPC_blocks,LDPC_dim+parity_bits_per_LDPC_block))

        if MSA == False:
            # count number of oligos with correct length after index removal
            num_correct_length = 0
            num_incorrect_length = 0
            index_set = set([])
            index_with_correct_length = set([])
            with open(tmp_data_file,'r') as infile_data, open(tmp_index_file,'r') as infile_index:
                    for line_data, line_index in zip(infile_data,infile_index):
                        index = int(line_index)
                        index_set.add(index)
                        # TODO: change this
                        if len(line_data) != num_bases_payload+1:
                            num_incorrect_length += 1
                            continue
                        index_with_correct_length.add(index)
                        num_correct_length += 1
                        block_number = index/num_oligos_per_LDPC_block
                        index_in_block = index%num_oligos_per_LDPC_block
                        if index_in_block < num_oligos_data_per_LDPC_block:
                                # data oligo
                                start_pos = index_in_block*bits_per_oligo
                                end_pos = (index_in_block+1)*bits_per_oligo
                                if end_pos > LDPC_dim:
                                        end_pos = LDPC_dim
                                total_counts[block_number][sys_pos[start_pos:end_pos]] += 1
                                bin_str = dna2bin_2bpb(line_data[:(end_pos-start_pos)/2])
                                bin_arr = np.array([int(c) for c in bin_str])
                                zero_counts[block_number][sys_pos[start_pos:end_pos]] += 1-bin_arr
                        else:
                                # parity oligo
                                start_pos = bits_per_oligo*(index_in_block-num_oligos_data_per_LDPC_block)
                                end_pos = start_pos + bits_per_oligo
                                if end_pos > parity_bits_per_LDPC_block:
                                        end_pos = parity_bits_per_LDPC_block
                                bin_str = dna2bin_2bpb(line_data[:num_bases_payload])
                                bin_arr = np.array([int(c) for c in bin_str])
                                total_counts[block_number][parity_pos[start_pos:end_pos]] += 1
                                zero_counts[block_number][parity_pos[start_pos:end_pos]] += 1-bin_arr[:end_pos-start_pos]

            print('Number of oligos with correct length after index removal', num_correct_length)
            print('Number of oligos with incorrect length after index removal', num_incorrect_length)
            print('Number of unique indices seen:', len(index_set))
            print('Number of unique indices seen with at least one correct length oligo:',len(index_with_correct_length))
        else:
            filename_index_success = 'outfile.success.index'
            filename_consensus_success = 'outfile.success.consensus'
            filename_failure = 'outfile.failure'
            filename_success = 'outfile.success'
            filename_index_failure = 'outfile.failure.index'
            fout_index_success = open(filename_index_success,'w')
            fout_consensus_success = open(filename_consensus_success,'w')
            fout_failure = open(filename_failure,'w')
            fout_success = open(filename_success,'w')
            fout_index_failure = open(filename_index_failure,'w')
            # Do MSA using Kalign
            index_set = set([])
            tmp_dir = infile+'tmp.splitted_read'
            os.mkdir(tmp_dir)
            with open(tmp_data_file,'r') as infile_data, open(tmp_index_file,'r') as infile_index:
                for line_data, line_index in zip(infile_data,infile_index):
                    index = int(line_index)
                    index_set.add(index)
                    fout = open(tmp_dir+'/'+str(index)+'.fasta','a+')
                    fout.write('>\n'+line_data)
                    fout.close()
            print('Number of unique indices:', len(index_set))
            suffix_fasta = '.fasta'
            suffix_kalign = '.kalign.fasta'
            print('Running Kalign for MSA for each index')
            for index in index_set:
                subprocess.call(['./kalign2_current/kalign '+tmp_dir+'/'+str(index)+suffix_fasta+' ' + tmp_dir+'/'+str(index)+suffix_kalign+' -quiet'],shell=True)
            num_indices_correct_len = 0
            num_reads_utilized = 0
            num_indices_wrong_len = 0
            num_reads_wasted = 0
            print('Computing counts from Kalign output')
            for index in index_set:
                if os.path.isfile(tmp_dir+'/'+str(index)+suffix_kalign):
                    reads = read_in_kalign_fasta(tmp_dir+'/'+str(index)+suffix_kalign)
                else:
                    # kalign does not create new file if only one read
                    reads = read_in_kalign_fasta(tmp_dir+'/'+str(index)+suffix_fasta)
                read_counts = get_read_counts(reads)
                ret = filter_counts(read_counts,num_bases_payload,sync,sync_pos-num_bases_index)
                if ret == None:
                    num_indices_wrong_len += 1
                    num_reads_wasted += len(reads)
                    fout_failure.write(str(index)+'\n')
                    fout_failure.write(''.join([reads[i]+'\n' for i in range(len(reads))]))
                    fout_failure.write('\n')
                    fout_index_failure.write(str(index)+'\n')
                else:
                    fout_success.write(''.join([reads[i]+'\n' for i in range(len(reads))]))
                    fout_success.write('Consensus:\n'+ret[2]+'\n\n')
                    fout_index_success.write(str(index)+'\n')
                    fout_consensus_success.write(ret[2]+'\n')
                    num_indices_correct_len += 1
                    num_reads_utilized += len(reads)
                    block_number = index/num_oligos_per_LDPC_block
                    index_in_block = index%num_oligos_per_LDPC_block
                    if index_in_block < num_oligos_data_per_LDPC_block:
                        # data oligo
                        start_pos = index_in_block*bits_per_oligo
                        end_pos = (index_in_block+1)*bits_per_oligo
                        if end_pos > LDPC_dim:
                                end_pos = LDPC_dim
                        total_counts[block_number][sys_pos[start_pos:end_pos]] += ret[0][:end_pos-start_pos]
                        zero_counts[block_number][sys_pos[start_pos:end_pos]] += ret[1][:end_pos-start_pos]
                    else:
                        # parity oligo
                        start_pos = bits_per_oligo*(index_in_block-num_oligos_data_per_LDPC_block)
                        end_pos = start_pos + bits_per_oligo
                        if end_pos > parity_bits_per_LDPC_block:
                                end_pos = parity_bits_per_LDPC_block
                        total_counts[block_number][parity_pos[start_pos:end_pos]] += ret[0][:end_pos-start_pos]
                        zero_counts[block_number][parity_pos[start_pos:end_pos]] += ret[1][:end_pos-start_pos]
            print('Number of indices with correct length consensus:',num_indices_correct_len)
            print('Number of indices with incorrect length consensus:',num_indices_wrong_len)
            print('Number of reads utilized in the counts',num_reads_utilized)
            print('Number of reads not utilized in the counts',num_reads_wasted)
            shutil.rmtree(tmp_dir)
	os.remove(tmp_index_file)
	os.remove(tmp_data_file)

	# log likelihood ratio
	llr = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
	# Now decode LDPC blocks one by one
	for i in range(num_LDPC_blocks):
		f_out = open(outfile+".tmp",'w')
		for j in range(LDPC_dim+parity_bits_per_LDPC_block):
			f_out.write(str(llr[i][j])+' ')
		f_out.close()
		subprocess.call(["./LDPC-codes/decode "+LDPC_prefix+".pchk "+outfile+'.tmp '+outfile+'.tmp.1 '+"misc 0.0 prprp "+str(LDPC_max_iter)], shell=True)
		f_in = open(outfile+'.tmp.1','r')
		LDPC_decoded_str = f_in.read().rstrip('\n')
		f_in.close()
		LDPC_decoded_bin_str = ''.join([LDPC_decoded_str[j] for j in sys_pos])
		decoded_data.append(LDPC_decoded_bin_str)

	os.remove(outfile+".tmp")
	os.remove(outfile+".tmp.1")
	decoded_str = ''.join(decoded_data)
	decoded_str = decoded_str[:data_len]
	f_out = open(outfile,'w')
	f_out.write(binary_string_to_bytes(decoded_str))
	f_out.close()
	return 0


def simulate_indelsubs(read, sub_prob = 0.0, del_prob = 0.0, ins_prob = 0.0):
    '''
    add iid indels and substitions to read (copied from flappie_new/viterbi/util.py)
    '''
    char_list = [c for c in read]
    pos_in_char_list = 0
    new_char_list = []
    alphabet = {}
    alphabet['all'] = ['A','C','G','T']
    alphabet['A'] = ['C','G','T']
    alphabet['C'] = ['A','G','T']
    alphabet['G'] = ['C','A','T']
    alphabet['T'] = ['C','G','A']
    while True:
        ins = (np.random.random_sample()<ins_prob)
        if ins:
            new_char_list.append(np.random.choice(alphabet['all']))
        else:
            if pos_in_char_list == len(char_list):# end of original read and not inserting
                break
            _del = (np.random.random_sample()<del_prob)
            if _del:
                pos_in_char_list += 1
            else:
                sub = (np.random.random_sample()<sub_prob)
                if sub:
                    new_char_list.append(np.random.choice(alphabet[char_list[pos_in_char_list]]))
                else:
                    new_char_list.append(char_list[pos_in_char_list])
                pos_in_char_list += 1
    return ''.join(new_char_list)

def sample_reads_indel(infile,outfile,num_reads,sub_prob = 0.0, del_prob = 0.0, ins_prob = 0.0):
	'''
	Sample num_reads (with indels and substitutions) from oligos in infile and write to outfile.
	'''
	dna2int = {'A':0,'C':1,'G':2,'T':3}
	int2dna = {0:'A',1:'C',2:'G',3:'T'}
	f_in = open(infile,'r')
	f_out = open(outfile,'w')
	input_lines = f_in.readlines();
	f_in.close()
        unique_reads = set([])
	for i in range(num_reads):
		clean_sample = random.choice(input_lines).rstrip('\n')
                unique_reads.add(clean_sample)
                output_sample = simulate_indelsubs(clean_sample, sub_prob=sub_prob, del_prob=del_prob, ins_prob=ins_prob)
		f_out.write("%s\n" %output_sample)
	f_out.close()
        print('Number of unique oligos = ', len(unique_reads))
	return

def find_min_coverage(infile_data,oligo_length,BCH_bits,LDPC_alpha,LDPC_prefix,file_size, sub_prob, eps_decode, num_experiments, mode = 'correct',ins_prob=0.0,del_prob=0.0,start_coverage=1.0,sync='',sync_pos=-1):
	'''
	Find minimum coverage (in steps of 0.2) when we have 100% successes in num_experiment trials with the given parameters.
	'''
	with open(infile_data,'r') as f:
		orig_str = f.read()
	# Encode data
	outfile_oligos = infile_data +'.tmp.oligo'
	outfile_reads = infile_data+".tmp.reads"
	outfile_dec = infile_data+".tmp.dec"
	encode_data_noRLL(infile_data,oligo_length,outfile_oligos,BCH_bits,LDPC_alpha,LDPC_prefix,sync=sync,sync_pos=sync_pos)
	coverage = start_coverage
	while True:
		num_reads = int(coverage*file_size*4.0/oligo_length)
		print 'coverage:',coverage
		print 'num_reads:',num_reads

		num_successes = 0
		for _ in range(num_experiments):
			sample_reads_indel(outfile_oligos,outfile_reads,num_reads,sub_prob=sub_prob,del_prob=del_prob,ins_prob=ins_prob)
			status = decode_data_noRLL(outfile_reads,oligo_length,outfile_dec,BCH_bits,LDPC_alpha,LDPC_prefix,file_size, eps_decode, mode,MSA=True,sync=sync,sync_pos=sync_pos)
			if status == 0:
				with open(outfile_dec,'r') as f:
					dec_str = f.read()
				if dec_str == orig_str:
					num_successes += 1
                                else:
                                    break
			else:
				break
		if num_successes == num_experiments:
			break
		coverage += 0.2
	os.remove(outfile_reads)
	os.remove(outfile_dec)
	os.remove(outfile_oligos)
	return coverage

def remove_barcodes_flexbar(infile_reads, start_barcode, end_barcode, outfile_reads):
    start_barcode_RC = reverse_complement(end_barcode)
    end_barcode_RC = reverse_complement(start_barcode)
    start_barcode_list = [start_barcode,start_barcode_RC]
    end_barcode_list = [end_barcode,end_barcode_RC]
    numreads = 0
    numreads_forward = 0
    numreads_reverse = 0
    numreads_failed = 0
    fout = open(outfile_reads,'w')
    flexbarinfile = infile_reads+'.flexbar.in'
    flexbaroutfile = infile_reads+'.flexbar.out'

    # first write all reads to a fasta format
    with open(infile_reads,'r') as f, open(flexbarinfile+'.fasta','w') as f1:
        for line in f:
            numreads += 1
            f1.write('>'+str(numreads)+'\n'+line)

    del_penalty = -1 # penalty for deletion/gap (note : mismatch penalty is -1, match gain is 1)
    min_adapter_overlap = 10 # minimum overlap required between adapter and read
    min_read_length = 50 # minimum read length after adapter removal
    max_error_rate = 0.2 # max error rate (indel+substitution) for adapter
    for i in range(2):
        # i = 0: forward, i = 1: reverse
        # call flexbar with the appropriate parameters
        # start barcode
        subprocess.call(['flexbar '+
                        '-r '+flexbarinfile+'.fasta '+
                        '-as '+start_barcode_list[i]+' '+
                        '-t '+flexbaroutfile+' '+
                        '-ag '+str(del_penalty)+' '+
                        '-ao '+str(min_adapter_overlap)+' '+
                        '-m '+str(min_read_length)+' '+
                        '-at '+str(max_error_rate)+' '+
                        '-ae LEFT '+
                        '-g' # to tag reads with barcode removed
                        ],
                        shell=True)
        os.rename(flexbaroutfile+'.fasta',flexbarinfile+'.fasta')
        # end barcode
        subprocess.call(['flexbar '+
                        '-r '+flexbarinfile+'.fasta '+
                        '-as '+end_barcode_list[i]+' '+
                        '-t '+flexbaroutfile+' '+
                        '-ag '+str(del_penalty)+' '+
                        '-ao '+str(min_adapter_overlap)+' '+
                        '-m '+str(min_read_length)+' '+
                        '-at '+str(max_error_rate)+' '+
                        '-ae RIGHT '+
                        '-g' # to tag reads with barcode removed
                        ],
                        shell=True)
        os.rename(flexbaroutfile+'.fasta',flexbarinfile+'.fasta')
        # now divide reads into reads that had both adapters removed and those that didn't
        with open(flexbarinfile+'.fasta','r') as f1, open(flexbaroutfile+'.fasta','w') as f2:
            index_line = f1.readline()
            tmp_counter = 0
            while True:
                index_line = index_line.rstrip('\n')
                # get entire read
                trimmed_read = ''
                while True:
                    trimmed_read_line = f1.readline()
                    if trimmed_read_line == '':
                        next_index_line = ''
                        break
                    if trimmed_read_line[0] == '>':
                        next_index_line = trimmed_read_line
                        break
                    trimmed_read += trimmed_read_line.rstrip('\n')
                if index_line.endswith('_Flexbar_removal_cmdline_Flexbar_removal_cmdline'):
                    if i == 0:
                        numreads_forward += 1
                        fout.write(trimmed_read+'\n')
                    else:
                        numreads_reverse += 1
                        fout.write(reverse_complement(trimmed_read)+'\n')
                else:
                    tmp_counter += 1
                    f2.write('>'+str(tmp_counter)+'\n'+trimmed_read+'\n')
                index_line = next_index_line
                if index_line == '': # end of file
                    break
        os.rename(flexbaroutfile+'.fasta',flexbarinfile+'.fasta')
    fout.close()
    print('numreads',numreads)
    print('numreads_forward',numreads_forward)
    print('numreads_reverse',numreads_reverse)
    print('numreads_failed',numreads-numreads_forward-numreads_reverse)

    os.remove(flexbarinfile+'.fasta')
    os.remove(flexbaroutfile+'.log')
    return
