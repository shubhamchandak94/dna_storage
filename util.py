import itertools
import math 
import numpy as np
import bchlib
import binascii
import random
import os
import subprocess

#some constants
LDPC_dim = 256000
rll2_block_len = 87
index_block_len = 16 
BCH_POLYNOMIAL = 67
BCH_bits_per_error = 6
LDPC_max_iter = 100

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

def bin2dna_rll2_low_rate(bin_string):
	'''
	Generate dna accordiing to 1.67 bpb scheme, mapping 5 bits to 3 nt.
	Generated dna string has no runs of 3 nt.
	One substitution in dna can cause at most 2 bits.
	Throw error if binary string length is not multiple of 5.	
	'''
	if len(bin_string)%5 != 0:
		raise Exception('binary string length not multiple of 5')
	# Generate dictionary
	d1 = {'00':'A','01':'C','10':'G','11':'T'}
	option_1 = {'000':'AA','001':'AC','010':'AG','011':'AT',\
		    '100':'CG','101':'TC','110':'TG','111':'GG'}
	option_2 = {'000':'CC','001':'GT','010':'CT','011':'GC',\
		    '100':'TA','101':'GA','110':'CA','111':'TT'}
	d = {}
	for b in itertools.product('01',repeat=5):
		b = ''.join(b)
		mid_base = d1[b[0:2]]
		s1 = option_1[b[2:]]
		s2 = option_2[b[2:]]
		if s1[0]!=mid_base and s1[1]!=mid_base:
			d[b] = s1[0]+mid_base+s1[1]
		else:
			d[b] = s2[0]+mid_base+s2[1]
	return ''.join([d[bin_string[5*i:5*(i+1)]] for i in range(len(bin_string)/5)])

def dna2bin_rll2_low_rate(dna_string):
	'''
	Decode dna data encoded according to bin2dna_rll2_low_rate and possibly passed through noisy channel.
	Throw error if dna string length is not multiple of 3.
	'''
	if len(dna_string)%3 != 0:
		raise Exception('dna string length not multiple of 3')
	# Generate dictionary
	d = {}
	d1 = {'A':'00','C':'01','G':'10','T':'11'}
	d2 = {'AA':'000','AC':'001','AG':'010','AT':'011',\
	      'CG':'100','TC':'101','TG':'110','GG':'111',	
    	      'CC':'000','GT':'001','CT':'010','GC':'011',\
	      'TA':'100','GA':'101','CA':'110','TT':'111'}
	for s in itertools.product('ACGT',repeat=3):
		s = ''.join(s)
		d[s] = d1[s[1]]+d2[s[0]+s[2]]
	return ''.join([d[dna_string[3*i:3*(i+1)]] for i in range(len(dna_string)/3)])

def generate_rll2_count_arrays(block_len):
	'''
	Generate arrays used for high rate encoding to dna without runs of 3.
	block_len is the number of dna symbols
	returns (bin_block_len,rll2_count)
	where bin_block_len is the number of binary bits converted to block_len dna symbols, rll2_count is array used for encoding and decoding.
	'''
	rll2_count = [3,12]
	for i in range(block_len-2):
		rll2_count.append(3*rll2_count[-1]+3*rll2_count[-2])
	bin_block_len = int(math.floor(math.log(rll2_count[-1],2)))
	# rll_count stores number of 0,1,2,3 strings of that don't start with a 3 and don't contain 33 as a substring
	return (bin_block_len,[1]+rll2_count)


def bin2dna_rll2_high_rate(bin_string, block_len, bin_block_len, rll2_count):
	'''
	Encode binary data to DNA using a high rate encoding, writing bin_block_len binary symbols to block_len dna symbols without any runs of 3.
	The last two parameters of this funciton can be obtained by callling generate_rll2_count_arrays(block_len).
	Throw error if binary string length is not multiple of the binary block length as returned by generate_rll2_count_arrays.
	'''
	if len(bin_string)%bin_block_len != 0:
		raise Exception('binary string length not multiple of bin_block length')
	delta_coded = np.zeros(len(bin_string)/bin_block_len*block_len,dtype=int)
	for i in range(len(bin_string)/bin_block_len):
		num = int(bin_string[i*bin_block_len:(i+1)*bin_block_len],2)
		for j in range(block_len):
			next_symbol = num/rll2_count[block_len-j-1]
			print next_symbol
			if next_symbol > 3:
				raise Exception('next_symbol>3')
			num -= next_symbol*rll2_count[block_len-j-1]
			delta_coded[i*block_len+j] = next_symbol
		if num != 0:
			raise Exception('num not zero')
	#swaps 0s with 3s for simpler delta coding
	ind_0 = np.nonzero(delta_coded == 0)
	ind_3 = np.nonzero(delta_coded == 3)
	delta_coded[ind_0] = 3
	delta_coded[ind_3] = 0
	quaternary_coded = np.mod(np.cumsum(delta_coded),4)
	d = {0:'A',1:'C',2:'G',3:'T'}
	return ''.join([d[j] for j in quaternary_coded])	

# TODO fix error
def dna2bin_rll2_high_rate(dna_string, block_len, bin_block_len, rll2_count):
	'''
	Decoder for bin2dna_rll2_high_rate.
	Return (status,bin_string) where status is 0 if success and 1 if any of the blocks fail to decode (we use delta coding internally so it's important that the entire string is correct).
	Throw error if dna string length is not multiple of the block length.
	'''
	if len(dna_string)%block_len != 0:
		raise Exception('dna string length not multiple of block length')
	d = {'A':0,'C':1,'G':2,'T':3}
	quaternary_coded = np.array([d[s] for s in dna_string],dtype=int)
	delta_coded = np.mod(np.diff(quaternary_coded),4)
	delta_coded = np.insert(delta_coded,0,quaternary_coded[0])
	#swaps 0s with 3s 
	ind_0 = np.nonzero(delta_coded == 0)
	ind_3 = np.nonzero(delta_coded == 3)
	delta_coded[ind_0] = 3
	delta_coded[ind_3] = 0
	bin_string_array = []
	delta_coded = delta_coded.tolist()# for inf precision arithmetic
	print delta_coded
	for i in range(len(dna_string)/block_len):
		# check if valid 
		if delta_coded[i*block_len] == 3 or sum([(delta_coded[i*block_len+j]+delta_coded[i*block_len+j+1]==6) for j in range(block_len-1)]) != 0:
			print "hello1"
			return (1,'')
		num = sum([delta_coded[i*block_len+j]*rll2_count[block_len-j-1] for j in range(block_len)])
		if num >= rll2_count[block_len]:
			#Error
			print "hello2"
			return (1,'')
		bin_string_array.append((bin(num)[2:]).zfill(bin_block_len))
	return (0,''.join([s for s in bin_string_array]))
	
def bin2dna_rll2_index(bin_string, block_len, bin_block_len, rll2_count):
	'''
	Similar to bin2dna_rll1_high_rate but operates on a single index block.
	Guarantees that the last 2 bases are different and no runs of 3 nt.
	Throw error if binary string length is not equal to binary block length
	'''
	if len(bin_string) != bin_block_len:
		raise Exception('binary string length not equal to bin_block length')
	delta_coded = np.zeros(block_len,dtype=int)
	num = int(bin_string,2)
	for j in range(block_len):
		next_symbol = num/rll2_count[block_len-j-1]
		if next_symbol > 3:
			raise Exception('next_symbol>3')
		num -= next_symbol*rll2_count[block_len-j-1]
		delta_coded[j] = next_symbol
	if num != 0:
		raise Exception('num not zero')
	#swaps 0s with 3s for simpler delta coding
	ind_0 = np.nonzero(delta_coded == 0)
	ind_3 = np.nonzero(delta_coded == 3)
	delta_coded[ind_0] = 3
	delta_coded[ind_3] = 0
	delta_coded = np.flipud(delta_coded) # so last two bases are unequal (since first one never 0)
	quaternary_coded = np.mod(np.cumsum(delta_coded),4)
	d = {0:'A',1:'C',2:'G',3:'T'}
	return ''.join([d[j] for j in quaternary_coded])	
	
def dna2bin_rll2_index(dna_string, block_len, bin_block_len, rll2_count):
	'''
	Decoder for bin2dna_rll2_index.
	Return (status,bin_string) where status is 0 if success and 1 if failure.
	Throw error if dna string length is not equal to the block length.
	'''
	if len(dna_string) != block_len:
		raise Exception('dna string length not equal to block length')
	d = {'A':0,'C':1,'G':2,'T':3}
	quaternary_coded = np.array([d[s] for s in dna_string],dtype=int)
	delta_coded = np.mod(np.diff(quaternary_coded),4)
	delta_coded = np.insert(delta_coded,0,quaternary_coded[0])
	delta_coded = np.flipud(delta_coded)
	#swaps 0s with 3s 
	ind_0 = np.nonzero(delta_coded == 0)
	ind_3 = np.nonzero(delta_coded == 3)
	delta_coded[ind_0] = 3
	delta_coded[ind_3] = 0
	# check if valid 
	if delta_coded[0] == 3 or sum([(delta_coded[i]+delta_coded[i+1]==6) for i in range(block_len-1)]) != 0:
		return (1,'')
	bin_string_array = []
	delta_coded = delta_coded.tolist()# for inf precision arithmetic
	num = sum([delta_coded[j]*rll2_count[block_len-j-1] for j in range(block_len)])
	if num >= rll2_count[block_len]:
		#Error
		return (1,'')
	bin_string_array = (bin(num)[2:]).zfill(bin_block_len)
	return (0,bin_string_array)
	
def binary_string_to_bytes(bin_string):
	'''
	Convert binary string to bytes (needed for BCH).
	Throw error if bin_string length not multiple of 8.
	'''
	if len(bin_string)%8 != 0:
		raise Exception('binary string length not multiple of 8')
	return binascii.unhexlify(((hex(int(bin_string,2)))[2:]).zfill(len(bin_string)/4))

def bytes_to_binary_string(byte_string):
	'''
	Convert bytes to binary string (needed for BCH).
	'''
	return bin(int(binascii.hexlify(byte_string), 16))[2:].zfill(len(byte_string)*8)

def add_index(num_oligos,BCH_bits,infile_name,outfile_name):
	'''
	Generate DNA encoding of all indexes from 0 to num_oligos-1, each protected with BCH_bits protection and RLL encoded. 
	The index DNA strings are appropriately concatenated with corresponding lines from infile_name and the resulting oligos are written to outfile_name, one line per index. 
	Throw error if num_oligos > 2**30 = 1073741824.
	'''
	if num_oligos > 2**30:
		raise Exception('Too many oligos')
	block_len = index_block_len
	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/5))*3
		num_bases_index = block_len + num_bases_BCH
	# dictionary to generate separator base b/w index and payload
	d = {'AA':'C','AC':'G','AG':'T','AT':'G',\
	      'CG':'T','TC':'A','TG':'C','GG':'A',	
    	      'CC':'A','GT':'C','CT':'A','GC':'T',\
	      'TA':'G','GA':'T','CA':'G','TT':'C'}
	index = 0
	bin_block_len, rll2_count = generate_rll2_count_arrays(block_len)
	with open(infile_name) as infile, open(outfile_name, 'w') as outfile:
		for line in infile:
			bin_string = bin(index)[2:].zfill(bin_block_len)
			dna = bin2dna_rll2_index(bin_string, block_len, bin_block_len, rll2_count)
			if BCH_bits != 0:
				bits_ecc = bytes_to_binary_string(bch.encode(binary_string_to_bytes(dna2bin_2bpb(dna))))
				bits_ecc = bits_ecc[:BCH_bits*BCH_bits_per_error]
				dna_ecc = bin2dna_rll2_low_rate(bits_ecc.zfill(num_bases_BCH*5/3))
				sep_base = d[dna_ecc[-1]+line[0]]
				outfile.write(dna+dna_ecc+sep_base+line)
			else:
				sep_base = d[dna[-1]+line[0]]
				outfile.write(dna+sep_base+line)
			index += 1
			if index == num_oligos:
				break
			
def remove_index(num_oligos,BCH_bits,infile_name,outfile_data,outfile_index,mode="correct"):
	'''
	Decode index from a collection of (noisy) reads in infile_name and write data and index to outfile_data and outfile_index, line by line, skipping positions where index failed to decode.
	Mode can be "correct" or "detect":
	Correct - erase if correction fails
	Detect - don't try to correct, erase if any error detected
	'''
	if num_oligos > 2**30:
		raise Exception('Too many oligos')
	block_len = index_block_len
	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/5))*3
		num_bases_index = block_len + num_bases_BCH
	else:
		num_bases_index = block_len	
	bin_block_len, rll2_count = generate_rll2_count_arrays(block_len)
	count = 0
	with open(infile_name) as infile, open(outfile_data, 'w') as f_data, open(outfile_index, 'w') as f_index:
		for line in infile:
			dna_data = line[num_bases_index+1:]
			dna_index = line[:num_bases_index]		
			if BCH_bits != 0:
				dna_ecc = dna_index[-num_bases_BCH:]
				bin_ecc = dna2bin_rll2_low_rate(dna_ecc)
				bin_ecc = bin_ecc[-BCH_bits*BCH_bits_per_error:]
				extra_len = 8*int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/8))-len(bin_ecc)
				bin_ecc = bin_ecc+'0'*extra_len
				(bitflips,cor_index,cor_ecc) = bch.decode(binary_string_to_bytes(dna2bin_2bpb(dna_index[:-num_bases_BCH])),binary_string_to_bytes(bin_ecc))
				if bitflips >= 0: #success in correction
					if mode == "detect" and bitflips > 0:
						continue
					(status,bin_index) = dna2bin_rll2_index(bin2dna_2bpb(bytes_to_binary_string(cor_index)),block_len, bin_block_len, rll2_count)
					if status == 0:
						if int(bin_index,2) < num_oligos:
							f_data.write(dna_data)
							f_index.write(str(int(bin_index,2))+'\n')
							count += 1
			else:
				(status,bin_index) = dna2bin_rll2_index(dna_index,block_len, bin_block_len, rll2_count)
				if status == 0:
					if int(bin_index,2) < num_oligos:
						f_data.write(dna_data)
						f_index.write(str(int(bin_index,2))+'\n')	
						count += 1
	print "Successfully decoded",count,"indices"		

def rll2_pad(dna_str,padded_len):
	'''
	Pad dna_str to padded_len such that the padding does not introduce any runs of >= 3
	'''	
	extra_len = padded_len - len(dna_str)
	if extra_len == 0:
		return dna_str
	d = {'A':'C','C':'G','G':'T','T':'A'}
	d1 = {'A':0,'C':1,'G':2,'T':3}
	d2 = {0:'A',1:'C',2:'G',3:'T'}
	delta_coded = np.random.choice([1,2,3],extra_len)
	delta_coded[0] = d1[d[dna_str[-1]]]
	quaternary_coded = np.mod(np.cumsum(delta_coded),4)
	return dna_str+''.join([d2[i] for i in quaternary_coded])		
	
def encode_data(infile,oligo_length,outfile,BCH_bits,LDPC_alpha,LDPC_prefix):
	'''
	Encode binary data in infile to oligos written to oufile.
	LDPC_prefix.pchk and LDPC_prefix.gen are the LDPC matrices. LDPC_prefix.systematic contains positions of systematic bits.
	int(LDPC_dim*LDPC_alpha) should be the number of parity check bits.
	infile is a bytestream file.
	'''
	f_in = open(infile,'r')
	data = f_in.read()	
	bin_data = bytes_to_binary_string(data)

	# calculate various parameters for encoding
	data_len = len(bin_data)
	block_len = rll2_block_len
	bin_block_len, rll2_count = generate_rll2_count_arrays(block_len)
	coded_blocks_per_LDPC_block = LDPC_dim/(2*block_len)
	uncoded_bits_per_LDPC_block = coded_blocks_per_LDPC_block*bin_block_len
	num_LDPC_blocks = int(math.ceil(1.0*data_len/uncoded_bits_per_LDPC_block))
	parity_bits_per_LDPC_block = int(LDPC_alpha*LDPC_dim)
	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/5))*3
		num_bases_index = index_block_len + num_bases_BCH
	else:
		num_bases_index = index_block_len
	num_bases_payload = oligo_length - num_bases_index - 1
	parity_bits_per_oligo = (num_bases_payload/3)*5
	num_oligos_data_per_LDPC_block = int(math.ceil(1.0*LDPC_dim/(2*num_bases_payload)))
	num_oligos_parity_per_LDPC_block = int(math.ceil(1.0*parity_bits_per_LDPC_block/parity_bits_per_oligo))
	num_oligos_per_LDPC_block = num_oligos_data_per_LDPC_block+num_oligos_parity_per_LDPC_block
	overall_rate = 1.0*data_len/(num_oligos_per_LDPC_block*num_LDPC_blocks*oligo_length)
	print 'overall rate:' ,overall_rate, 'bpb'
	print 'num oligos:', num_LDPC_blocks*num_oligos_per_LDPC_block
	print 'oligo length:', oligo_length,'bases'
	print 'bases per oligo for index + index parity + 1:', num_bases_index + 1
	print 'fraction of oligos used for parity check:', 1.0*num_oligos_parity_per_LDPC_block/num_oligos_per_LDPC_block
	print 'number of LDPC blocks:', num_LDPC_blocks
	f_out = open(outfile+'.tmp','w')
	# pad bin_data to multiple of uncoded_bits_per_LDPC_block
	extra_len = uncoded_bits_per_LDPC_block*num_LDPC_blocks - len(bin_data)
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

	for i in range(num_LDPC_blocks):
		bin_data_block = bin_data[i*uncoded_bits_per_LDPC_block:(i+1)*uncoded_bits_per_LDPC_block]
		# encode RLL high rate
		rll2_dna_data_block = bin2dna_rll2_high_rate(bin_data_block, block_len, bin_block_len, rll2_count)
		f_LDPC_input = open(outfile+'.tmp.1','w')
		f_LDPC_input.write(dna2bin_2bpb(rll2_pad(rll2_dna_data_block,LDPC_dim/2)))
		f_LDPC_input.close()
		# encode LDPC
		subprocess.call(["./LDPC-codes/encode "+LDPC_prefix+".pchk "+LDPC_prefix+".gen "+outfile+'.tmp.1 '+outfile+'.tmp.2 '], shell=True)
		# read parity bits 
		f_LDPC_output = open(outfile+'.tmp.2','r')
		encoded_LDPC = f_LDPC_output.read().rstrip('\n')
		f_LDPC_output.close()
		parity_bits_block = ''.join([encoded_LDPC[j] for j in parity_pos])
		# pad to multiple of parity_bits_per_oligo
		extra_len = parity_bits_per_oligo*num_oligos_parity_per_LDPC_block - len(parity_bits_block)
		parity_bits_block = parity_bits_block + ''.join(np.random.choice(['0','1'],extra_len))
		# Write dna oligos (without index)	
		# First write for data
		for j in range(num_oligos_data_per_LDPC_block-1):
			f_out.write(rll2_dna_data_block[j*num_bases_payload:(j+1)*num_bases_payload]+'\n')
		# for last oligo, might need to pad
		f_out.write(rll2_pad(rll2_dna_data_block[(num_oligos_data_per_LDPC_block-1)*num_bases_payload:],num_bases_payload)+'\n')
		# Now write parity bits
		for j in range(num_oligos_parity_per_LDPC_block):
			f_out.write(rll2_pad(bin2dna_rll2_low_rate(parity_bits_block[j*parity_bits_per_oligo:(j+1)*parity_bits_per_oligo]),num_bases_payload)+'\n')

	rll2_dna_data = bin2dna_rll2_high_rate(bin_data, block_len, bin_block_len, rll2_count)
	f_out.close()
	add_index(num_oligos_per_LDPC_block*num_LDPC_blocks,BCH_bits,outfile+'.tmp',outfile)
	os.remove(outfile+'.tmp')
	os.remove(outfile+'.tmp.1')
	os.remove(outfile+'.tmp.2')
	return

def decode_data(infile,oligo_length,outfile,BCH_bits,LDPC_alpha,LDPC_prefix,file_size, eps, mode = 'correct'):
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
	# calculate various parameters for decoding
	block_len = rll2_block_len
	bin_block_len, rll2_count = generate_rll2_count_arrays(block_len)
	coded_blocks_per_LDPC_block = LDPC_dim/(2*block_len)
	uncoded_bits_per_LDPC_block = coded_blocks_per_LDPC_block*bin_block_len
	num_LDPC_blocks = int(math.ceil(1.0*data_len/uncoded_bits_per_LDPC_block))
	parity_bits_per_LDPC_block = int(LDPC_alpha*LDPC_dim)
	if BCH_bits != 0:
		bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
		# calculate number of bases used for index
		num_bases_BCH = int(math.ceil(1.0*BCH_bits*BCH_bits_per_error/5))*3
		num_bases_index = index_block_len + num_bases_BCH
	else:
		num_bases_index = index_block_len
	num_bases_payload = oligo_length - num_bases_index - 1
	parity_bits_per_oligo = (num_bases_payload/3)*5
	num_oligos_data_per_LDPC_block = int(math.ceil(1.0*LDPC_dim/(2*num_bases_payload)))
	num_oligos_parity_per_LDPC_block = int(math.ceil(1.0*parity_bits_per_LDPC_block/parity_bits_per_oligo))
	num_oligos_per_LDPC_block = num_oligos_data_per_LDPC_block+num_oligos_parity_per_LDPC_block
	overall_rate = 1.0*data_len/(num_oligos_per_LDPC_block*num_LDPC_blocks*oligo_length)
	print 'overall rate:' ,overall_rate, 'bpb'
	print 'num oligos:', num_LDPC_blocks*num_oligos_per_LDPC_block
	print 'oligo length:', oligo_length,'bases'
	print 'bases per oligo for index + index parity + 1:', num_bases_index + 1
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
	
	# first decode index	
	remove_index(num_LDPC_blocks*num_oligos_per_LDPC_block,BCH_bits,infile,infile+'.tmp.data',infile+'.tmp.index')
	
	# Now store counts in an array
	total_counts = np.zeros((num_LDPC_blocks,LDPC_dim+parity_bits_per_LDPC_block))
	zero_counts = np.zeros((num_LDPC_blocks,LDPC_dim+parity_bits_per_LDPC_block))
	with open(infile+'.tmp.data','r') as infile_data, open(infile+'.tmp.index','r') as infile_index:
		for line_data, line_index in zip(infile_data,infile_index):
			index = int(line_index)
			block_number = index/num_oligos_per_LDPC_block
			index_in_block = index%num_oligos_per_LDPC_block
			if index_in_block < num_oligos_data_per_LDPC_block:
				# data oligo
				start_pos = index_in_block*num_bases_payload*2
				end_pos = (index_in_block+1)*num_bases_payload*2
				if end_pos > LDPC_dim:
					end_pos = LDPC_dim
				total_counts[block_number][sys_pos[start_pos:end_pos]] += 1
				bin_str = dna2bin_2bpb(line_data[:(end_pos-start_pos)/2])
				bin_arr = np.array([int(c) for c in bin_str])
				zero_counts[block_number][sys_pos[start_pos:end_pos]] += 1-bin_arr
			else:
				# parity oligo
				start_pos = parity_bits_per_oligo*(index_in_block-num_oligos_data_per_LDPC_block)
				end_pos = start_pos + parity_bits_per_oligo
				if end_pos > parity_bits_per_LDPC_block:
					end_pos = parity_bits_per_LDPC_block
				bin_str = dna2bin_rll2_low_rate(line_data[:3*(num_bases_payload/3)])
				bin_arr = np.array([int(c) for c in bin_str])
				total_counts[block_number][parity_pos[start_pos:end_pos]] += 1
				zero_counts[block_number][parity_pos[start_pos:end_pos]] += 1-bin_arr[:end_pos-start_pos]
	os.remove(infile+'.tmp.data')	
	os.remove(infile+'.tmp.index')
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
		LDPC_decoded_str = ''.join([LDPC_decoded_str[j] for j in sys_pos])
		LDPC_decoded_str = LDPC_decoded_str[:2*coded_blocks_per_LDPC_block*block_len]
		(status,LDPC_decoded_bin_str) = dna2bin_rll2_high_rate(bin2dna_2bpb(LDPC_decoded_str),block_len, bin_block_len, rll2_count)
		if status != 0:
			# failure
			os.remove(outfile+".tmp")
			os.remove(outfile+".tmp.1")
			return 1	
		decoded_data.append(LDPC_decoded_bin_str)

	os.remove(outfile+".tmp")
	os.remove(outfile+".tmp.1")
	decoded_str = ''.join(decoded_data)
	decoded_str = decoded_str[:data_len] 
	f_out = open(outfile,'w')
	f_out.write(binary_string_to_bytes(decoded_str))
	f_out.close()
	return 0
