import distance
fastq_name = "Pool_S1_L001_R1_001.fastq"
outfile_name = "reads_GATACTATCGAGATTACTCCAAGTC"
start_barcode = "GATACTATCGAGATTACTCCAAGTC"
start_barcode_len = len(start_barcode)
hamming_dist_hist = [0]*(start_barcode_len+1)
linenum = 0
readlen = 100
fout = open(outfile_name,'w')
with open(fastq_name) as f:
	for line in f:
		if linenum%4 == 1:
#			hamming_dist_hist[distance.hamming(start_barcode,line[:start_barcode_len])] += 1 
			if start_barcode == line[:start_barcode_len]:
				fout.write(line[start_barcode_len:start_barcode_len+readlen]+'\n')
		linenum += 1

#for i in range(start_barcode_len+1):
#	print(i,hamming_dist_hist[i])
