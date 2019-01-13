import distance
fastq_name = "/raid/nanopore/shubham/illumina_fastq_1/Pool_S1_L001_R1_001.fastq"
outfile_name = "reads_GATACTATCGAGATTACTCCAAGTC_GCTAGTCGATCCTCTGCTGCAATCG_thresh_5"
start_barcode = "GATACTATCGAGATTACTCCAAGTC"
start_barcode_len = len(start_barcode)
end_barcode = "GCTAGTCGATCCTCTGCTGCAATCG"
end_barcode_len = len(end_barcode)
hamming_dist_hist = [0]*(start_barcode_len+1)
hamming_dist_end_hist = [0]*(end_barcode_len+1)
linenum = 0
readlen = 100
end_barcode_hamming_thresh = 5
fout = open(outfile_name,'w')
with open(fastq_name) as f:
	for line in f:
		if linenum%4 == 1:
#			hamming_dist_hist[distance.hamming(start_barcode,line[:start_barcode_len])] += 1 
			if start_barcode == line[:start_barcode_len]:
#                            hamming_dist_end_hist[distance.hamming(end_barcode,line[start_barcode_len+readlen:start_barcode_len+readlen+end_barcode_len])] += 1    
                            if distance.hamming(end_barcode,line[start_barcode_len+readlen:start_barcode_len+readlen+end_barcode_len]) < end_barcode_hamming_thresh:
                                fout.write(line[start_barcode_len:start_barcode_len+readlen]+'\n')
		linenum += 1

#for i in range(start_barcode_len+1):
#	print(i,hamming_dist_hist[i])
#for i in range(end_barcode_len+1):
#	print(i,hamming_dist_end_hist[i])
