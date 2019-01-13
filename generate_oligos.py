infile_reads = "reads_myfile_160K_9_2"
start_barcode = "TCCTGTGCTGCCTGTAATGAGCCAA"
end_barcode = "AGCATAGAACTGAGACCACGGATTG"
outfile_oligos = "oligos_myfile_160K_9_2.fa"
prefix = "20190112_myfile_160K_9_2"

fout = open(outfile_oligos,'w')
with open(infile_reads,'r') as f:
    i = 0
    for line in f:
        read = line.rstrip('\n')
        fout.write('>'+prefix+'_'+start_barcode+'_'+end_barcode+'_'+str(i)+'\n'+start_barcode+read+end_barcode+'\n')
        i += 1
fout.close()
