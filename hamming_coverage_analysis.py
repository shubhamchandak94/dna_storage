import distance

read_file = 'tmp.reads'#'reads_TGAACCGTATGCGACTCATACCGCC_ATCGGCGCGCCATGATGGAGCGAGA_thresh_5_8000'
code_file = 'ill_code_3_33_2'

f = open(code_file)
code_list = [l.rstrip('\n') for l in f.readlines()]
f.close()
oligos_seen = set([])
hamming = 0
with open(read_file) as f:
    for line in f:
        l = line.rstrip('\n')
        hamming_list = [distance.hamming(l,l_1) for l_1 in code_list]
        min_hamming = min(hamming_list)
        hamming += min_hamming
        oligos_seen.add(hamming_list.index(min_hamming))
print('Hamming:',hamming)
print('Oligos seen:',len(oligos_seen))
