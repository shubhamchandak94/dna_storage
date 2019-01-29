import distance

infile_correct = 'ill_code_3_33_2'
infile_index = 'coded_15.read.1.shuf.200000.trimmed.tmp.index'
infile_data = 'coded_15.read.1.shuf.200000.trimmed.tmp.data'
outfile_error_free = 'coded_15.read.1.shuf.200000.trimmed.error_free'

payload_length = 74
correct_payloads = []
correct_index = []
with open(infile_correct) as f:
    for line in f:
        l = line.rstrip('\n')
        correct_payloads.append(l[-payload_length:])
        correct_index.append(l[:-payload_length])

fout = open(outfile_error_free,'w')
hamming_distances = []
levenshtein_distances = []
hamming_distances_per_index = [[] for _ in range(len(correct_payloads))]
print('hamming','levenshtein')
with open(infile_index) as f_index, open(infile_data) as f_data:
    for line_index, line_data in zip(f_index,f_data):
        index = int(line_index)
        l = line_data.rstrip('\n')
        if len(l) != payload_length:
            continue
        ham = distance.hamming(l,correct_payloads[index])
        lev = distance.levenshtein(l,correct_payloads[index])
        hamming_distances.append(ham)
        hamming_distances_per_index[index].append(ham)
        levenshtein_distances.append(lev)
        print(ham,lev)
        fout.write(correct_index[index]+correct_payloads[index]+'\n')
print('average hamming (%)',sum(hamming_distances)/(len(hamming_distances)*payload_length)*100)
print('average levenshtein (%)',sum(levenshtein_distances)/(len(levenshtein_distances)*payload_length)*100)
print('index, average hamming (%))')
for i in range(len(correct_payloads)):
    if len(hamming_distances_per_index[i]) != 0:
        print(i,sum(hamming_distances_per_index[i])/(len(hamming_distances_per_index[i])*payload_length)*100)
