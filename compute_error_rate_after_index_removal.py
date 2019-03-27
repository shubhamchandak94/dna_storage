import distance

infile_correct = 'reads_myfile_160K_9_2'
infile_index = 'analysis_9_39000_flexbar/outfile.success.index'
infile_data = 'analysis_9_39000_flexbar/outfile.success.consensus'
outfile_error_free = 'L007406_S1_L001_R1_001.reads.shuf.50000.trimmed.new.error_free'

payload_length = 82
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
        if lev >= 1:
            print('Correct:\n'+correct_payloads[index])
            print('Read:\n'+l)
            print('Diff:\n')
            print(''.join([str(int(l[i] == correct_payloads[index][i])) for i in range(payload_length)]))
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
