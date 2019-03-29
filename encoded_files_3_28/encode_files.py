import sys
sys.path.append('/raid/nanopore/shubham/dna_storage')
import util_final as util
import filecmp
import os
num_files = 7

start_barcodes = [
'GCTACATGTATACTGCGAGACAGAC',
'TCTATCTACTCGTGCTCGCTAGCTG',
'CATCAGCAGTAGAGAGTAGCGCGAT',
'TCACGAGATAGTCACGCTGTCACGT',
'GACGCGTCTAGCTGCATCTCGTACA',
'AGCGATGCTCTCTCAGTAGCATCTA',
'TAGTGTCGTCTGAGCGCAGAGATAT',
]

end_barcodes = [
'CGATAGTCGCAGTCGCACATCACTC',
'TGAGATCACAGCTACATAGTGAGAG',
'GTGTCACTATATCGCTCTACTGTGA',
'AGCTCGATGCGTGAGTGACTGTGAG',
'GTGTAGTGCGCATATGTCTAGACGT',
'ACACTACGAGCTATAGATCTACTGC',
'TCGCGCTAGACTCTGTATGAGTCGC',
]

infile_name = [
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',
'288.tar.gz.encrypted',
]

LDPC_code_prefix = [
'../LDPC-codes/matrices/ldpc_9_new',
'../LDPC-codes/matrices/ldpc_13_new',
'../LDPC-codes/matrices/ldpc_33_new',
'../LDPC-codes/matrices/ldpc_9_new',
'../LDPC-codes/matrices/ldpc_9_new',
'../LDPC-codes/matrices/ldpc_9_new',
'../LDPC-codes/matrices/ldpc_33_new',
]

LDPC_alpha = {
'../LDPC-codes/matrices/ldpc_9_new': 0.5,
'../LDPC-codes/matrices/ldpc_13_new': 0.3,
'../LDPC-codes/matrices/ldpc_33_new': 0.1,
}
BCH_bits = [
2,
2,
2,
3,
1,
2,
1,
]

sync = [
'AGT',
'AGT',
'AGT',
'AGT',
'AGT',
'',
'',
]

bin_index_len = 14
oligo_len = 100

for i in range(num_files):
    index_len = (bin_index_len+6*BCH_bits[i])//2
    sync_pos = (oligo_len-index_len)//2+index_len
    print('sync_pos',sync_pos)
    file_size = os.path.getsize(infile_name[i])
    util.encode_data(infile_name[i],oligo_len,'reads.'+str(i),BCH_bits[i],LDPC_alpha[LDPC_code_prefix[i]],LDPC_code_prefix[i],bin_index_len,sync = sync[i], sync_pos = sync_pos)
    # test decoding to see that there is no strange issue
    tmpfile_decoded = 'tmpfile_decoded'
    util.decode_data('reads.'+str(i),oligo_len,tmpfile_decoded,bin_index_len,BCH_bits[i],LDPC_alpha[LDPC_code_prefix[i]],LDPC_code_prefix[i],file_size,0.01,sync=sync[i],sync_pos=sync_pos)
    assert filecmp.cmp(infile_name[i],tmpfile_decoded)
    os.remove(tmpfile_decoded)
    with open('reads.'+str(i)) as f_reads, open('oligos.'+str(i),'w') as f_oligos:
        for j, line in enumerate(f_reads):
            f_oligos.write('>oligos_'+str(i)+'_'+start_barcodes[i]+'_'+end_barcodes[i]+'_'+str(j)+'\n')
            f_oligos.write(start_barcodes[i]+line.rstrip('\n')+end_barcodes[i]+'\n')
