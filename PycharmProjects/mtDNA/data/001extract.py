from Bio import Entrez, SeqIO
Entrez.email = 'julie.trolle@nyumc.org'

import pandas as pd
import xlrd

df = pd.read_excel('MtAccNo.xlsx', index_col=0, header=1)
Mito_AccNo = df['Mito']
#print(Mito_AccNo)

list_of_Mito_AccNo = []
for id in Mito_AccNo:
    list_of_Mito_AccNo.append(id)
print(len(list_of_Mito_AccNo))

# strain_dict = {}
# for AccNo in list_of_Mito_AccNo:
#     handle = Entrez.efetch(db='nucleotide', id=AccNo, rettype='gb', retmode='text')
#     record = SeqIO.read(handle, 'genbank')
#

for AccNo in list_of_Mito_AccNo:
    out_handle = open('%s.gb' %AccNo, 'w')
    net_handle = Entrez.efetch(db='nucleotide', id=AccNo, rettype='gb', retmode='text')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()