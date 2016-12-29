from Bio import SeqIO
from collections import OrderedDict
import os
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import IUPAC

g = open('/Users/julietrolle/PycharmProjects/mtDNA/data/final_mito_len.gb', 'r')
han = SeqIO.read(g, 'genbank')

intergenic_features = []
longest_intergenic_features = []

for (i, f) in enumerate(han.features):
    record = han.features[i]
    if record.type == 'intergenic':
        intergenic_features.append(record)
        if len(record) > 3900:
            print(record.qualifiers['name'], len(record))
            longest_intergenic_features.append(record)
    else:
        continue

print(longest_intergenic_features)

#print length of these for all gb files, send to Suds

list_dir_files = os.listdir('/Users/julietrolle/PycharmProjects/mtDNA/data/genbank_files')

list_gb_files = [] #to ensure that only gb files are opened
for element in list_dir_files:
    if element.endswith('gb') == True:
        list_gb_files.append(element)

len_indv_int_feat = {}

for gb in list_gb_files:
    f = open('/Users/julietrolle/PycharmProjects/mtDNA/data/genbank_files/%s' % gb, 'r')
    handle = SeqIO.read(f, 'genbank')

    # for (i, f) in enumerate(handle.features):
    #     rec = handle.features[i]
    #     for rec in longest_intergenic_features:
    #         len_indv_int_feat[rec.qualifiers['name']][gb] = len(rec)

print(len_indv_int_feat)
#read up on alignment, N/C terminus tagging etc, mt transport

