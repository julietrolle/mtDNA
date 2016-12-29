from Bio import SeqIO
from collections import OrderedDict
import os
from Bio.SeqUtils import GC

class mtDNA_SequenceFeature(object):
    'sequence features for mtDNA project CDS,tRNA, rRNA and ncRNA genes'
    def __init__(self, name, seq, location, locus_tag, db_xref, genome_name, whole_genome_record):
        self.name = name
        self.seq = seq
        self.locus_tag = locus_tag
        self.location = location
        self.db_xref = db_xref
        self.genome_name = genome_name
        self.whole_genome_record = whole_genome_record
    def __repr__(self):
        return (str(self.name) + "-" + str(self.genome_name))


class mtDNA_IntergenicFeature(object):
    'sequence features for mtDNA project for intergenic spaces'
    def __init__(self, name, seq, location_start, location_end, genome_name, whole_genome_record):
        self.name = name
        self.seq = seq
        self.location_start = location_start
        self.location_end = location_end
        self.genome_name = genome_name
        self.whole_genome_record = whole_genome_record
    def __repr__(self):
        return (str(self.name) + "-" + str(self.genome_name))

all_gene_features_dict = OrderedDict()
all_intergenic_features_dict = OrderedDict()
all_no_features_dict = OrderedDict()
parse = open('parsing_stats.txt', 'w')
problematic = open('problematic_files.txt', 'w')
count = 0

f = open('/Users/julietrolle/PycharmProjects/mtDNA/data/genbank_files/CP004115.gb')
genome = SeqIO.read(f, 'genbank')

list_genes = [['COX1'], ['ATP8'], ['ATP6'], ['COB'], ['OLI1'], ['VAR1'], ['COX3'], ['COX2']]
tRNA_list = [['tRNA-Pro']]
CDS_count = 0
tRNA_count = 0
rRNA_count = 0
ncRNA_count = 0
everything_else_count = 0
no1 = []
no2 = []
no_final = []

gene_features = OrderedDict()

sep = '.'
genome_name = 'CP004115'  #gb.split(sep, 1)[0] #splits gb into two elements and extracts the first element - just gives the genome a name
list_intergenic_locations = []
list_intergenic_name = []

for (i, f) in enumerate(genome.features):
    record = genome.features[i]
    if record.type == 'CDS':
        try: #if no try statement, a KeyError is raised as COX2 has no qualifier labelled gene
            if record.qualifiers['gene'] in list_genes:
                gene_name = record.qualifiers['gene'][0]
                feature = mtDNA_SequenceFeature(gene_name, record.extract(genome.seq), record.location,
                                                record.qualifiers['locus_tag'][0],
                                                record.qualifiers['db_xref'][0], genome_name, genome)
                print(gene_name)
                gene_features[gene_name] = feature #appends gene_name-genome.name to ordered dict

                list_intergenic_locations.append(record.location)
                list_intergenic_name.append(gene_name)
                CDS_count = CDS_count + 1
        except KeyError:
            #print('passed over a KeyError', gene_name)
            continue

    elif record.type == 'tRNA':
        seq = record.extract(genome.seq)
        gene_name = record.qualifiers['product'][0]
        print(gene_name)
        if record.qualifiers['product'][0] in tRNA_list: #adds a '2' to the gene_name if it has already appeared
            gene_name = str(record.qualifiers['product'][0] + '2')
            print(gene_name)
            feature = mtDNA_SequenceFeature(gene_name, record.extract(genome.seq), record.location,
                                            record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                            genome_name, genome)
            gene_features[gene_name] = feature
            list_intergenic_locations.append(record.location)
            list_intergenic_name.append(gene_name)
            tRNA_count = tRNA_count + 1
            continue

        feature = mtDNA_SequenceFeature(gene_name, record.extract(genome.seq), record.location,
                                        record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                        genome_name, genome)
        gene_features[gene_name] = feature
        print(record.qualifiers['product'][0])
        tRNA_list.append(record.qualifiers['product'][0])
        list_intergenic_locations.append(record.location)
        list_intergenic_name.append(gene_name)
        tRNA_count = tRNA_count + 1

    elif record.type == 'rRNA':
        gene_name = record.qualifiers['product'][0]
        feature = mtDNA_SequenceFeature(gene_name, record.extract(genome.seq), record.location,
                                        record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                        genome_name, genome)
        gene_features[record.qualifiers['product'][0]] = feature
        print(record.qualifiers['product'][0])
        list_intergenic_locations.append(record.location)
        list_intergenic_name.append(gene_name)
        rRNA_count = rRNA_count + 1

    elif record.type == 'ncRNA':
        gene_name = record.qualifiers['product'][0]
        feature = mtDNA_SequenceFeature(gene_name, record.extract(genome.seq), record.location,
                                        record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                        genome_name, genome)
        gene_features[record.qualifiers['product'][0]] = feature
        print(record.qualifiers['product'][0])
        ncRNA_count = ncRNA_count + 1
        list_intergenic_locations.append(record.location)
        list_intergenic_name.append(gene_name)

    else:
        no1.append(record.location) #maybe filter these for things that have already appeared?

print(list_intergenic_locations)
print(100*'-')
#intergenic features
intergenic_seq = []
intergenic_features = OrderedDict()
for i in range(len(list_intergenic_locations)):
    j = i + 1 #to subtract one gene location from another

    if i == len(list_intergenic_locations):
        break #avoids errors i assume

    elif i == 0:
        intergenic = genome.seq[:list_intergenic_locations[i].start]
        intergenic_seq.append(intergenic)
        intergenic_begin = intergenic
        intergenic = genome.seq[list_intergenic_locations[i].end:list_intergenic_locations[j].start]
        intergenic_seq.append(intergenic)
        intergenic_object = mtDNA_IntergenicFeature(
            str(str(list_intergenic_name[i])) + "-" + str(list_intergenic_name[j]), intergenic,
            list_intergenic_locations[i].end, list_intergenic_locations[j].start, genome_name, genome)
        intergenic_features[
            str(list_intergenic_name[i]) + "-" + str(list_intergenic_name[j])] = intergenic_object

    elif i == len(list_intergenic_locations) - 1: #for later creation of intergenic feature between last gene and first gene
        intergenic = genome.seq[list_intergenic_locations[i].end:]
        intergenic_seq.append(intergenic)
        intergenic_end = intergenic

        continue

    else:
        intergenic = genome.seq[list_intergenic_locations[i].end:list_intergenic_locations[j].start]
        intergenic_seq.append(intergenic)
        intergenic_object = mtDNA_IntergenicFeature(
            str(str(list_intergenic_name[i])) + "-" + str(list_intergenic_name[j]), intergenic,
            list_intergenic_locations[i].end, list_intergenic_locations[j].start, genome_name, genome)
        intergenic_features[
            str(list_intergenic_name[i]) + "-" + str(list_intergenic_name[j])] = intergenic_object

intergenic_fin = intergenic_end + intergenic_begin #last intergenic features (last - first)
intergenic_seq.append(intergenic_fin)

intergenic_object = mtDNA_IntergenicFeature(
    str(str(list_intergenic_name[len(list_intergenic_name) - 1]) + "-" + str(list_intergenic_name[0])),
    intergenic_fin, list_intergenic_locations[len(list_intergenic_locations) - 1], list_intergenic_locations[0],
    genome_name, genome)

intergenic_features[str(list_intergenic_name[len(list_intergenic_name) - 1]) + "-" + str(
    list_intergenic_name[0])] = intergenic_object

# # picking out things in no1 list, that are associated with a feature record
# no1_location = []
# list_location1 = []
# for i in range(len(no1)):
#     start = no1[i].location.start
#     end = no1[i].location.end
#     location = str(str(start) + "-" + str(end))
#     no1_location.append(location)
#
# for j in range(len(list_intergenic_locations)):
#     start1 = list_intergenic_locations[j].start
#     end1 = list_intergenic_locations[j].end
#     location1 = str(str(start1) + "-" + str(end1))
#     list_location1.append(location1)
#
# for i in no1_location:
#     if i in list_location1:
#         continue
#     else:
#         sep = '-'
#         start = i.split(sep, 1)[0]
#         end = i.split(sep, 1)[1]
#         no2.append(start)
#         no2.append(end)
#
# for i in range(len(no2)):
#     if i == len(no2):
#         break
#
#     elif i % 2 == 0:
#         start = no2[i]
#         end = no2[i + 1]
#         for (n, x) in enumerate(genome.features):
#             record = genome.features[n]
#             if str(record.location.start) == start and str(record.location.end) == end:
#                 no_final.append(record)
#
#     else:
#         continue
# count = count + 1

# master dictionary with all genomes
all_gene_features_dict[genome_name] = gene_features
all_intergenic_features_dict[genome_name] = intergenic_features
all_no_features_dict[genome_name] = no_final

parse.write("\n\n" + str(count) + ")" + str(genome_name) + ":\n" + "Number of protein coding genes: " + str(
    CDS_count) + "\nNumber of tRNA genes: " + str(tRNA_count) + "\nNumber of rRNA genes: " + str(
    rRNA_count) + "\nNumber of ncRNA genes: " + str(ncRNA_count) + "\nTotal number of features: " + str(
    len(genome.features)) + "\nNumber of gene features: " + str(
    len(gene_features)) + "\nNumber of features not in any of the categories: " + str(
    len(no_final)) + "\nNumber of intergenics: " + str(len(intergenic_features)))

parse.close()
problematic.close()

count1 = 0
for k, v in all_gene_features_dict.items():
    print(k)
    for j in all_gene_features_dict[k]:
        print(str(all_gene_features_dict[k][j].name) + ":" + str(GC(all_gene_features_dict[k][j].seq)))
        count1 = count1 + 1

print(count1)

for k, v in all_intergenic_features_dict.items():
    print("\n" + str(k))
    for j in all_intergenic_features_dict[k]:
        print(str(all_intergenic_features_dict[k][j].name) + ":" + str(GC(all_intergenic_features_dict[k][j].seq)))
        count1 = count1 + 1

print(all_gene_features_dict)
print(all_intergenic_features_dict)