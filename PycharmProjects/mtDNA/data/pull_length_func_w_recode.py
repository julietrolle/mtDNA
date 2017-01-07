from Bio import SeqIO
from collections import OrderedDict
import os
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

class mtDNA_SequenceChar(object):
    'sequence features for mtDNA project CDS,tRNA, rRNA and ncRNA genes'
    def __init__(self, name, type,seq, location, locus_tag, db_xref, genome_name, whole_genome_record):
        self.name = name
        self.seq = seq
        self.type = type
        self.locus_tag = locus_tag
        self.location = location
        self.db_xref = db_xref
        self.genome_name = genome_name
        self.whole_genome_record = whole_genome_record
    def __repr__(self):
        return(str(self.name)+"-"+str(self.genome_name))

class mtDNA_IntergenicChar(object):
    'sequence features for mtDNA project for intergenic spaces'
    def __init__(self, name,seq, location_start,location_end, genome_name, whole_genome_record):
        self.name = name
        self.seq = seq
        self.location_start = location_start
        self.location_end = location_end
        self.genome_name = genome_name
        self.whole_genome_record = whole_genome_record
    def __repr__(self):
        return(str(self.name)+"-"+str(self.genome_name))

def gene_char_assignment(gb_file):
    CDS_count = 0
    tRNA_count = 0
    rRNA_count = 0
    ncRNA_count = 0
    everything_else_count = 0
    for (i, f) in enumerate(gb_file.features):
        record = gb_file.features[i]
        if record.type == 'CDS':
            try:
                if record.qualifiers['gene'] in gene_list:
                    gene_name = record.qualifiers['gene'][0]
                    char = mtDNA_SequenceChar(gene_name, record.type, record.extract(gb_file.seq), record.location,
                                                    record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                                    genome_name, handle)

                    gene_chars[gene_name] = char

                    list_gene_locations.append(record.location)
                    list_gene_name.append(gene_name)
                    CDS_count = CDS_count + 1
            except:
                if record.type == 'ncRNA':
                    continue

        elif record.type == 'tRNA':
            seq = record.extract(gb_file.seq)
            gene_name = record.qualifiers['product'][0]

            if record.qualifiers['product'][0] in tRNA_list:
                gene_name = str(record.qualifiers['product'][0] + '2')

                char = mtDNA_SequenceChar(gene_name, record.type, record.extract(gb_file.seq), record.location,
                                                record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                                genome_name, handle)
                gene_chars[gene_name] = char
                list_gene_locations.append(record.location)
                list_gene_name.append(gene_name)
                tRNA_count = tRNA_count + 1
                continue

            feature = mtDNA_SequenceChar(gene_name, record.type, record.extract(gb_file.seq), record.location,
                                            record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                            genome_name, handle)
            gene_chars[gene_name] = feature

            tRNA_list.append(record.qualifiers['product'][0])
            list_gene_locations.append(record.location)
            list_gene_name.append(gene_name)
            tRNA_count = tRNA_count + 1

        elif record.type == 'rRNA':
            gene_name = record.qualifiers['product'][0]
            feature = mtDNA_SequenceChar(gene_name, record.type, record.extract(gb_file.seq), record.location,
                                            record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                            genome_name, handle)
            gene_chars[record.qualifiers['product'][0]] = feature

            list_gene_locations.append(record.location)
            list_gene_name.append(gene_name)
            rRNA_count = rRNA_count + 1

        elif record.type == 'ncRNA':
            gene_name = record.qualifiers['product'][0]
            feature = mtDNA_SequenceChar(gene_name, record.type, record.extract(gb_file.seq), record.location,
                                            record.qualifiers['locus_tag'][0], record.qualifiers['db_xref'][0],
                                            genome_name, handle)
            gene_chars[record.qualifiers['product'][0]] = feature

            ncRNA_count = ncRNA_count + 1
            list_gene_locations.append(record.location)
            list_gene_name.append(gene_name)

        else:
            no1.append(record)
    return gene_chars

def intergen_char_assignment(gb_file):
    for i in range(len(list_gene_locations)):
        j = i + 1

        if i == len(list_gene_locations): #stops loop so it doesn't attempt with non-existent j's
            break

        elif i == 0:
            intergenic = handle.seq[:list_gene_locations[i].start] #from 0 through beginning location of gene1
            intergenic_seq.append(intergenic)
            intergenic_begin = intergenic

            intergenic = handle.seq[list_gene_locations[i].end:list_gene_locations[j].start] #from gene a thru b
            intergenic_seq.append(intergenic)
            intergenic_object = mtDNA_IntergenicChar(
                str(str(list_gene_name[i])) + "-" + str(list_gene_name[j]), intergenic,
                list_gene_locations[i].end, list_gene_locations[j].start, genome_name, handle)

            intergenic_chars[str(list_gene_name[i]) + "-" + str(list_gene_name[j])] = intergenic_object

        elif i == len(list_gene_locations) - 1:
            intergenic = handle.seq[list_gene_locations[i].end:]
            intergenic_seq.append(intergenic)
            intergenic_end = intergenic
            continue


        else:
            intergenic = handle.seq[list_gene_locations[i].end:list_gene_locations[j].start]
            intergenic_seq.append(intergenic)
            intergenic_object = mtDNA_IntergenicChar(
                str(str(list_gene_name[i])) + "-" + str(list_gene_name[j]), intergenic,
                list_gene_locations[i].end, list_gene_locations[j].start, genome_name, handle)

            intergenic_chars[str(list_gene_name[i]) + "-" + str(list_gene_name[j])] = intergenic_object

    intergenic_fin = intergenic_end + intergenic_begin #joins first and last
    intergenic_seq.append(intergenic_fin)
    intergenic_object = mtDNA_IntergenicChar(
        str(str(list_gene_name[len(list_gene_name) - 1]) + "-" + str(list_gene_name[0])),
        intergenic_fin, list_gene_locations[len(list_gene_locations) - 1], list_gene_locations[0],
        genome_name, handle)

    intergenic_chars[str(list_gene_name[len(list_gene_name) - 1]) + "-" + str(
        list_gene_name[0])] = intergenic_object

    return intergenic_chars

def func_for_no1(gb_file):
    count = 0
    for i in range(len(no1)):
        start = no1[i].location.start
        end = no1[i].location.end
        location = str(str(start) + "-" + str(end))
        no1_location.append(location)

    for j in range(len(list_gene_locations)):
        start1 = list_gene_locations[j].start
        end1 = list_gene_locations[j].end
        location1 = str(str(start1) + "-" + str(end1))
        list_location1.append(location1)

    for i in no1_location:
        if i in list_location1:
            continue
        else:
            sep = '-'
            start = i.split(sep, 1)[0]
            end = i.split(sep, 1)[1]
            no2.append(start)
            no2.append(end)

    for i in range(len(no2)):
        if i == len(no2):
            break

        elif i % 2 == 0:
            start = no2[i]
            end = no2[i + 1]
            for (n, x) in enumerate(handle.features):
                record = handle.features[n]
                if str(record.location.start) == start and str(record.location.end) == end:
                    no_final.append(record)

        else:
            continue
    count = count + 1
    return count #not too sure about this

def split_into_codons(gene_seq, n):
    list_of_codons = []
    for i in range(0, len(gene_seq), n):
        list_of_codons.append(gene_seq[i:i+n])
    return list_of_codons

def high_gc_recode(codon_list):
    new_codon_list = []
    for codon in codon_list:
        if codon == 'TTT': #Phe
            new_codon_list.append(Seq('TTC', IUPAC.ambiguous_dna))
        elif codon == 'TTA': #Leu
            new_codon_list.append(Seq('TTG', IUPAC.ambiguous_dna))
        elif codon == 'ACT': #Thr
            new_codon_list.append(Seq('ACC', IUPAC.ambiguous_dna))
        elif codon == 'ACA': #Thr
            new_codon_list.append(Seq('ACG', IUPAC.ambiguous_dna))
        elif codon == 'CTT': #Thr
            new_codon_list.append(Seq('CTC', IUPAC.ambiguous_dna))
        elif codon == 'CTA': #Thr
            new_codon_list.append(Seq('CTG', IUPAC.ambiguous_dna))
        elif codon == 'ATT': #Ile
            new_codon_list.append(Seq('ATC', IUPAC.ambiguous_dna))
        elif codon == 'ATA': #Met / Start
            new_codon_list.append(Seq('ATG', IUPAC.ambiguous_dna))
        elif codon == 'GTT': #Val
            new_codon_list.append(Seq('GTC', IUPAC.ambiguous_dna))
        elif codon == 'GTA': #Val
            new_codon_list.append(Seq('GTG', IUPAC.ambiguous_dna))
        elif codon == 'AGT': #Ser
            new_codon_list.append(Seq('AGC', IUPAC.ambiguous_dna))
        elif codon == 'TCT': #Ser
            new_codon_list.append(Seq('TCC', IUPAC.ambiguous_dna))
        elif codon == 'TCA': #Ser
            new_codon_list.append(Seq('TCG', IUPAC.ambiguous_dna))
        elif codon == 'CCT': #Pro
            new_codon_list.append(Seq('CCC', IUPAC.ambiguous_dna))
        elif codon == 'CCA': #Pro
            new_codon_list.append(Seq('CCG', IUPAC.ambiguous_dna))
        elif codon == 'GCT': #Ala
            new_codon_list.append(Seq('GCC', IUPAC.ambiguous_dna))
        elif codon == 'GCA': #Ala
            new_codon_list.append(Seq('GCG', IUPAC.ambiguous_dna))
        elif codon == 'TAT': #Tyr
            new_codon_list.append(Seq('TAC', IUPAC.ambiguous_dna))
        elif codon == 'TAA': #Stop
            new_codon_list.append(Seq('TAG', IUPAC.ambiguous_dna))
        elif codon == 'CAT': #His
            new_codon_list.append(Seq('CAC', IUPAC.ambiguous_dna))
        elif codon == 'CAA': #Gln
            new_codon_list.append(Seq('CAG', IUPAC.ambiguous_dna))
        elif codon == 'AAT': #Asn
            new_codon_list.append(Seq('AAC', IUPAC.ambiguous_dna))
        elif codon == 'AAA': #Lys
            new_codon_list.append(Seq('AAG', IUPAC.ambiguous_dna))
        elif codon == 'GAT': #Asp
            new_codon_list.append(Seq('GAC', IUPAC.ambiguous_dna))
        elif codon == 'GAA': #Glu
            new_codon_list.append(Seq('GAG', IUPAC.ambiguous_dna))
        elif codon == 'TGT': #Cys
            new_codon_list.append(Seq('TGC', IUPAC.ambiguous_dna))
        elif codon == 'TGA': #Trp
            new_codon_list.append(Seq('TGG', IUPAC.ambiguous_dna))
        elif codon == 'CGT': #Arg
            new_codon_list.append(Seq('CGC', IUPAC.ambiguous_dna))
        elif codon == 'CGA': #Arg
            new_codon_list.append(Seq('CGG', IUPAC.ambiguous_dna))
        elif codon == 'AGG': #Arg (EXCEPTION)
            new_codon_list.append(Seq('CGC', IUPAC.ambiguous_dna))
        elif codon == 'AGA': #Arg (EXCEPTION)
            new_codon_list.append(Seq('CGG', IUPAC.ambiguous_dna))
        elif codon == 'GGT': #Gly
            new_codon_list.append(Seq('GGC', IUPAC.ambiguous_dna))
        elif codon == 'GGA': #Gly
            new_codon_list.append(Seq('GGG', IUPAC.ambiguous_dna))
        else:
            new_codon_list.append(codon)
    return new_codon_list

def codons_to_seq_obj(codon_list):
    seq = Seq('', IUPAC.ambiguous_dna)
    for codon in codon_list:
        seq = seq + codon
    return seq

parse = open('parsing_stats.txt','w')
problematic = open('problematic_files.txt','w')
all_length = open ('all_lengths.txt','w')
all_GC = open ('all_GCs.txt','w')

all_gene_chars_dict = {}
all_intergen_chars_dict = {}
all_no_chars_dict = {}
count = 0

list_dir_files = os.listdir('/Users/julietrolle/PycharmProjects/mtDNA/data/genbank_files')

list_gb_files = [] #to ensure that only gb files are opened
for element in list_dir_files:
    if element.endswith('gb') == True:
        list_gb_files.append(element)

for gb in list_gb_files:
    try:
        f = open('/Users/julietrolle/PycharmProjects/mtDNA/data/genbank_files/%s' % gb, 'r')
        handle = SeqIO.read(f, 'genbank')
        genome_name = gb.split('.', 1)[0]

        no1 = []
        no2 = []
        no_final = []

        gene_list = [['COX1'], ['ATP8'], ['ATP6'], ['COB'], ['OLI1'], ['VAR1'], ['COX2'], ['COX3']]
        tRNA_list = []
        list_gene_locations = []  # to extract intergenic locations
        list_gene_name = []  # to make intergenic names
        gene_chars = OrderedDict()  # store all gene chars by key(gene_name)

        gene_func = gene_char_assignment(handle)

        intergenic_seq = []
        intergenic_chars = OrderedDict()

        intergen_func = intergen_char_assignment(handle)

        no1_location = []
        list_location1 = []

        func_for_no1(handle)

        all_gene_chars_dict[genome_name] = gene_func
        all_intergen_chars_dict[genome_name] = intergen_func
        all_no_chars_dict[genome_name] = no_final

        count = count + 1

        # parse.write("\n\n" + str(count) + ")" + str(genome_name) + ":\n" + "Number of protein coding genes: " + str(
        #     CDS_count) + "\nNumber of tRNA genes: " + str(tRNA_count) + "\nNumber of rRNA genes: " + str(
        #     rRNA_count) + "\nNumber of ncRNA genes: " + str(ncRNA_count) + "\nTotal number of features: " + str(
        #     len(handle.features)) + "\nNumber of gene features: " + str(
        #     len(gene_chars)) + "\nNumber of features not in any of the categories: " + str(
        #     len(no_final)) + "\nNumber of intergenics: " + str(len(intergenic_chars)))

    except ValueError:
         problematic.write("\n"+str(genome_name) + " ValueError")
         print(ValueError)
         continue
    except KeyError:
        problematic.write("\n" + str(genome_name) +" KeyError")
        print(ValueError)
        continue

# parse.close()
problematic.close()

top_len_chars_list = []
top_len_intergenic_chars_list = []
problem_len = open("len_problem.txt",'w')
file_len = open("len_stats.txt",'w')

#parses thru to find shortest gene
for i in all_gene_chars_dict['KP263414']:  # for gene feature i in this file (but all files have same features so works)
    len_top_seq = 1000000000000  # initialize loop
    for k, v in all_gene_chars_dict.items():  # for file, dict with each feature
        try:
            current_seq_len = len(all_gene_chars_dict[k][i].seq) #[genome][gene] --> gets you features for that gene
            if current_seq_len < len_top_seq:
                len_top_seq = current_seq_len
                top_len_record = all_gene_chars_dict[k][i] #for the smallest feature, input into top len rec list

        except KeyError:
            problem_len.write("\n" + str(v) + "KeyError-feature")
            continue

    top_len_chars_list.append(top_len_record)
    file_len.write("\n" + str(top_len_record) + " " + str(len(top_len_record.seq))) #puts name and len in lenstats

#recoding CDS
recoded_top_len_chars_list = []
for feature in top_len_chars_list:
    if feature.type == 'CDS':
        gene_CDS = feature.seq
        codons_gene_CDS = split_into_codons(gene_CDS, 3)
        recoded_codons_gene_CDS = high_gc_recode(codons_gene_CDS)
        recoded_gene_CDS = codons_to_seq_obj(recoded_codons_gene_CDS)
        recoded_top_len_chars_list.append(recoded_gene_CDS)
    else:
        recoded_top_len_chars_list.append(feature.seq)
        continue

#parses thru to find shortest intergenic region

file_len.write("\nIntergenic")

for i in all_intergen_chars_dict['KP263414']:
    len_top_seq_int = 100000000000000
    for k, v in all_intergen_chars_dict.items():
        try:
            current_seq_int_len = len(all_intergen_chars_dict[k][i].seq)
            if current_seq_int_len < len_top_seq_int:
                len_top_seq_int = current_seq_int_len
                top_len_record_int = all_intergen_chars_dict[k][i]

        except KeyError:
            problem_len.write("\n" + str(v) + "KeyError-intergenic")
            continue

    top_len_intergenic_chars_list.append(top_len_record_int)
    file_len.write("\n" + str(top_len_record_int) + " " + str(len(top_len_record_int.seq)) + ' ' + str(GC(top_len_record_int.seq))) #puts name and GC in lenstats, I'll add len too

file_len.close()
problem_len.close()

#starting point: list with seq+features of shortest intergenic sequences

final_record_len = SeqRecord('')
final_len = ''
final_list_len = []

for i in range(len(top_len_chars_list)):
    diction_len = {}
    diction_intergenic_len = {}  # {name : intergenic + number, genome name : genome of that feature}
    seq_feature_list_len = []
    seq_intergenic_feature_list_len = []

    g_record = str(top_len_chars_list[i].seq)
    gene_record = Seq(g_record, IUPAC.unambiguous_dna)
    name_gene = str("gene " + str(i))
    diction_len['note'] = [name_gene]
    genome_name = str(top_len_chars_list[i].genome_name)
    diction_len['genome_name'] = [genome_name]
    diction_len['product'] = top_len_chars_list[i].name
    diction_len['locus_tag'] = top_len_chars_list[i].name

    seq_feature = SeqFeature(FeatureLocation(0, len(top_len_chars_list[i].seq)), type=top_len_chars_list[i].type,
                             id=str(top_len_chars_list[i].name) + "-" + str(top_len_chars_list[i].genome_name),
                             qualifiers=diction_len)

    seq_feature_list_len.append(seq_feature)

    seq_record_gene = SeqRecord(gene_record, id='SC2mitoV2.1', name=name_gene,
                                description="Redesign of Saccharomyces cerevisiae mitochondrial genome using 100 yeast genomes resource",
                                features=seq_feature_list_len)

    final_list_len.append(seq_record_gene)

    record_intergenic = str(top_len_intergenic_chars_list[i].seq)
    record_intergenic = Seq(record_intergenic, IUPAC.unambiguous_dna) #makes a seq object out of each intergen feat
    name_intergenic = str("intergenic " + str(i)) #numbers each feature
    diction_intergenic_len['note'] = [name_intergenic]
    genome_name_intergenic = str(top_len_intergenic_chars_list[i].genome_name)
    diction_intergenic_len['genome_name'] = [genome_name_intergenic] #how does it distinguish one feature from another? dict not Ordered
    diction_intergenic_len['locus_tag'] = top_len_intergenic_chars_list[i].name

    seq_feature_intergenic = SeqFeature(FeatureLocation(0, len(top_len_intergenic_chars_list[i].seq)),
                                        type="intergenic", id=str(top_len_intergenic_chars_list[i].name) + "-" + str(
            top_len_intergenic_chars_list[i].genome_name), qualifiers=diction_intergenic_len)

    #why FeatureLocation(0, ...) ?

    seq_intergenic_feature_list_len.append(seq_feature_intergenic)
    seq_record_intergenic = SeqRecord(record_intergenic, id='SC2mitoV2.1', name=name_intergenic,
                                      description="Redesign of Saccharomyces cerevisiae mitochondrial genome using 100 yeast genomes resource",
                                      features=seq_intergenic_feature_list_len)
    final_list_len.append(seq_record_intergenic)
    final_len = final_len + gene_record + record_intergenic #str
    final_record_len = final_record_len + seq_record_gene + seq_record_intergenic #seqRec

output_handle = open("final_mito_len.gb",'w')
SeqIO.write(final_record_len, output_handle,'genbank') #SeqRecord
output_handle.close()

#NB: lots of repeating variable names

final_record_len = SeqRecord('')
final_len = ''
final_list_len = []

for i in range(len(top_len_chars_list)):
    diction_len = {}
    diction_intergenic_len = {}  # {name : intergenic + number, genome name : genome of that feature}
    seq_feature_list_len = []
    seq_intergenic_feature_list_len = []

    g_record = str(recoded_top_len_chars_list[i])
    gene_record = Seq(g_record, IUPAC.unambiguous_dna)
    name_gene = str("gene " + str(i))
    diction_len['note'] = [name_gene]
    genome_name = str(top_len_chars_list[i].genome_name)
    diction_len['genome_name'] = [genome_name]
    diction_len['product'] = top_len_chars_list[i].name
    diction_len['locus_tag'] = top_len_chars_list[i].name

    seq_feature = SeqFeature(FeatureLocation(0, len(top_len_chars_list[i].seq)), type=top_len_chars_list[i].type,
                             id=str(top_len_chars_list[i].name) + "-" + str(top_len_chars_list[i].genome_name),
                             qualifiers=diction_len)
    seq_feature_list_len.append(seq_feature)
    seq_record_gene = SeqRecord(gene_record, id='SC2mitoV2.1', name=name_gene,
                                description="Redesign of Saccharomyces cerevisiae mitochondrial genome using 100 yeast genomes resource",
                                features=seq_feature_list_len)
    final_list_len.append(seq_record_gene)

    record_intergenic = str(top_len_intergenic_chars_list[i].seq)
    record_intergenic = Seq(record_intergenic, IUPAC.unambiguous_dna) #makes a seq object out of each intergen feat
    name_intergenic = str("intergenic " + str(i)) #numbers each feature
    diction_intergenic_len['note'] = [name_intergenic]
    genome_name_intergenic = str(top_len_intergenic_chars_list[i].genome_name)
    diction_intergenic_len['genome_name'] = [genome_name_intergenic] #how does it distinguish one feature from another? dict not Ordered
    diction_intergenic_len['locus_tag'] = top_len_intergenic_chars_list[i].name


    seq_feature_intergenic = SeqFeature(FeatureLocation(0, len(top_len_intergenic_chars_list[i].seq)),
                                        type="intergenic", id=str(top_len_intergenic_chars_list[i].name) + "-" + str(
            top_len_intergenic_chars_list[i].genome_name), qualifiers=diction_intergenic_len)

    #why FeatureLocation(0, ...) ?

    seq_intergenic_feature_list_len.append(seq_feature_intergenic)
    seq_record_intergenic = SeqRecord(record_intergenic, id='SC2mitoV2.1', name=name_intergenic,
                                      description="Redesign of Saccharomyces cerevisiae mitochondrial genome using 100 yeast genomes resource",
                                      features=seq_intergenic_feature_list_len)
    final_list_len.append(seq_record_intergenic)
    final_len = final_len + gene_record + record_intergenic #str
    final_record_len = final_record_len + seq_record_gene + seq_record_intergenic #seqRec

output_handle = open("recoded_final_mito_len.gb",'w')
SeqIO.write(final_record_len, output_handle,'genbank') #SeqRecord
output_handle.close()

#intergenic variability

g = open('/Users/julietrolle/PycharmProjects/mtDNA/data/final_mito_len.gb', 'r')
han = SeqIO.read(g, 'genbank')

intergenic_features = []
longest_intergenic_features = []

for (i, f) in enumerate(han.features): #id longest intergenic features
    record = han.features[i]
    if record.type == 'intergenic':
        intergenic_features.append(record)
        if len(record) > 3900:
            print(record.qualifiers['note'], len(record))
            longest_intergenic_features.append(record)
    else:
        continue

# print(longest_intergenic_features) #longest features in final mito len: 0, 5, 10, 32

intergenic_0_seq = {}
intergenic_5_seq = {}
intergenic_10_seq = {}
intergenic_32_seq = {}

for i in all_intergen_chars_dict['KP263414']: #Create dict with all intergenic_xx
    for k, v in all_intergen_chars_dict.items():
        try:
            if all_intergen_chars_dict[k][i].name == 'tRNA-Pro-15S ribosomal RNA':
                intergenic_0_seq[k] = all_intergen_chars_dict[k][i].seq
            if all_intergen_chars_dict[k][i].name == 'ATP6-tRNA-Glu':
                intergenic_5_seq[k] = all_intergen_chars_dict[k][i].seq
            if all_intergen_chars_dict[k][i].name == 'VAR1-21S ribosomal RNA':
                intergenic_10_seq[k] = all_intergen_chars_dict[k][i].seq
            if all_intergen_chars_dict[k][i].name == 'COX3-tRNA-Met2':
                intergenic_32_seq[k] = all_intergen_chars_dict[k][i].seq
        except KeyError:
            continue

# for k, v in intergenic_32_seq.items(): #create gb file for all intergenic x, edit dict.items, seqrec name and file name
#     seq_feat_list = []
#     seq_feat = SeqFeature(FeatureLocation(0, len(v)), type='intergenic', id=k)
#     seq_feat_list.append(seq_feat)
#     seq_rec = SeqRecord(v, name='intergenic_32' + '-' + str(k), features=seq_feat_list)
#     outp_han = open('%s-intergenic_32.gb' % k, 'w')
#     SeqIO.write(seq_rec, outp_han, 'genbank')
#     outp_han.close()

#compare GC to length
#transcriptome annotations
