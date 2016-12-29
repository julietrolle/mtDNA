from Bio import SeqIO
from collections import OrderedDict
import os
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import IUPAC

f = open('/Users/julietrolle/PycharmProjects/mtDNA/data/final_mito_len.gb', 'r')
handle = SeqIO.read(f, 'genbank')

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

recoded_genes_list = []
for feature in handle.features:
    if feature.type == 'CDS':
        gene_CDS = handle.seq[feature.location.start.position:feature.location.end.position]
        codons_gene_CDS = split_into_codons(gene_CDS, 3)
        recoded_codons_gene_CDS = high_gc_recode(codons_gene_CDS)
        recoded_gene_CDS = codons_to_seq_obj(recoded_codons_gene_CDS)
        recoded_genes_list.append(recoded_gene_CDS)
    else:
        # recoded_genes_list.append(feature)
        continue

print(recoded_genes_list)