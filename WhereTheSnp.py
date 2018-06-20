## WhereTheSNP.py ##
## Where the snp? Which gene it located in? Does it cause AA change? ##
## Author: Zhang Fei ##
## Please contact: zhangfei-123@foxmail.com ##
## 2018-06-17 ##

## The output of this script are examed and is same with that from VEP

import sys, re

fl_snpquery = sys.argv[1]
fl_gff = sys.argv[2]
fl_cds = sys.argv[3]
prefix = sys.argv[4]

## build the information dic set
## Please notice the coding type of cds fasta sequence. It may cause bug.

def my_gffdic_build(fl_gff, fl_cds):
    with open(fl_cds, 'r') as fh_cds:
        dic_cds_sequences, transcript_num = {}, 0
        for line in fh_cds:
            if line[0] == '>' and transcript_num == 0:
                transcript = line.split(' ')[0][1:]
                sequences_temp = ''
                ## 0 means this is the initial transcript, and 1 means this is not the initial transcript
                transcript_num = 1
            elif line[0] != '>':
                sequences_temp = sequences_temp + line.strip('\n')
            elif line[0] == '>' and transcript_num == 1:
                dic_cds_sequences[transcript] = sequences_temp
                transcript = line.split(' ')[0][1:]
                sequences_temp = ''
        dic_cds_sequences[transcript] = sequences_temp

    dic_geneinfo, dic_gene_transcript, dic_transcript_info, dic_transcript_cds, dic_cds_transcript, dic_cds_info = {}, {}, {}, {}, {}, {}
    with open(fl_gff, 'r') as fh_gff:
        for line in fh_gff:
            line = line.strip('\n').split('\t')
            if line[2] == 'gene':
                chrom_gene, start_gene, end_gene, direction_gene = line[0], int(line[3]), int(line[4]), line[6]
                gene = re.split(';|=', line[8])[1]
                dic_geneinfo[gene] = [chrom_gene, start_gene, end_gene, direction_gene]
            elif line[2] == 'mRNA':
                chrom_mRNA, start_mRNA, end_mRNA, direction_mRNA = line[0], int(line[3]), int(line[4]), line[6]
                mRNA = re.split(';|=', line[8])[1]
                if mRNA[0:2] == 'GR':
                    motherGene = mRNA.split('_')[0]
                    if dic_gene_transcript.get(motherGene) is None:
                        dic_gene_transcript[motherGene] = [mRNA]
                    else:
                        dic_gene_transcript[motherGene].append(mRNA)
                else:
                    motherGene = mRNA.replace('T', '')
                    if dic_gene_transcript.get(motherGene) is None:
                        dic_gene_transcript[motherGene] = [mRNA]
                    else:
                        dic_gene_transcript[motherGene].append(mRNA)
                dic_transcript_info[mRNA] = [chrom_mRNA, start_mRNA, end_mRNA, int(line[4]) - int(line[3]), direction_mRNA]
            elif line[2] == 'CDS':
                chrom_cds, start_cds, end_cds = line[0], line[3], line[4]
                motherTranscript = re.split(';|=', line[8])[1]
                cdsID = re.split(';|=', line[8])[3]
                if dic_transcript_cds.get(motherTranscript) is None:
                    dic_transcript_cds[motherTranscript] = [cdsID]
                else:
                    dic_transcript_cds[motherTranscript].append(cdsID)
                if dic_cds_transcript.get(cdsID) is None:
                    dic_cds_transcript[cdsID] = [motherTranscript]
                else:
                    dic_cds_transcript[cdsID].append(motherTranscript)
                dic_cds_info[cdsID] = chrom_cds + ":" + start_cds + ":" + end_cds

    dic_result = {'dic_geneinfo':dic_geneinfo, 'dic_gene_transcript':dic_gene_transcript, 'dic_transcript_info': dic_transcript_info, 'dic_cds_sequences': dic_cds_sequences, \
            'dic_transcript_cds': dic_transcript_cds, 'dic_cds_transcript':dic_cds_transcript, 'dic_cds_info':dic_cds_info}
    return dic_result

def find_longest_transcript(gene, dic_gene_transcript, dic_transcript_info):
    transcript_list = dic_gene_transcript[gene]
    if len(transcript_list) == 1:
        return transcript_list[0]
    else:
        length_list = [dic_transcript_info[i][3] for i in transcript_list]
        index_longest = length_list.index(max(length_list))
        return transcript_list[index_longest]

def codon_detect(snp_coord, cds_sequences):
    ## cds_sequences should all from ATG to stop codon
    ## in this situation, we ingore of the situation of pharse and consider all protein start at ATG
    ## in the further we can add the parse infomation
    ## python consider the start at 0, but common case consider as 1, so we should do the corresponded adjustment
    snp_coord = snp_coord - 1
    condition_snp = snp_coord % 3
    if condition_snp == 0:
        coord_snp = [snp_coord, snp_coord + 1, snp_coord + 2]
        pharse = 0
    elif condition_snp == 1:
        coord_snp = [snp_coord - 1, snp_coord, snp_coord + 1]
        pharse = 1
    elif condition_snp == 2:
        coord_snp = [snp_coord - 2, snp_coord - 1, snp_coord]
        pharse = 2
    codon = cds_sequences[coord_snp[0]] + cds_sequences[coord_snp[1]] + cds_sequences[coord_snp[2]]
    return [codon, pharse]

def promoter_detect(snpchrom, snplocation, dic_gff):
    record = [['-', '-']]
    for line_gene in dic_gff['dic_geneinfo'].keys():
        entry = dic_gff['dic_geneinfo'][line_gene]
        gene_chrom, gene_start, gene_end, gene_direction = entry[0], entry[1], entry[2], entry[3]
        if gene_chrom == snpchrom and gene_direction == '+':
            if gene_start - snplocation >= 0 and gene_start - snplocation <= 2000:
                distance = -(gene_start - snplocation + 1)
                record.append(['promoter_region', line_gene, distance])
        elif gene_chrom == snpchrom and gene_direction == '-':
            if snplocation - gene_end >= 0 and snplocation - gene_end <= 2000:
                distance = -(snplocation - gene_end + 1)
                record.append(['promoter_region', line_gene, distance])
    return record

## This following  part gives the location of SNP and has passed the test
## if we start at 0, ensure the snp coord by the step following steps and find the nucleotides by the SNP_coord (start at 01) of course, and we will find the real location where the SNP are
## if gene in minus chain, the SNP will be the complementary one of hmp data record. But we are sure to follow the nucleotide we've found to calculate AA
## However, we need a input of SNP type in this locus to test whether AA change has happended, condiser the plus/minus chain the gene located
## All SNP input should be record as the type in plus chain, this function could be achieved in a special function

def cds_detect(snpchrom, snplocation, gene, dic_gff):
    transcript = find_longest_transcript(gene, dic_gff['dic_gene_transcript'], dic_gff['dic_transcript_info'])
    cds_sequences = dic_gff['dic_cds_sequences'][transcript]
    #print cds_sequences
    entry_transcript = dic_gff['dic_transcript_info'][transcript]
    transcript_chrom, transcript_start, transcript_end, transcript_direction = entry_transcript[0], entry_transcript[2], entry_transcript[3], entry_transcript[4]
    cds_list = dic_gff['dic_transcript_cds'][transcript]
    whetherincds, cds_start_list, cds_end_list = 'F', [], []
    for line_cd in cds_list:
        cds_entry = dic_gff['dic_cds_info'][line_cd].split(':')
        #print cds_entry
        cds_chrom, cds_start, cds_end = cds_entry[0], int(cds_entry[1]), int(cds_entry[2])
        if snpchrom == cds_chrom and snplocation >= cds_start and snplocation <= cds_end:
            whetherincds = 'T'
        cds_start_list.append(cds_start)
        cds_end_list.append(cds_end)
    cds_start_list.sort()
    cds_end_list.sort()
    cds_length_list = [cds_end_list[i] - cds_start_list[i] for i in range(len(cds_end_list))]
    #print cds_start_list, cds_end_list, cds_length_list
    ## further optimization: package this part as a function()
    #print whetherincds
    if whetherincds == 'F':
        record = [gene, transcript_direction, 'not in cds']
        return record
    elif whetherincds == 'T' and transcript_direction == '+':
        cds_start_list, cds_end_list, cds_length_list, snp_coord = cds_start_list, cds_end_list, cds_length_list, 0
        for i in range(len(cds_start_list)):
            if snplocation >= cds_start_list[i] and snplocation <= cds_end_list[i]:
                SNP_relative_locus = snplocation - cds_start_list[i] + 1
                snp_coord = snp_coord + SNP_relative_locus
                SNP_reference = cds_sequences[snp_coord - 1]
                codon_info = codon_detect(snp_coord, cds_sequences)
                codon_reference = codon_info[0]
                codon_pharse = codon_info[1]
                break
            else:
                ## NOTICE: Off-by-one error
                snp_coord = snp_coord + cds_length_list[i] + 1
        #print snp_coord, SNP_reference, codon_reference
    
    elif whetherincds == 'T' and transcript_direction == '-':
        cds_start_list.reverse()
        cds_end_list.reverse()
        cds_length_list.reverse()
        cds_start_list, cds_end_list, snp_coord = cds_end_list, cds_start_list, 0
        for i in range(len(cds_start_list)):
            if snplocation <= cds_start_list[i] and snplocation >= cds_end_list[i]:
                SNP_relative_locus = cds_start_list[i] - snplocation + 1
                snp_coord = snp_coord + SNP_relative_locus
                ## snp_coord do give the right location of SNP in coding region
                SNP_reference = cds_sequences[snp_coord - 1]
                codon_info = codon_detect(snp_coord, cds_sequences)
                codon_reference = codon_info[0]
                codon_pharse = codon_info[1]
                break
            else:
                ## NOTICE: Off-by-one error
                snp_coord = snp_coord + cds_length_list[i] + 1
        #print snp_coord, SNP_reference, codon_reference
    record = [gene, transcript_direction, 'cds', snp_coord, SNP_reference, codon_reference, codon_pharse]
    return record

## this function is only meaningful when apply to the SNP in cds
## input:
## record: the output put function 'cds_detect'
## snprecor: a list of SNP type in this loci. e.g [T, C]
def AAchange_detect(record, snprecord):
    dic_complementary = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    dic_codon_DNA={'TTT':'F','CTT':'L','ATT':'I','GTT':'V','TTC':'F','CTC':'L','ATC':'I','GTC':'V','TTA':'L', \
            'CTA':'L','ATA':'I','GTA':'V','TTG':'L','CTG':'L','ATG':'M','GTG':'V','TCT':'S','CCT':'P','ACT':'T', \
            'GCT':'A','TCC':'S','CCC':'P','ACC':'T','GCC':'A','TCA':'S','CCA':'P','ACA':'T','GCA':'A', \
            'TCG':'S','CCG':'P','ACG':'T','GCG':'A','TAT':'Y','CAT':'H','AAT':'N','GAT':'D', \
            'TAC':'Y','CAC':'H','AAC':'N','GAC':'D','TAA':'Stop','CAA':'Q','AAA':'K','GAA':'E', \
            'TAG':'Stop','CAG':'Q','AAG':'K','GAG':'E','TGT':'C','CGT':'R','AGT':'S','GGT':'G', \
            'TGC':'C','CGC':'R','AGC':'S','GGC':'G','TGA':'Stop','CGA':'R','AGA':'R','GGA':'G', \
            'TGG':'W','CGG':'R','AGG':'R','GGG':'G'
            }
    RefSNP, RefCodon, snpPharse = record[4], record[5], record[6]
    RefAA = dic_codon_DNA[RefCodon]
    if record[1] == '+':
        AltSNP = list(set(snprecord) - set(RefSNP))[0]
        if snpPharse == 0:
            AltCodon = AltSNP + RefCodon[1] + RefCodon[2]
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = AltSNP + RefCodon[1].lower() + RefCodon[2].lower()
            RefCodon_out = RefCodon[0] + RefCodon[1].lower() + RefCodon[2].lower()
        elif snpPharse == 1:
            AltCodon = RefCodon[0] + AltSNP + RefCodon[2]
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = RefCodon[0].lower() + AltSNP + RefCodon[2].lower()
            RefCodon_out = RefCodon[0].lower() + RefCodon[1] + RefCodon[2].lower()
        elif snpPharse == 2:
            AltCodon = RefCodon[0] + RefCodon[1] + AltSNP
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = RefCodon[0].lower() + RefCodon[1].lower() + AltSNP
            RefCodon_out = RefCodon[0].lower() + RefCodon[1].lower() + RefCodon[2]
    elif record[1] == '-':
        RefSNP, RefCodon, snpPharse = record[4], record[5], record[6]
        snprecord = [dic_complementary[i] for i in snprecord]
        AltSNP = list(set(snprecord) - set(RefSNP))[0]
        if snpPharse ==0:
            AltCodon = AltSNP + RefCodon[1] + RefCodon[2]
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = AltSNP + RefCodon[1].lower() + RefCodon[2].lower()
            RefCodon_out = RefCodon[0] + RefCodon[1].lower() + RefCodon[2].lower()
        elif snpPharse == 1:
            AltCodon = RefCodon[0] + AltSNP + RefCodon[2]
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = RefCodon[0].lower() + AltSNP + RefCodon[2].lower()
            RefCodon_out = RefCodon[0].lower() + RefCodon[1] + RefCodon[2].lower()
        elif snpPharse == 2:
            AltCodon = RefCodon[0] + RefCodon[1] + AltSNP
            AltAA = dic_codon_DNA[AltCodon]
            AltCodon_out = RefCodon[0].lower() + RefCodon[1].lower() + AltSNP
            RefCodon_out = RefCodon[0].lower() + RefCodon[1].lower() + RefCodon[2]
    result = record + ['/'.join(snprecord)] + [RefCodon_out, RefAA, AltSNP, AltCodon_out, AltAA]
    return result

#### main function ####
dic_gff = my_gffdic_build(fl_gff, fl_cds)
result_promoter_region = [['snpname', 'snpchrom', 'snplocation', 'snptype', 'region', 'gene', 'distance']]
result_not_in_cds = [['snpname', 'snpchrom', 'snplocation', 'snptype', 'region', 'gene', 'gene_direction']]
result_cds = [['snpname', 'snpchrom', 'snplocation', 'snptype', 'region', 'gene', 'gene_direction', 'snpcoord', 'Ref.codon', 'Alt.codon', 'Ref.aa', 'Alt.aa']]
result_left = [['snpname', 'snpchrom', 'snplocation', 'snptype']]

with open(fl_snpquery, 'r') as fh_snpquery:
    for snp in fh_snpquery:
        snp = snp.strip('\r\n').split('\t')
        snpname, snpchrom, snplocation, gene, snpdirection, snprecord = snp[0], snp[1], int(snp[2]), 'NA', snp[4], re.split('/', snp[3])
        for line_gene in dic_gff['dic_geneinfo'].keys():
            entry = dic_gff['dic_geneinfo'][line_gene]
            gene_chrom, gene_start, gene_end, gene_direction = entry[0], entry[1], entry[2], entry[3]
            if snpchrom == gene_chrom and snplocation >= gene_start and snplocation <= gene_end:
                gene = line_gene
                continue
        if gene != 'NA':
#            print '###################'
#            print snpname, "start"
#            print gene
            record = cds_detect(snpchrom, snplocation, gene, dic_gff)
            if record[2] == 'cds':
                AA_detect_result = AAchange_detect(record, snprecord)
                result = [snpname, snpchrom, snplocation, snp[3], 'CDS'] + AA_detect_result[0:2] + [AA_detect_result[3], AA_detect_result[8], AA_detect_result[9], AA_detect_result[11], AA_detect_result[12]]
                result = [str(i) for i in result]
                result_cds.append(result)
            elif record[2] == 'not in cds':
                result = [snpname, snpchrom, snplocation, snp[3], 'not in cds'] + record[0:2]
                result = [str(i) for i in result]
                result_not_in_cds.append(result)
        if gene == 'NA':
#            print '##################'
#            print snpname, "start"
            record = promoter_detect(snpchrom, snplocation, dic_gff)
            for line in record:
                if line[0] == 'promoter_region':
                    result = [snpname, snpchrom, snplocation, snp[3]] + line
                    result = [str(i) for i in result]
                    result_promoter_region.append(result)
                elif line[0] == '-':
                    result = [snpname, snpchrom, snplocation, snp[3]]
                    result = [str(i) for i in result]
                    result_left.append(result)

with open(prefix + '_promoter_region.txt', 'w') as fh_promoter_region:
    result_promoter_region = ['\t'.join(i) for i in result_promoter_region]
    fh_promoter_region.write('\n'.join(result_promoter_region))
with open(prefix + '_not_in_cds.txt', 'w') as fh_not_in_cds:
    result_not_in_cds = ['\t'.join(i) for i in result_not_in_cds]
    fh_not_in_cds.write('\n'.join(result_not_in_cds))
with open(prefix + '_cds.txt', 'w') as fh_cds:
    result_cds = ['\t'.join(i) for i in result_cds]
    fh_cds.write('\n'.join(result_cds))
with open(prefix + '_others.txt', 'w') as fh_others:
    result_left = ['\t'.join(i) for i in result_left]
    fh_others.write('\n'.join(result_left))
