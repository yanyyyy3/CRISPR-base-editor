#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 16:59:30 2019

@author: yanying
"""
#import sys
from Bio import SeqIO
# from Bio.Seq import Seq
import argparse
import logging
import os
import csv
import time
from scipy.spatial.distance import hamming
start_time=time.time()
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="database gRNA file (fasta)")
parser.add_argument("-r1", default=None, help="Read1")
parser.add_argument("-l", type=int,default=20, help="Length of gRNAs")
parser.add_argument("-f", default='fastq', help="parse format")
parser.add_argument("-o", "--output", default="results", help="output file name")
Keio_essential_genes=['accA', 'accB', 'accC', 'accD', 'acpP', 'acpS', 'adk', 'alaS', 'alsK', 'argS', 'asd', 'asnS', 'aspS', 'bamA', 'bamD', 'bcsB', 'birA', 'btuB', 'can', 'cca', 'cdsA', 'chpS', 'coaA', 'coaD', 'coaE', 'ymfK', 'csrA', 'cydA', 'cydC', 'cysS', 'dapA', 'dapB', 'dapD', 'dapE', 'def', 'degS', 'der', 'dfp', 'dicA', 'djlB', 'dnaA', 'dnaB', 'dnaC', 'dnaE', 'dnaG', 'dnaN', 'dnaX', 'dut', 'dxr', 'dxs', 'eno', 'entD', 'era', 'erpA', 'fabA', 'fabB', 'fabD', 'fabG', 'fabI', 'fabZ', 'fbaA', 'ffh', 'fldA', 'fmt', 'folA', 'folC', 'folD', 'folE', 'folK', 'folP', 'frr', 'ftsA', 'ftsB', 'ftsE', 'ftsH', 'ftsI', 'ftsK', 'ftsL', 'ftsN', 'ftsQ', 'ftsW', 'ftsX', 'ftsY', 'ftsZ', 'fusA', 'gapA', 'glmM', 'glmS', 'glmU', 'glnS', 'gltX', 'glyQ', 'glyS', 'gmk', 'gpsA', 'groL', 'groS', 'grpE', 'gyrA', 'gyrB', 'hemA', 'hemB', 'hemC', 'hemD', 'hemE', 'hemG', 'hemH', 'hemL', 'hisS', 'holA', 'holB', 'ileS', 'infA', 'infB', 'infC', 'ispA', 'ispB', 'ispD', 'ispE', 'ispF', 'ispG', 'ispH', 'ispU', 'kdsA', 'kdsB', 'lepB', 'leuS', 'lexA', 'lgt', 'ligA', 'lnt', 'lolA', 'lolB', 'lolC', 'lolD', 'lolE', 'lptA', 'lptB', 'lptC', 'lptD', 'lptE', 'lptF', 'lptG', 'lpxA', 'lpxB', 'lpxC', 'lpxD', 'lpxH', 'lpxK', 'lspA', 'map', 'mazE', 'metG', 'metK', 'minD', 'minE', 'mlaB', 'mqsA', 'mraY', 'mrdA', 'mrdB', 'mreB', 'mreC', 'mreD', 'msbA', 'mukB', 'mukE', 'mukF', 'murA', 'murB', 'murC', 'murD', 'murE', 'murF', 'murG', 'murI', 'murJ', 'nadD', 'nadE', 'nadK', 'nrdA', 'nrdB', 'nusA', 'nusG', 'obgE', 'orn', 'parC', 'parE', 'pgk', 'pgsA', 'pheS', 'pheT', 'plsB', 'plsC', 'polA', 'ppa', 'prfA', 'prfB', 'priB', 'prmC', 'proS', 'prs', 'psd', 'pssA', 'pth', 'purB', 'pyrG', 'pyrH', 'racR', 'rho', 'ribA', 'ribB', 'ribC', 'ribD', 'ribE', 'ribF', 'rnc', 'rne', 'rnpA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplJ', 'rplK', 'rplL', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ', 'rplR', 'rplS', 'rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmH', 'rpoA', 'rpoB', 'rpoC', 'rpoD', 'rpoE', 'rpoH', 'rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE', 'rpsG', 'rpsH', 'rpsJ', 'rpsK', 'rpsL', 'rpsN', 'rpsP', 'rpsR', 'rpsS', 'rpsU', 'rseP', 'rsmI', 'secA', 'secD', 'secE', 'secF', 'secM', 'secY', 'serS', 'spoT', 'ssb', 'suhB', 'tadA', 'tdcF', 'thiL', 'thrS', 'tilS', 'tmk', 'topA', 'trmD', 'trpS', 'tsaB', 'tsaC', 'tsaD', 'tsaE', 'tsf', 'tyrS', 'ubiA', 'ubiB', 'ubiD', 'valS', 'waaA', 'waaU', 'wzyE', 'yagG', 'yceQ', 'ydfB', 'ydiL', 'yefM', 'yejM', 'yhbV', 'yhhQ', 'yiaD', 'yidC', 'yihA', 'yqgF', 'yrfF', 'zipA']
args = parser.parse_args()
db=args.db
read1=args.r1
l=args.l
parse_format=args.f
output_file_name = args.output
#memory=args.available_memory
if os.path.isdir(output_file_name)==False:
    os.mkdir(output_file_name)
else:
    os.system("rm -r %s"%output_file_name)
    os.mkdir(output_file_name)
def Substrings(sequence,length):
    substrings=[]
    for i in range(24,len(sequence)-length):
        substrings.append({'seq':sequence[i:i+length],'left':sequence[i-24:i],'right':sequence[i+length:],'pos':i})
        
    return substrings
def main():    
    # load guide library fasta file
    duplicated_guides={}
    nonduplicated_guides={}
    counts_read={}
    LIBRARY=SeqIO.parse(open(db),"fasta")
    for fasta in LIBRARY:
        counts_read[fasta.seq]={'read_count':0,'fwd':0,'rev':0,'pos':[]}
        if fasta.seq not in nonduplicated_guides.keys() and fasta.seq not in duplicated_guides.keys() :
            nonduplicated_guides.update({fasta.seq:fasta.id})
        else:
            
            if fasta.seq not in duplicated_guides.keys():
                duplicated_guides[fasta.seq]=[]
                duplicated_guides[fasta.seq].append(nonduplicated_guides[fasta.seq])
                del nonduplicated_guides[fasta.seq]
            duplicated_guides[fasta.seq].append(fasta.id)
    print(len(nonduplicated_guides.keys()),len(duplicated_guides.keys()))
    
    logging.info("Number of nonduplicated guides: %s, Number of duplicated guides: %s"%(len(nonduplicated_guides.keys()),len(duplicated_guides.keys())))
    ###read fastq files
    READ1=open(read1,'r')
    READ1=READ1.readlines()
    logging.info("Number of Reads: %s"%(int(len(READ1)/4)))
    
    for fasta in range(0,len(READ1),4):
        sequence=str(READ1[fasta+1].replace("\n",""))
        for i in range(20,min(28,len(sequence)-l)):
            substring={'seq':sequence[i:i+l],'left':sequence[i-20:i],'right':sequence[i+l:],'pos':i}
            if substring['left'][-5:]== "TCTCC"  and substring['right'][:6]=="GTTTTA"[:min(len(substring['right']),6)]:
                if substring['seq'] in duplicated_guides.keys() or substring['seq'] in nonduplicated_guides.keys():# only count perfect match
                    counts_read[substring['seq']]['read_count']+=1
                    counts_read[substring['seq']]['pos'].append(substring['pos'])
                    if substring['left'][-5:]=="TCTCC":
                        counts_read[substring['seq']]['fwd']+=1
                    if substring['right'][:6]=="GTTTTA"[:min(len(substring['right']),6)]:
                        counts_read[substring['seq']]['rev']+=1
                    break
        
            
    with open(output_file_name+"/counts.csv","w") as output:
        header=["sequence_20nt","guide_name","essentiality","NC",output_file_name.split("/")[-1],'fwd','rev','pos']
        dict_writer = csv.DictWriter(output, header, lineterminator='\n',delimiter='\t') 
        dict_writer.writeheader()
        for key in nonduplicated_guides.keys():
            if nonduplicated_guides[key][:6]=='random':
                NC=1
                essentiality=0
            else:
                NC=0
                if nonduplicated_guides[key].split("_")[1] in Keio_essential_genes:
                    essentiality=1
                else:
                    essentiality=0
            dict_writer.writerow({"sequence_20nt":key,"guide_name":nonduplicated_guides[key],"essentiality":essentiality,"NC":NC,output_file_name.split("/")[-1]:counts_read[key]['read_count'],\
                                  'fwd':counts_read[key]['fwd'],'rev':counts_read[key]['rev'],'pos':list(set(counts_read[key]['pos']))})
        for key in duplicated_guides.keys():
            if duplicated_guides[key][:6]=='random':
                NC=1
                essentiality=0
            else:
                NC=0
                for gene in duplicated_guides[key]:
                    if gene.split("_")[1] in Keio_essential_genes:
                        essentiality=1
                        break
                    else:
                        essentiality=0
            dict_writer.writerow({"sequence_20nt":key,"guide_name":",".join(duplicated_guides[key]),"essentiality":essentiality,"NC":NC,output_file_name.split("/")[-1]:counts_read[key]['read_count'],\
                                  'fwd':counts_read[key]['fwd'],'rev':counts_read[key]['rev'],'pos':list(set(counts_read[key]['pos']))})

if __name__ == '__main__':
    logging_file= output_file_name + '/log.txt'
    logging.basicConfig(filename=logging_file,format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info('db file: %s, Read1 file: %s' %(db,read1))
    main()
    logging.info("Execution Time: %s seconds" %('{:.2f}'.format(time.time()-start_time)))    