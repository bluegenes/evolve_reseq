# equalize references (dumb it down)
#!/usr/bin/env python
#######################################
### Tessa Pierce ###
### 11/18/2016 ###
#######################################

import sys
import math
import argparse

#read one fasta into memory
def import_fasta(ref1_fasta):
    fastaDict = {}
    bp = ''
    with open(ref1_fasta, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if len(bp) > 0:
                    fastaDict[name] = bp
                    bp = ''
                name = line
            else:
                bp = bp + line
        fastaDict[name] = bp
    return fastaDict

def write_contig(name,bp,outF):
    outF.write(name + '\n' + bp + '\n')

def findNs(bases, ch='N'):
    return [i for i, ltr in enumerate(bases.upper()) if ltr == ch]

def equalize_contigs(name,bp1,bp2,out1,out2):
    shortest = min(len(bp1),len(bp2))
    bp1 =bp1[:shortest] # truncate to shortest
    bp2 =bp2[:shortest]
    all_Ns = set(findNs(bp1) + findNs(bp2)) # get all indices with N's
    bp1_list=list(bp1)
    bp2_list=list(bp2) # make list, so is mutable
    for i in all_Ns:
        bp1_list[i] = 'n' #change all N indices to 'n' in both bp lists
        bp2_list[i] = 'n'
    eq_bp1 = ''.join(bp1_list) # change list back to string
    eq_bp2 = ''.join(bp2_list)
    write_contig(name,eq_bp1,out1) # write out new contigs in appropriate files
    write_contig(name,eq_bp2,out2)

def equalize_refs(ref1,ref2,o1,o2):
    ref1Dict = import_fasta(ref1)
    out1 = open(o1, 'w')
    out2 = open(o2, 'w')
    with open(ref2, 'r') as f:
        bp=''
        for line in f:
            line=line.rstrip()
            if line.startswith('>'):
                if len(bp) > 0:
                    if name in ref1Dict.keys():
                        ref1contig = ref1Dict[name]
                        equalize_contigs(name, ref1contig, bp, out1, out2)
                    bp = ''
                name = line
            else:
                bp = bp + line
        if name in ref1Dict.keys():
            ref1contig = ref1Dict[name]
            equalize_contigs(name, ref1contig, bp, out1, out2)
    out1.close()
    out2.close()


if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="Equalize References")
    parser.add_argument('-i1', '--firstRef', help='first reference fasta')
    parser.add_argument('-i2', '--secondRef', help='second reference fasta')
    parser.add_argument('-o1', '--out1', help='output equalized reference 1')
    parser.add_argument('-o2', '--out2', help='output equalized reference 2')
    args = parser.parse_args()
    equalize_refs(args.firstRef,args.secondRef,args.out1,args.out2)
