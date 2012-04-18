#!/usr/bin/env python2.6

# Mariya Semenikhina
# FALL09 reverseORF

import  sys
import  optparse
import math
import gzip
import string
from string import maketrans
import random

def eachfastagz(fileobj):
    sofar = fileobj.readline()
    sofar = ''#sofar.upper()
    print sofar
    for line in fileobj:
        line = line.rstrip('\r\n')  # strip newline and MS-DOS ^M
        if line.startswith('>'):
            yield sofar
            sofar = ''#line
        else:
            sofar += line.upper()
    yield sofar

def get_counts(counts, seq, k):
    for start in range(len(seq)-k):
	if seq[start:start+k] not in counts:
	    counts[(seq[start:start+k])] = 1
	else:
	    counts[(seq[start:start+k])] += 1
    return counts
def reverse_comp (dna):    # define reverse complement
    complement_table = maketrans("ACGTU", "TGCAA")
    return dna[::-1].translate(complement_table)

def reverse_kmers(kmers, rev_kmers):
    comp_key = ''
    for key in kmers.keys():
	comp_key = reverse_comp(key)
	if comp_key in rev_kmers:
	    rev_kmers[comp_key] += int(kmers[key])
	else:
	    rev_kmers[comp_key] = int(kmers[key])
    return rev_kmers

def getORFs (ORF, seq, stop_codons):
    k=len(stop_codons[0])
    i = len(ORF.keys())
    ORF = ''#ORF[i]=''
    for start in range(len(seq)-k):
	if seq[start:start+k] not in stop_codons:
	    ORF+=(seq[start:start+k])#ORF[i]+=(seq[start:start+k])
	else:
	    i+=1
	    yield ORF
	    ORF = ''#ORF[i]=''
    yield ORF #return ORF

def maxORF(ORF, seq, stop_codons):
    k=len(stop_codons[0])
    ORFnew = ''
    for start in range(len(seq)-k):
        if seq[start:start+k] not in stop_codons:
            ORFnew+=(seq[start:start+k])
        else:
            ORF = max(ORFnew,ORF)
	    ORFnew = ''
    return ORF

def get_GC(seq, counts):
    for letter in seq:
	if letter in counts:
	    counts[letter]+= 1
	else:
	    counts[letter] = 1
    return counts
def get_codon(fileobj, codon, translate_tbl):
    for line in fileobj:
	line = line.rstrip('\r\n')  # strip newline and MS-DOS ^M
	codon[line.split()[0]] = line.split()[2]
	translate_tbl[line.split()[0]] =  (line.split()[1]).strip('()')
    return codon

def main():
    parser = optparse.OptionParser()
    parser.add_option('-c', '--codon', 
                  action="store", type="string",
		  dest="codon_filename", 
                  help="Read the codon bias table from the specified file."
                  )
    parser.add_option('-a', '--alphabet',
                  dest="alphabet",
                  default='ACGTU',
                  type = 'str'
                  )
    parser.add_option('-n','--num_sequences',
                  dest="num_sequences",
                  default = 10000,
                  type='int',
		  help="The default number of sequences generated in the sample should be 10 000."
                  )
    parser.add_option("-f", "--file",
                  type="string",
                  dest="filename",
                  help="FASTA FILE")
    parser.add_option("-p", "--protein",  
		  type="string", 
		  dest="protein",
                  help="Read the protein sequence for the forward strand (in FASTA format) from the specified file."
		  )
    options, remainder = parser.parse_args()

	#######Handles gzipped files default stdin#######
    if options.filename:
        if options.filename.endswith(".gz"):
            file = gzip.open(options.filename,'r')
        else:
            file = open(options.filename,'r')
    else:
	file = sys.stdin
    bases = 0
    GC = dict()
    ORF = dict()
    kmers = dict()
    rev_kmers = dict()
    stop_codons = ('TAG','TGA', 'TAA')

    for i in eachfastagz(file):
	bases+=len(i)
	get_GC(i, GC)
	ORF = getORFs (ORF, i, stop_codons)
	kmers = get_counts(kmers, i, 3)
	kmers = get_counts(kmers, reverse_comp(i),3)
    file.close()
    #print GC
    #rev_kmers = reverse_kmers(kmers, rev_kmers)
    sum_bases = 0
    sum_kmers = 0
    bases +=bases
    for value in GC.values():
	sum_bases += float(value)
    for value in kmers.values():
        sum_kmers += float(value)
    GC_freq = (float(GC['G'])+float(GC['C']))/sum_bases
    TA_freq = 1-GC_freq
    prob_ATG = ((float(GC['A'])/sum_bases)*(float(GC['T'])/sum_bases)*(float(GC['G'])/sum_bases))
    count_ATG = float(kmers['ATG'])/(sum_kmers)
    prob_stop = ((float(GC['T'])/sum_bases)*(float(GC['A'])/sum_bases)*(float(GC['G'])/sum_bases))*2
    prob_stop +=((float(GC['T'])/sum_bases)*(float(GC['A'])/sum_bases)*(float(GC['A'])/sum_bases))
    count_stop = float(kmers['TAG'])/sum_kmers+float(kmers['TGA'])/sum_kmers+float(kmers['TAA'])/sum_kmers
    expected_GC_orf = prob_ATG*(1-(prob_stop)/2)**387
    expected_kmer_orf = count_ATG*(1-(count_stop)/3)**387
        
    print'# Reading from %s'%options.filename
    print'# The double-stranded genome contains %d bases, and %d  3-mers.'%(bases, sum_kmers)
    print'# probablility of G+C is %f'%GC_freq
    print'# probability of ATG start coon is %f based on GC conteot, %f based on dipect count'%(prob_ATG, count_ATG)
    print'# probabilityofsstop coeon is %f based on GC conteot, %f based on diprect count'%(prob_stop,count_stop)
    print'# Expected number of 388 or longer orfs based on GC frequency is %.4g'%expected_GC_orf
    print'# Expected number of 388 or longer orfs based on 3-mer frequency is %.4g' %expected_kmer_orf

    codon = dict()
    translate_tbl = dict()
    if options.codon_filename:
        if options.codon_filename.endswith(".gz"):
            file = gzip.open(options.codon_filename,'r')
        else:
            file = open(options.codon_filename,'r')
	codon =  get_codon(file, codon, translate_tbl)

	prob = list()
	ORF_data = dict()
	ORF=dict()

	for key in codon:
	    for i in range(0, int(float(codon[key])*10)):
		prob.append(key)
	for i in range(1,options.num_sequences):
	    seq = 'ATG'
	    for j in range(1,559):
		seq+=prob[random.randint(0, len(prob)-1)]
	    #ORF = getORFs (ORF, seq, stop_codons)
	    #ORF = getORFs (ORF, seq[1::], stop_codons)
	    #ORF = getORFs (ORF, seq[2::], stop_codons)
	    for orfy in getORFs(ORF, seq, stop_codons): #for key in ORF.keys():
	    	temp = len(orfy) #temp = str(len(key))
	        if temp in ORF_data:
		    ORF_data[temp] += 1
	        else:
		    ORF_data[temp] = 1
	keys = ORF_data.keys()
	keys.sort()
	size_test =sum(ORF_data.values())
	for key in keys:
	    var =  float(ORF_data[key]) * (1-(float(ORF_data[key])/size_test))
	    sigma = math.sqrt(var)
	    p_value = (float(ORF_data[key])/size_test-key/560)/sigma
	    print '%-s \t %d \t %f \t %f\t' % (key, ORF_data[key], float(ORF_data[key])/size_test,p_value)
	#if options.protein:
	#    if options.protein.endswith(".gz"):
        #	file = gzip.open(options.protein,'r')
        #    else:
	#	file = open(options.protein,'r')
	    #for i in eachfastagz(file):
		#for aa in i

if __name__ == "__main__":
    main()
