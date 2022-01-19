import sys
import pysam
import getopt
import os.path
from Bio.Seq import Seq
import re
import argparse
import gzip
import string

## NOTE: this may require python 2.7.5...

AB_LEN = 8 # antibody barcode length. At the beginning of the R1 & R2.
TM_LEN = 42 # antibody + adaptor length removed from 5' end.

# Argument parsing:
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', type=str, required=True, metavar='<R1 file>',
						help='input R1 fastq file')
	parser.add_argument('-odir', type=str, required=True,
						help='output folder path')
	parser.add_argument('-b1', nargs='+', required=True,
						help='list of antibody barcode for read 1')
	parser.add_argument('-b2', nargs='+', required=True,
						help='list of antibody barcode for read 2')
	# Use like: python arg.py -l 1234 2345 3456 4567
	# nargs='+' takes 1 or more arguments, nargs='*' takes zero or more
	parser.add_argument('-b', type=str, required=True, metavar='<validBcFile',
						help='valid cell barcode list file')
	args = parser.parse_args()
	return args

def hammingDist(x, y):
	# find the Hamming distance between two input strings:
	if len(x)!=len(y):
		hd = len(x)
	else:
		hd = 0
		for i in range(len(x)):
			if x[i]!=y[i]:
				hd+=1	# count mismatches
	return hd

### Barcode error-correction:
def correctBC(bc, bcset):
	bcNew = None	# default (uncorrectable)
	for k in bcset:
		if hammingDist(bc, k)<=1:
			bcNew = k # corrected barcode
		else:
			bcNew = None # barcode can't be unambiguously corrected
	return bcNew # return corrected barcode

def parseBarcodes(bc1, bc2, bc, bc1set, bc2set, bcset): # Need to add Cell barcode later.
	# returned flags
	bc1Valid = 0
	bc2Valid = 0
	bcValid = 0
	bc1Corr = 0
	bc2Corr = 0
	bcCorr = 0
	bailout = 0

	# test bc1
	if bc1 in bc1set:
		bc1Valid = 1
	else:
		bc1 = correctBC(bc1, bc1set) # corrected barcode or None if not correctable
		if bc1 != None:
			bc1Corr = 1
		else:
			bailout = 1 # exit without testing bc2 and bc

	# test bc2
	if not bailout:
		if bc2 in bc2set:
			bc2Valid = 1
		else:
			bc2 = correctBC(bc2, bc2set)
			if bc2 != None:
				bc2Corr = 1
			else:
				bailout = 1

	# test bc
	if not bailout:
		if bc in bcset:
			bcValid = 1
		else:
			bc = None

	return (bc1, bc2, bc, bc1Valid, bc2Valid, bcValid, bc1Corr, bc2Corr, bcCorr)

tab = string.maketrans("ACTG", "TGAC")
def reverse_complement(seq):
	return seq.translate(tab)[::-1]

def fastqWrite(f, r, rName, lTrim):
	# write each field of the read:
	f.write('@%s %s\n' % (rName, r.comment))	 # read name and comment
	f.write('%s\n' % r.sequence[lTrim:])				  # the sequence
	f.write('+%s %s\n' % (rName, r.comment))	 # read name and comment (filler?)
	f.write('%s\n' % r.quality[lTrim:])				   # the quality string
	return

def run(args):
	R1file = args.i # R1 fastq: read 1
	R2file = args.i.replace('R1','R2') # R2 fastq: index 2 10x cell barcode
	R3file = args.i.replace('R1','R3') # R3 fastq: read 2
	odir = args.odir # out put directory
	if odir==None:
		odir = os.path.dirname(R1file)  # if not provided, put the output file in the same directory as the input
	bc1set = args.b1 # antibody barcodes in R1
	bc2set = args.b2 # antibody barcodes in R2
	bcfile = open(args.b, 'r')
	bcset = bcfile.read().splitlines()

	# open the input file:
	fq1 = pysam.FastqFile(R1file)
	fq2 = pysam.FastqFile(R2file)
	fq3 = pysam.FastqFile(R3file)
	R1name = os.path.basename(R1file)
	ofile = os.path.join(odir, R1name.replace('.fastq.gz','_valid.fastq.gz'))
	fq1out = gzip.open(ofile,'wb')  # open outut file for read 1
	fq2out = gzip.open(ofile.replace('R1','R2'),'wb') # save read 2 'R3' to 'R2'.
	bc1outname = ofile.replace('_valid.fastq.gz','_bc1.txt')
	bc1out = open(bc1outname,'wb')
	bc2out = open(bc1outname.replace('R1','R2').replace('bc1','bc2'),'wb')

	# counters:
	nCount = 0 # total input reads number
	nValid = 0 # final valid counts (bc1 & bc2 & bc all valid)
	nbc1Valid = 0 # valid antibody barcodes in R1
	nbc2Valid = 0 # valid antibody barcodes in R3
	nbcValid = 0 # valid cell barcodes in R2
	nbc1Corr = 0 # Corrected antiboty barcodes in R1
	nbc2Corr = 0 # Corrected antiboty barcodes in R3
	nbcCorr = 0 # Corrected cell barcodes in R2

	# loop over all reads:
	while 1:
		try:
			r1 = fq1.next() # read 1 (P7 side)
			r2 = fq2.next() # 10x cell barcode
			r3 = fq3.next() # read 2 (P5 side)
			nCount+=1 # read counter
		except StopIteration:
			break # last item
		
		# get the two antibody barcodes
		bc1 = r1.sequence[:AB_LEN]
		bc2 = r3.sequence[:AB_LEN]
		bc = reverse_complement(r2.sequence)

        # write the original barcode before correcting
		bc1out.write('%s\n' % bc1)
		bc2out.write('%s\n' % bc2)

		# check three barcodes and update the counts:
		(bc1, bc2, bc, bc1Valid, bc2Valid, bcValid, bc1Corr, bc2Corr, bcCorr) = parseBarcodes(bc1, bc2, bc, bc1set, bc2set, bcset)
		nbc1Valid += bc1Valid
		nbc2Valid += bc2Valid
		nbcValid += bcValid
		nbc1Corr += bc1Corr
		nbc2Corr += bc2Corr
		nbcCorr += bcCorr

		# write out the sequence read if bc1, bc2 and bc are all valid:
		if bc1!=None and bc2!=None and bc!=None:
			nValid += 1
			# create the new read name
			rName1 = '%s:%s:%s:%s' % (r1.name, bc, bc1, bc2) # paired reads require the same name.
			rName3 = '%s:%s:%s:%s' % (r3.name, bc, bc1, bc2) # always add the barcode -> cell:read1:read2
			fastqWrite(fq1out, r1, rName1, TM_LEN)
			fastqWrite(fq2out, r3, rName3, TM_LEN)

	# close the input files:
	fq1.close()
	fq2.close()
	fq3.close()
	# close the output files:
	fq1out.close()
	fq2out.close()
	bc1out.close()
	bc2out.close()

	# print counts:
	print 'total input reads: %d' % nCount
	print 'final valid reads: %d' % nValid
	print 'valid antibody barcodes in R1: %d' % nbc1Valid
	print 'valid antibody barcodes in R3: %d' % nbc2Valid
	print 'valid cell barcodes in R2: %d' % nbcValid
	print 'Corrected antiboty barcodes in R1: %d' % nbc1Corr
	print 'Corrected antiboty barcodes in R3: %d' % nbc2Corr
	print 'Corrected cell barcodes in R2: %d' % nbcCorr

	return

if __name__ == "__main__":
	args = get_args()
	run(args)











