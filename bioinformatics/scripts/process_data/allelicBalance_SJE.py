#! /usr/bin/env python

# Author: Joana
# Modified by Simon Jacobsen Ellerstrand
# Original script found at https://github.com/joanam/scripts/blob/master/allelicBalance.py

from sys import *
import os, time, argparse, re, scipy, numpy
from collections import defaultdict
from scipy import stats

parser = argparse.ArgumentParser(description='Filter out genotypes with strong allic disbalance in the given VCF file')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-r1', '--ratiohom', type=float, dest='r1', help="hard cutoff for allelic ratio to recode genotype to homozygous [default: 0]", default=0)
parser.add_argument('-r2', '--ratioexcl', type=float, dest='r2', help="hard cutoff for allelic ratio to exclude genotype. Must be higher than -r1 if specified [default: 0]", default=0)

args = parser.parse_args()

readbytes=0
ratiohom=args.r1
ratioexcl=args.r2

if ratioexcl>0 and ratioexcl<=ratiohom:
	import sys
	sys.exit('r2 must be higher than r1')

inputF=open(args.i,'r')
outputF=open(args.o, 'w')


for Line in inputF:
	# DATA SECTION: clause checks if the header section is over
	if re.match('^#',Line) is None:

		# Get the columns of that line
		columns=Line.strip("\n").split("\t")

		# Only check SNPs
		if columns[5]!=".":
			# Add the info to the site
			result=columns[0:9]

			# Get the genotypes
			genotypecolumns=range(9,len(columns))

			# Check each individual if it is a heterozygote
			for ind in genotypecolumns:
				genotype=columns[ind]
				genotype=genotype.split(":")

				if "/" in genotype[0]:
					alleles=genotype[0].split("/") # counts the occurrences of each allele (up to 3 alt alleles -> LIMITATION)
				elif "|" in genotype[0]:
					alleles=genotype[0].split("|") # counts the occurrences of each allele (up to 3 alt alleles -> LIMITATION)

				# If the genotype is heterozygous check the allelic balance
				if alleles[0]!=alleles[1]:
					reads=genotype[1].split(",")

					# if one of the alleles has no reads (weirdly this happens)
					if reads[0]=="." or reads[1]=="." or int(reads[0])==0 or int(reads[1])==0:
						result.append("./.")

					else:
						# replace the genotype by the more common allele if the allelic ratio is very small
						if float(reads[0])/float(reads[1])<ratiohom or float(reads[1])/float(reads[0])<ratiohom:

							if reads[0]>reads[1]:
								genotype[0]=alleles[0]+"/"+alleles[0]
							else:
								genotype[0]=alleles[1]+"/"+alleles[1]

							result.append(":".join(genotype))


						# exclude if the allelic disbalance is strong
						elif float(reads[0])/float(reads[1])<ratioexcl or float(reads[1])/float(reads[0])<ratioexcl:

							result.append("./.")


						# If there is no allelic disbalance
						else:
							result.append(":".join(genotype))


				# If the genotype is homozygous, just append as is
				else:
					result.append(":".join(genotype))

			outputF.write('\t'.join(result)+"\n")

		# If it is not a SNP just write the line out
		else:
			outputF.write(Line)

	# If it is a header line, just write it out
	else:
		outputF.write(Line)



inputF.close()
outputF.close()
