#########################################
#										#
# Author: Kar-Tong Tan 					#
# Email: ktan@g.harvard.edu 			#
#										#
# #######################################

# Version 1.4
# - Made base cleanup return original case instead of upper case characters


# Version 1.5
# - gapped alignments could cause wrong base quality to be used. This version
#   deals with this issue properly.

# Version 1.5.1
# - Add in information about the base counts

# Version 1.5.2
# - Allow score cutoff for trinucleotide vs. binucleotide case

# Version 1.5.3
# - Ignore positions where the reference nucleotide is a 'N'
#   which can raise an error.

# Version 1.5.3
# - Ignore positions where the reference nucleotide is a 'N'
#   which can raise an error.

# Version 1.5.4
# - Reports the number of types of nucleotides detected
# - Also reports the monoallelic logP value
# - Rearrange and made main calling method a function


# Version 1.7
# - Parallel processing integrated in
# - Argument parsing 


import sys
import os
import subprocess
#sys.path.insert(0, "/home/kartong/Code/functions/python/NGS/")
#sys.path.insert(0, "/home/users/nus/csitkt/Code/functions/python/NGS/")
from lib.processMpileup import *
from lib.stats import *
from lib.calculateVariantStateLikelihoodModels_v1_1 import *
import contextlib
import argparse
import multiprocessing
import traceback
import time


# Initialize some global variables with some default values
MINBASEQUAL = 33 + 5
BASESTOTRIMFROMEDGE = 6
score_cutoff = -5
REPORT_BASE_COUNT = 1



def runProcess(exe):
	#p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	for line in p.stdout:
		yield line
	
	# This is necessary to tell OS to wait for exit status or we will
	# get a zombie process
	os.wait()



def excludeElementsFromList(inputList, indexToExclude):
	filteredList = []
	for index in xrange(len(inputList)):
		if index not in indexToExclude:
			filteredList.append(inputList[index])

	return filteredList


def processMpileupLine(pileupLine):
	'''
	Process the mpileup line
	'''
	lineArr   = pileupLine.strip().split("\t")
	position  = lineArr[1]
	refNt 	  = lineArr[2]
	interpretPileupResult = InterpretPileupStringMod2(lineArr[4])
	qualChars = lineArr[5]
	baseReadPositions = lineArr[6].split(',')

	# Get the base characters and the indices of sites which are 
	# gapped alignments
	baseChars = interpretPileupResult[0]
	indexOfGappedAlignment = interpretPileupResult[2]

	# Remove the base qualities and read positions for gapped alignment (i.e. > and <) sites if 
	# there are such sites (i.e. the list for the index is not empty)
	qualChars_excluded = []
	baseReadPositions_excluded = []

	if indexOfGappedAlignment:
		qualChars_excluded = excludeElementsFromList(qualChars, indexOfGappedAlignment)
		baseReadPositions_excluded = excludeElementsFromList(baseReadPositions, indexOfGappedAlignment)
	else:
		qualChars_excluded = qualChars
		baseReadPositions_excluded = baseReadPositions


	if len(qualChars_excluded) != len(baseChars) :
		sys.exit("Error: number of quality characters not equal to base characters")

	
	return [position, refNt, baseChars, qualChars_excluded, baseReadPositions_excluded]



# def cleanUpBaseArr(baseArr, qualArr, refNt):
# 	'''
# 	Convert all the base characters to ATGC format
# 	and remove all unusual an irrelevant characters.
# 	Also, convert all characters to uppercase
# 	'''
# 	baseArrNew = []
# 	qualArrNew = []
# 	refNtUpper = refNt.upper()


# 	for i in range(0, len(baseArr)):
# 		baseCurrent = baseArr[i]
# 		qualCurrent = qualArr[i]

# 		if ord(qualCurrent) < int(MINBASEQUAL):
# 			# Skip if base qual is too low
# 			continue

# 		if str(baseCurrent) == '.' or str(baseCurrent) == ',':
# 			baseArrNew.append(refNtUpper)
# 			qualArrNew.append(qualCurrent)

# 		elif str(baseCurrent) == 'N' or str(baseCurrent) == 'n':
# 			continue
		
# 		elif str(baseCurrent) == '*':
# 			continue

# 		else:
# 			baseArrNew.append(baseCurrent.upper())
# 			qualArrNew.append(qualCurrent)

# 	return(baseArrNew, qualArrNew, refNtUpper)





def cleanUpBaseArr(baseArr, qualArr, refNt, baseReadPosition):
	'''
	Convert all the base characters to ATGC format
	and remove all unusual an irrelevant characters.
	Also, all characters reported based as upper/lower
	case based on strand. Modified to also take in 
	read position.
	'''
	baseArrNew = []
	qualArrNew = []
	refNtUpper = refNt.upper()
	refNtLower = refNt.lower()
	baseReadPositionNew = []
	delCount 	= 0


	for i in range(0, len(baseArr)):
		baseCurrent = baseArr[i]
		qualCurrent = qualArr[i]
		baseReadPositionCurrent = baseReadPosition[i]

		if str(baseCurrent) == '*':
			'''
			Have to put this step before base qual filter
			as del ususally represented by low base qual.
			'''
			delCount += 1
			continue
		
		if ord(qualCurrent) < int(MINBASEQUAL):
			# Skip if base qual is too low
			continue

		if str(baseCurrent) == '.':
			baseArrNew.append(refNtUpper)
			qualArrNew.append(qualCurrent)
			baseReadPositionNew.append(baseReadPositionCurrent)

		elif str(baseCurrent) == ',':
			baseArrNew.append(refNtLower)
			qualArrNew.append(qualCurrent)
			baseReadPositionNew.append(baseReadPositionCurrent)

		elif str(baseCurrent) == 'N' or str(baseCurrent) == 'n':
			continue

		else:
			baseArrNew.append(baseCurrent)
			qualArrNew.append(qualCurrent)
			baseReadPositionNew.append(baseReadPositionCurrent)

	return(baseArrNew, qualArrNew, refNtUpper, baseReadPositionNew, delCount)


def filterBasePositionOnReads(baseArr, qualArr, baseReadPosition,
							  readlength, lengthFromEdgeToTrim):
	'''
	Trims and removes data that is a certain number of
	bases from the edge of the read.
	'''
	leftEdgeToKeep 	= int(lengthFromEdgeToTrim)
	rightEdgeToKeep	= int(readlength) - int(lengthFromEdgeToTrim)

	baseArrNew = []
	qualArrNew = []
	baseReadPositionNew = []

	for i in range(0, len(baseReadPosition)):
		currentBasePosn = baseReadPosition[i]
		
		if(int(currentBasePosn) >= int(leftEdgeToKeep) and int(currentBasePosn) <= int(rightEdgeToKeep)):
			baseArrNew.append(baseArr[i])
			qualArrNew.append(qualArr[i])
			baseReadPositionNew.append(currentBasePosn)

	return(baseArrNew, qualArrNew, baseReadPositionNew)




def calculateRDDvalue(baseChars, refNt):
	'''
	Calculate the RDD value of the site
	'''
	totalNtCount = len(baseChars)
	refNtCount = 0

	for nt in baseChars:
		if str(nt) == str(refNt):
			refNtCount += 1
	
	try:
		RDDval = float(1) - float(refNtCount) / totalNtCount
	except:
		RDDval = 0

	return RDDval


def checkStrandBias(baseChars, refNt):
	'''
	Check if there is strand bias in the non-reference 
	bases.
	'''
	refNtLower = refNt.lower()
	refUpperCount = 0
	refLowerCount = 0
	altUpperCount = 0
	altLowerCount = 0

	for character in baseChars:
		if character.isupper():
			if str(character) == str(refNt):
				refUpperCount += 1
			else:
				altUpperCount += 1
		else:
			if str(character) == str(refNtLower):
				refLowerCount += 1
			else:
				altLowerCount += 1

	refPlusProp = -1
	altPlusProp = -1

	try:
		refPlusProp = float(refUpperCount) / (refUpperCount + refLowerCount)
	except:
		pass

	try:
		altPlusProp = float(altUpperCount) / (altUpperCount + altLowerCount)
	except:
		pass


	return [refUpperCount, refLowerCount, altUpperCount, 
			altLowerCount, refPlusProp, altPlusProp]


def check_number_of_types_of_nucleotides(base_count, refNt):
	has_ref_nt = 0
	number_of_types_of_nt = 0

	nucleotide_types = ['A', 'T', 'G', 'C']
	if base_count[refNt.upper()] + base_count[refNt.lower()] > 0:
		has_ref_nt = 1

	for Nt in nucleotide_types:
		if(base_count[Nt.upper()] + base_count[Nt.lower()] > 0):
			number_of_types_of_nt += 1

	return(has_ref_nt, number_of_types_of_nt)




@contextlib.contextmanager
def smart_open(filename=None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()



def calculate_read_position_bias_scores(refNt, baseCharacters, baseReadPositions):
	number_of_base_characters = len(baseCharacters)
	number_of_base_positions = len(baseReadPositions)

	if(number_of_base_characters != number_of_base_positions):
		sys.exit("Error: number of base characters not equal to base positions")

	variant_nt_positions_fwd = []
	variant_nt_positions_rev = []
	for i in range(0, number_of_base_characters):
		currBaseChar = baseCharacters[i]
		currBasePosition = baseReadPositions[i]

		if(currBaseChar.upper() != refNt.upper()):
			if(currBaseChar.isupper()):
				variant_nt_positions_fwd.append(currBasePosition)
			else:
				variant_nt_positions_rev.append(currBasePosition)


	variant_readposn_median_fwd = "NA"
	variant_readposn_median_rev = "NA"
	median_abs_dev_fwd = "NA"
	median_abs_dev_rev = "NA"

	if len(variant_nt_positions_fwd) != 0:
		variant_readposn_median_fwd = calculate_median(variant_nt_positions_fwd)
		median_abs_dev_fwd = calculate_median_absolute_deviation(variant_readposn_median_fwd, variant_nt_positions_fwd)

	if len(variant_nt_positions_rev) != 0:
		variant_readposn_median_rev = calculate_median(variant_nt_positions_rev)
		median_abs_dev_rev = calculate_median_absolute_deviation(variant_readposn_median_rev, variant_nt_positions_rev)


	return(variant_readposn_median_fwd, variant_readposn_median_rev, median_abs_dev_fwd, median_abs_dev_rev)








def calculate_modification_score_for_region(bamFile, refGenome, chrom,	posnStart, posnEnd,	RDDpropFilter,
		READLENGTH,	deletionCountCutoff, minDepth, outputFile = None):


	posnString = str(chrom) + ':' + str(posnStart) + '-' + str(posnEnd)
	pileupCmd = "samtools mpileup -f %s -Q 0 -d 1000000000 -O -r %s %s" %(refGenome, posnString, bamFile)
	#print pileupCmd
	#print pileupCmd.split()

	if outputFile:
		sys.stdout = open(outputFile, 'w')

	# Iterate over all lines in pileupFile
	for mpileupLine in runProcess(pileupCmd.split()):
		# Extract the baseArr, baseQual and refNt
		[position, refNt, baseChars, qualChars, baseReadPositions] = processMpileupLine(mpileupLine)
		#print str(baseReadPositions) + "crp"
		[baseCharsCleanMixcase, qualCharsClean, refNtUpper, baseReadPositionsClean, deletionCount] = cleanUpBaseArr(baseChars, qualChars, refNt, baseReadPositions)

		# Ignore positions where the reference nucleotide is 'N' or other types like 'M'
		# which causes error in rest of the script
		allowed_Nt = ['A', 'T', 'G', 'C']
		#if str(refNtUpper) == str("N"):
		if str(refNtUpper) not in allowed_Nt:
			continue

		# Check if we want to trim bases
		if int(BASESTOTRIMFROMEDGE) > 0:
			[baseCharsCleanMixcase, qualCharsClean, baseReadPositionsClean] = \
			filterBasePositionOnReads(baseCharsCleanMixcase, qualCharsClean, 
								 	  baseReadPositionsClean, READLENGTH, 
								 	  BASESTOTRIMFROMEDGE)

		# Calc strand bias
		strandBiasResults = checkStrandBias(baseCharsCleanMixcase, refNtUpper)
		


		baseCharsClean = [x.upper() for x in baseCharsCleanMixcase]

		# Calculate RDD value
		RDDprop = calculateRDDvalue(baseCharsClean, refNtUpper)

		# Get Depth and deletion proportion
		depth = len(baseCharsClean) + int(deletionCount)
		try:
			deletionProp = float(deletionCount) / depth
		except:
			deletionProp = 'NA'

		# Do log-likelihood calculations of the 3 models
		[monoalleleLogP, bialleleLogP, trialleleLogP, quartalleleLogP] = calculateLikelihoodForDifferentModels(baseCharsClean, qualCharsClean, refNtUpper)

		# Calculate score difference between tri and bi case
		score_diff_tri_vs_bi = trialleleLogP - bialleleLogP

		# Generate results
		result = [chrom, position, refNt, monoalleleLogP, bialleleLogP, trialleleLogP, quartalleleLogP, RDDprop, depth, deletionCount, deletionProp] + strandBiasResults

		# Report the base counts as well	
		if REPORT_BASE_COUNT:
			base_Counts = {'A':0, 'T':0, 'G':0, 'C':0, 'a':0, 't':0, 'g':0, 'c':0}
			
			for base in baseCharsCleanMixcase:
				base_Counts[base] += 1


			[has_ref_nt , types_of_nt] = check_number_of_types_of_nucleotides(base_Counts, refNt)
			baseCountResults = [base_Counts['A'], base_Counts['T'], base_Counts['G'], base_Counts['C'], \
						base_Counts['a'], base_Counts['t'], base_Counts['g'], base_Counts['c']]
			
			result = result + baseCountResults + [has_ref_nt , types_of_nt]

		result.append(score_diff_tri_vs_bi)


	
		#print score_diff_tri_vs_bi
			
		if float(score_diff_tri_vs_bi) >= float(score_cutoff) and float(RDDprop) >= float(RDDpropFilter) and int(deletionCount) >= int(deletionCountCutoff) and depth >= 5:

			# Do some extra calcs if meet cutoff
			[variant_readposn_median_fwd, variant_readposn_median_rev, median_abs_dev_fwd, median_abs_dev_rev] = calculate_read_position_bias_scores(refNt, baseCharsCleanMixcase, baseReadPositionsClean)
			result = result + [variant_readposn_median_fwd, variant_readposn_median_rev, median_abs_dev_fwd, median_abs_dev_rev]


			print "\t".join(map(str, result))
			#yield "\t".join(map(str, result))
	
	if outputFile:
		sys.stdout.flush()
		#sys.stdout.close()
		#pass
	
	#time.sleep(2)
	pass


def concatenate_files(fileList, outputFile):
	'''
	Concatenate multiple small files into one big files.
	Preferred when there are too many files to concatenate
	in a single command.
	'''
	os.system("touch %s" %outputFile)
	#touchCmd = "touch %s" %outputFile
	#subprocess.call(touchCmd.split())
	for fileCurr in fileList:
		fileCurr = fileCurr.replace("|", "\|") # deal with file names with pipe char
		os.system("cat %s >> %s" %(fileCurr, outputFile))



def calculate_modification_score_for_region_parallel(data):
	[bamFile, refGenome, chrom,	posnStart, posnEnd,	RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth, outputFile] = data
	#print data
	try:
		starttime = time.time()
		calculate_modification_score_for_region(bamFile, refGenome, chrom,	posnStart, posnEnd,	RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth, outputFile = outputFile)
		endtime = time.time()
		startStr = 'Start\t'+ str(starttime) + '\t' + '\t'.join([chrom, posnStart, posnEnd])
		endStr = 'End\t' + str(endtime) + '\t' + '\t'.join([chrom, posnStart, posnEnd])
		os.system("echo %s >> outputFile" %startStr)
		os.system("echo %s >> outputFile" %endStr)
		#print startStr
		#print endStr
	except:
		print data
		print("Exception in worker:")
		traceback.print_exc()
		raise
		sys.exit()


def calculate_modification_score_parallel_wrapper(label, regionsFile, bamFile, refGenome, 
	number_of_processes, RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth):


	# label 		= sys.argv[1]
	# regionsFile = sys.argv[2]
	# bamFile 	= sys.argv[3]
	# refGenome 	= sys.argv[4]
	# number_of_processes = int(sys.argv[5])
	# RDDpropFilter 		= float(sys.argv[6])
	# READLENGTH 			= int(sys.argv[7])
	# deletionCountCutoff = int(sys.argv[8])

	final_combined_file = label + ".combined.txt"
	tmp_dir	 	= label + "_tmp/"


	try:
		os.mkdir(tmp_dir)
	except:
		print >> sys.stderr, '[Warning] tmp folder path exists'
		pass

	multiProc = multiprocessing.Pool(number_of_processes)

	regions = open(regionsFile, 'r')
	multiproc_input = []
	outputFile_list = []
	for region_line in regions:
		#print region_line.strip
		region_line_list = region_line.strip().split("\t")
		chrom = region_line_list[0]
		posnStart = region_line_list[1]
		posnEnd = region_line_list[2]

		outputFile = tmp_dir + chrom + ":" + posnStart + "-" + posnEnd
		outputFile_list.append(outputFile)
		input_params = [bamFile, refGenome, chrom,	posnStart, posnEnd,	RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth, outputFile]

		multiproc_input.append(input_params)

	regions.close()



	# Perform the multi processing
	# Wrap in try and finally because of zombine processes
	# (https://stackoverflow.com/questions/30506489/python-multiprocessing-leading-to-many-zombie-processes)
	try:
		multiProc.map(calculate_modification_score_for_region_parallel, multiproc_input)
		multiProc.close()
		multiProc.join()
	except:
		multiProc.close()
		multiProc.join()
		raise
	
	#concatenate_args = ['cat'] + outputFile_list + ['>', final_combined_file]
	#print concatenate_args
	#concatenate_cmd = " ".join(concatenate_args)
	#os.system(concatenate_cmd)
	concatenate_files(outputFile_list, final_combined_file)


def printHeader():
	header = ["chrom", "position", "reference_nt", "mono-alleleLogP", "bi-alleleLogP", "tri-alleleLogP", "tetra-alleleLogP", "variant_proportion", "depth", "deletionCount", "deletion_proportion", 
"refUpperCount", "refLowerCount", "altUpperCount", "altLowerCount", "refPlusProp", "altPlusProp",
"A_count", "T_count", "G_count", "C_count", "a_count", "t_count", "g_count", "c_count", "has_reference_nt", "types_of_nt", "ModTect_score", "variant_readposn_median_fwd", "variant_readposn_median_rev", "median_abs_dev_fwd", "median_abs_dev_rev"]
	
	print "\t".join(header)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(prog ='prog', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("bamfile", type=str, help="bamfile")
	parser.add_argument("genome", type=str, help="reference genome file")
	parser.add_argument("chrom", type=str, help="chomosome")
	parser.add_argument("posnStart", type=str, help="start position, 1 based")
	parser.add_argument("posnEnd", type=str, help="end position, 1 based")
	#parser.add_argument('--minBaseQual', nargs='?', const=1, type=int, default=5, help='minimum base quality required for analysis')
	parser.add_argument('--minBaseQual', type=int, default=5, help='minimum base quality required')
	parser.add_argument('--basesToTrimFromEdge', type=int, default=6, help='bases to trim from edge of read')
	parser.add_argument('--scoreCutoff', type=float, default=1, help='modification score cutoff')
	parser.add_argument('--minDepth', type=int, default=5, help='min sequencing depth required')
	parser.add_argument('--readlength', type=float, default=101, help='length of read for input data')
	parser.add_argument('--mismatchPropFilter', type=float, default=0, help='minimum mismatch proportion required')
	parser.add_argument('--deletionCountCutoff', type=int, default=0, help='minimum number of deletions required')
	parser.add_argument('--threads', type=int, default=1, help='number of parallel processes to run at one time')
	parser.add_argument('--regionFile', type=str, help='A file containing the regions to analyze. Causes positional arguments chrom, posnStart and posnEnd to be ignored')
	parser.add_argument('--label', type=str, default="testingMod", help='The label for the outputfile. Required if regionFile is specified')
	args = parser.parse_args()


	MINBASEQUAL = 33 + args.minBaseQual
	BASESTOTRIMFROMEDGE = args.basesToTrimFromEdge
	score_cutoff = args.scoreCutoff
	REPORT_BASE_COUNT = 1


	# bamFile   = sys.argv[1]
	# refGenome = sys.argv[2]
	# chrom     = sys.argv[3]
	# posnStart = sys.argv[4]
	# posnEnd   = sys.argv[5]
	# RDDpropFilter = sys.argv[6]
	# READLENGTH = sys.argv[7]
	# deletionCountCutoff = sys.argv[8]


	bamFile   = args.bamfile
	refGenome = args.genome
	chrom     = args.chrom
	posnStart = args.posnStart
	posnEnd   = args.posnEnd
	RDDpropFilter = args.mismatchPropFilter
	READLENGTH = args.readlength
	deletionCountCutoff = args.deletionCountCutoff
	minDepth = args.minDepth
	label = args.label
	regionsFile = args.regionFile
	number_of_processes = args.threads

	if(args.regionFile):
		'''
		Performs parallel analysis of the region file
		'''

		calculate_modification_score_parallel_wrapper(label, regionsFile, bamFile, refGenome, 
		number_of_processes, RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth)

	else:
		'''
		Performs analysis based on the defined region
		'''
		printHeader()
		calculate_modification_score_for_region(bamFile, refGenome, chrom,	posnStart, posnEnd,	RDDpropFilter, READLENGTH, deletionCountCutoff, minDepth)
	
