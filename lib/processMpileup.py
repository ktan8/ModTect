from mpileup_mod import *


def readFile(filename):
	f = open(filename, 'r')
        line      = f.read()
	return line


def processMpileupFile(filename, outputFile):
	pileupLine = readFile(filename)
	processMpileupLine(pileupLine, outputFile)


def processMpileupLine(pileupLine, outputFile, outputFormat='byCol'):
	'''
	byCol puts per read information into a column
	'''
	OFFSET = 33
	
	# Process the mpileup line
	lineArr   = pileupLine.strip().split("\t")
	baseChars = InterpretPileupStringMod(lineArr[4])[0]

	
	# Proces the base qual characters
	baseQual = lineArr[5]
	baseQualVal = []
	for char in baseQual:
		baseVal = ord(char) - OFFSET
		baseQualVal.append(baseVal)

	
	# Process the read position
	readPosn = lineArr[6].split(',')
	

	# Write the data out as stand
	if str(outputFormat) == str('byCol'):
		f = open(outputFile, 'w')
		f.write('\t'.join(baseChars) + "\n")
		f.write('\t'.join(map(str, baseQualVal)) + "\n")
		f.write('\t'.join(readPosn) + "\n")
		f.close()
	elif str(outputFormat) == str('byRow'):
		f = open(outputFile, 'w')
		for i in range(0, len(baseChars)):
			resultLine = [baseChars[i], baseQualVal[i], readPosn[i]]
			f.write('\t'.join(map(str, resultLine)) + "\n")
		f.close()
	else:
		raise NameError('unknown outputFormat was indicated')
	
