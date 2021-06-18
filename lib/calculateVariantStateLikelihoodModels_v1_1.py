from collections import Counter
from math import log10

BASEOFFSET = 33
MINBASEQUAL = 5

def quartAllelicCase(baseArr, qualArr, refNt):
	# it must be able to tell if the nucleotide is m1, m2 or m3
	# and then which epsilon to use for computation
	# inputs: two arrays of base seq and base qual and the reference base


	(varNtOrder, varNtProportion) = identifyMutFractionAndOrdering(baseArr, refNt)

	# Set default values of eps and varNT in case there are less
	# than 3 variant nucleotides in the data
	eps1 = 0
	eps2 = 0
	eps3 = 0
	varNt1 = 'N'
	varNt2 = 'N'
	varNt3 = 'N'

	try:
		'''
		Need to put the 6 variables in the right order.
		if varNtOrder[2] is not defined, all 3 eps values
		would not be defined causing score difference
		to give a negative value
		'''
		#varNt1 = varNtOrder[0]
		#varNt2 = varNtOrder[1]
		#varNt3 = varNtOrder[2]
		#eps1 = varNtProportion[0]
		#eps2 = varNtProportion[1]
		#eps3 = varNtProportion[2]

		varNt1 = varNtOrder[0]
		eps1 = varNtProportion[0]

		varNt2 = varNtOrder[1]
		eps2 = varNtProportion[1]

		varNt3 = varNtOrder[2]
		eps3 = varNtProportion[2]

	except:
		pass


	sumLogLikelihood = 0

	for i in range(0,len(baseArr)):
		currentBase = baseArr[i]
		currentQual = qualArr[i]
		actualQualVal = ord(currentQual) - BASEOFFSET
		error = 10 ** (-1 * (float(actualQualVal) / 10) )
		prob = None

		if str(currentBase) == str(refNt):
			prob = (eps1 + eps2 + eps3) * (error / 3) + (1 - eps1 - eps2 - eps3) * (1 - error)

		elif str(currentBase) == str(varNt1):
			prob = eps1 * (1 - error) + (1 - eps1) * (error / 3)

		elif str(currentBase) == str(varNt2):
			prob = eps2 * (1 - error) + (1 - eps2) * (error / 3)

		elif str(currentBase) == str(varNt3):
			prob = eps3 * (1 - error) + (1 - eps3) * (error / 3)

		else:
			print "shit" + str(currentBase)
			raise Exception("Variant nucleotides in none of the cases %s" %currentBase)

		#print prob
		logProb = log10(prob)
		sumLogLikelihood += logProb

	return sumLogLikelihood



	# Need 3 variant bases for this model. What if there
	# isn't 3 variant bases and only 1 or 2???
	# take care of this later. Let's just work on simple case first


def triAllelicCase(baseArr, qualArr, refNt):

	(varNtOrder, varNtProportion) = identifyMutFractionAndOrdering(baseArr, refNt)

	# Set default values of eps and varNT in case there are less
	# than 2 variant nucleotides in the data
	eps1 = 0
	eps2 = 0
	varNt1 = 'N'
	varNt2 = 'N'

	try:
		varNt1 = varNtOrder[0]
		eps1 = varNtProportion[0]
		varNt2 = varNtOrder[1]
		eps2 = varNtProportion[1]
	except:
		pass


	sumLogLikelihood = 0

	for i in range(0,len(baseArr)):
		currentBase = baseArr[i]
		currentQual = qualArr[i]
		actualQualVal = ord(currentQual) - BASEOFFSET
		error = 10 ** (-1 * (float(actualQualVal) / 10) )
		prob = None

		if str(currentBase) == str(refNt):
			prob = (eps1 + eps2) * (error / 3) + (1 - eps1 - eps2) * (1 - error)

		elif str(currentBase) == varNt1:
			prob = eps1 * (1 - error) + (1 - eps1) * (error / 3)

		elif str(currentBase) == varNt2:
			prob = eps2 * (1 - error) + (1 - eps2) * (error / 3)

		else:
			prob = error / 3

		logProb = log10(prob)
		sumLogLikelihood += logProb

	return sumLogLikelihood



def biAllelicCase(baseArr, qualArr, refNt):

	(varNtOrder, varNtProportion) = identifyMutFractionAndOrdering(baseArr, refNt)

	# Set default values of eps and varNT in case there are less
	# than 2 variant nucleotides in the data
	eps1 = 0
	varNt1 = 'N'

	try:
		varNt1 = varNtOrder[0]
		eps1 = varNtProportion[0]
	except:
		pass


	sumLogLikelihood = 0

	for i in range(0,len(baseArr)):
		currentBase = baseArr[i]
		currentQual = qualArr[i]
		actualQualVal = ord(currentQual) - BASEOFFSET
		error = 10 ** (-1 * (float(actualQualVal) / 10) )
		prob = None

		if str(currentBase) == str(refNt):
			prob = (eps1) * (error / 3) + (1 - eps1) * (1 - error)

		elif str(currentBase) == varNt1:
			prob = eps1 * (1 - error) + (1 - eps1) * (error / 3)

		else:
			prob = error / 3

		logProb = log10(prob)
		sumLogLikelihood += logProb

	return sumLogLikelihood


def monoAllelicCase(baseArr, qualArr, refNt):
	
	sumLogLikelihood = 0
	
	for i in range(0,len(baseArr)):
		currentBase = baseArr[i]
		currentQual = qualArr[i]
		actualQualVal = ord(currentQual) - BASEOFFSET
		error = 10 ** (-1 * (float(actualQualVal) / 10) )
		prob = None

		if str(currentBase) == str(refNt):
			prob = (1 - error)
		else:
			prob = error / 3

		logProb = log10(prob)
		sumLogLikelihood += logProb
	
	return sumLogLikelihood



def identifyMutFractionAndOrdering(baseArr, refNt):
	'''
	Identify which is the highest proportion variant bases and
	also calculate the fraction of each base. Input should all
	be given in uppercase characters
	'''
	baseCountEachBase = Counter(baseArr)
	totalDepth = sum(baseCountEachBase.values())
	varBaseProportion = []
	varBaseOrder = []
	
	for (baseChar, baseCount) in baseCountEachBase.most_common():
		if str(baseChar) != str(refNt):
			currBaseProportion = float(baseCount) / totalDepth
			varBaseProportion.append(currBaseProportion)
			varBaseOrder.append(baseChar)

	return(varBaseOrder, varBaseProportion)


def calculateLikelihoodForDifferentModels(baseArr, qualArr, refNt):
	monoalleleLogP  = monoAllelicCase(baseArr, qualArr, refNt)
	bialleleLogP 	= biAllelicCase(baseArr, qualArr, refNt)
	trialleleLogP 	= triAllelicCase(baseArr, qualArr, refNt)
	quartalleleLogP = quartAllelicCase(baseArr, qualArr, refNt)

	return (monoalleleLogP, bialleleLogP, trialleleLogP, quartalleleLogP)



if __name__ == "__main__":
	print identifyMutFractionAndOrdering(['A','A','A','A','T','T','T','G','G','C'], 'C')
	print identifyMutFractionAndOrdering(['A','A','A','A','T','T'], 'C')
	print identifyMutFractionAndOrdering(['A','A','A','A','T','T'], 'A')
	print identifyMutFractionAndOrdering(['A','A','A','A'], 'A')


	quartAllelicCase(['A','A','A','A','T','T','T','G','G','C'], ['I','I','I','I','I','I','I','I','I','I'], 'C')

	quartAllelicCase(['A','A','A','A','T','T','T','G','G'], ['I','I','I','I','I','I','I','I','I'], 'A')


	# We can see the following case reduces to the same results
	print "===================="
	print triAllelicCase(['A','A','A','A'], ['I','I','I','I'], 'A')
	print quartAllelicCase(['A','A','A','A'], ['I','I','I','I'], 'A')
	print "===================="

	# We can see in the following 4 nucleotide case, the quart
	# allelic model explains the data the best. about 10^7 better
	print biAllelicCase(['A','C','G','T'], ['I','I','I','I'], 'A') # -9.68127984305
	print triAllelicCase(['A','C','G','T'], ['I','I','I','I'], 'A') # -5.98230018697
	print quartAllelicCase(['A','C','G','T'], ['I','I','I','I'], 'A') # -2.40823996531





# It seems that if there is only two types of nucleotides, the quart case will
# reduce to the bi-case and give the same likelihood values as the bi-case
# even if that subroutine was used for calculation
