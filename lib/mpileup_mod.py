import operator
from collections import defaultdict

def InterpretPileupString(alnstring):
	'''
	Interprets a samtools mpileup string and converts it into
	a list of nts characters and another list of indels
	'''

	chariter = iter(alnstring)
	indels = []
	nts = []
	for char in chariter:
		# Ignore the character following the ^ sign that represents the
 		# mapping quality of the read
 		if char == '^':
 			chariter.next()
 			continue

 		# Deal with character indicating end of read
 		if char == '$':
 			continue
		
		# Ignore gapped alignment for RNA
		if char == '>' or char == '<':
			continue

		# Skip asterisk for deleted base
		if char == '*':
			continue


		# Check for the character flag indicating an indel
		if char == '+' or char == '-':
			indelsign = char
			indelsize = ""
			char = chariter.next()
			
			# Get the size of the indel
			while(char.isdigit() == True):
				indelsize = indelsize + char
				char = chariter.next()

			# Extract the seq of the indel
			indelstring = char # first character has been previously accessed
			for x in range(int(indelsize) - 1):
				indelstring = indelstring + chariter.next()

			# Generate the final indel string
			indelstring = indelsign + indelsize + indelstring

			indels.append(indelstring)			
			continue

		nts.append(char)
	
	return nts, indels




def InterpretPileupStringMod(alnstring):
		'''
		Interprets a samtools mpileup string and converts it into
		a list of nts characters and another list of indels. This
		is modified to also count asterisk as part of the result.
		'''

		chariter = iter(alnstring)
		indels = []
		characters = []
	
		for char in chariter:
			# Ignore the character following the ^ sign that represents the
			# mapping quality of the read
			if char == '^':
				chariter.next()
				continue

			# Deal with character indicating end of read
			if char == '$':
				continue

			# Ignore gapped alignment for RNA
			if char == '>' or char == '<':
				continue

			# Include asterisk for deleted base
			if char == '*':
				characters.append('*')
				continue


			# Check for the character flag indicating an indel
			if char == '+' or char == '-':
				indelsign = char
				indelsize = ""
				char = chariter.next()

				# Get the size of the indel
				while(char.isdigit() == True):
					indelsize = indelsize + char
					char = chariter.next()

				# Extract the seq of the indel
				indelstring = char # first character has been previously accessed
				for x in range(int(indelsize) - 1):
					indelstring = indelstring + chariter.next()

				# Generate the final indel string
				indelstring = indelsign + indelsize + indelstring

				indels.append(indelstring)
				continue


			characters.append(char)

		return characters, indels




def InterpretPileupStringMod2(alnstring):
		'''
		Interprets a samtools mpileup string and converts it into
		a list of nts characters and another list of indels. This
		is modified to also count asterisk as part of the result.

		Modified to report position of the gapped alignments. This
		information can then be used to omit the base quality and 
		read positions of these gapped alignments (i.e. > and <).
		'''

		chariter = iter(alnstring)
		indels = []
		characters = []
		gappedAlignPosn = []
		counter = 0
	
		for char in chariter:
			# Ignore the character following the ^ sign that represents the
			# mapping quality of the read
			if char == '^':
				chariter.next()
				continue

			# Deal with character indicating end of read
			if char == '$':
				continue

			# Ignore gapped alignment for RNA
			if char == '>' or char == '<':
				gappedAlignPosn.append(counter)
				counter += 1
				continue

			# Include asterisk for deleted base
			if char == '*':
				characters.append('*')
				counter += 1
				continue


			# Check for the character flag indicating an indel
			if char == '+' or char == '-':
				indelsign = char
				indelsize = ""
				char = chariter.next()

				# Get the size of the indel
				while(char.isdigit() == True):
					indelsize = indelsize + char
					char = chariter.next()

				# Extract the seq of the indel
				indelstring = char # first character has been previously accessed
				for x in range(int(indelsize) - 1):
					indelstring = indelstring + chariter.next()

				# Generate the final indel string
				indelstring = indelsign + indelsize + indelstring

				indels.append(indelstring)
				continue


			characters.append(char)
			counter += 1

		return characters, indels, gappedAlignPosn






class pileup(object):
	'''
	An object type to process the samtools pileup type of data.
	This is capable also of returning some useful information for
	the pileup data that is fed into the algorithm.
	'''
	def __init__(self, refnt, pileupstr, qualstr="NA", mapQstr="NA", readPosnStr="NA"):
		self.refnt 			= refnt.upper()
		self.refnt_lower 	= self.refnt.lower()
		self.pileupstr 		= pileupstr
		self.qualstr 		= qualstr
		self.mapQstr 		= mapQstr
		self.readPosnStr 	= readPosnStr
		self.ntchars, self.indels 	= InterpretPileupString(pileupstr)
		self.pileupstrconv 			= self.convert2character()
		self.ntcounts 				= self.GetNtcounts()
		self.varnt 					= self.GetVarNT()
		self.vaf 					= self.calcVAF()



	def convert2character(self):
		'''
		Convert all the . and , characters into the actual bases.
		Returns a list of all the nt characters.
		'''
		charlist = []
		for nt in self.ntchars:
			if nt == '.':
				charlist.append(self.refnt)
			elif nt == ',':
				charlist.append(self.refnt_lower)
			else:
				charlist.append(nt)
		return charlist


	def GetNtcounts(self):
		'''
		Get counts for each nt
		'''
		counts = {'A': 0, 'T':0, 'G':0, 'C':0, 'N':0}
		for nt in self.pileupstrconv:
			counts[nt.upper()] += 1
		return counts
		

	def GetNtcountsStranded(self):
		'''
		Get counts for each nt
		'''
		counts = {'A': 0, 'T':0, 'G':0, 'C':0, 'N':0, \
		'a':0, 't':0, 'g':0, 'c':0, 'n':0}
		for nt in self.pileupstrconv:
			counts[nt] += 1
		return counts


	def NTcounts(self):
		return self.ntcounts


	def GetVarNT(self):
		'''
		Identify the variant nucleotide of interest in 
		the dataset
		'''
		# sorted_counts is in format of tuple of tuple with first val in first tuple
		# representing the nt and second value representing the count
		sorted_counts = sorted(self.ntcounts.items(), key=operator.itemgetter(1), reverse=True)

		varnt = ""
		if sorted_counts[0][0] != self.refnt:
			varnt = sorted_counts[0][0]
		elif sorted_counts[1][1] != 0:
			varnt = sorted_counts[1][0]
		else:
			varnt = "NA"

		return varnt


	def VarNT(self):
		'''
		Return the variant nucleotide for the position
		'''
		return self.varnt
		

	def calcVAF(self):
		'''
		Calculate the variant NT frequency
		'''
		VAF = float()
		totalbases = float(sum(self.ntcounts.values()))
		if self.varnt == "NA":
			VAF = 0
		elif totalbases == 0:
			VAF = 0
		else:
			VAF = float(self.ntcounts[self.varnt]) / totalbases
		return VAF


	def VAF(self):
		return self.vaf


	def ATGCNcount(self):
		'''
		Returns the ATGCN count of the pileup object
		'''
		return [self.ntcounts['A'], self.ntcounts['T'], self.ntcounts['G'], self.ntcounts['C'], self.ntcounts['N']]


	def depth(self):
		return sum(self.ntcounts.values())


	def varNtCount(self):
		try:
			return self.ntcounts[self.varnt]
		except:
			return 0


	def refNtCount(self):
		try:
			return self.ntcounts[self.refnt]
		except:
			return 0


	def getIndelObject(self):
		return self.indels


	def getIndelCount(self):
		'''
		NOTE that indel in pileup format indicates that indel
		starts in the NEXT position
		'''
		indelCounts = {'ins': 0, 'del': 0}

		if self.indels: # Check if indel list is empty	
			for indel in self.indels:
				if(indel[0] == "+"):
					indelCounts['ins'] += 1
				elif(indel[0] == "-"):
					indelCounts['del'] += 1

		return indelCounts


	def getIndelDetails(self):
		insertionResult = defaultdict(int)
		deletionResult 	= defaultdict(int)
		indelResult 	= {}
		
		if self.indels: # Check if indel list is empty	
			for indel in self.indels:
				if(indel[0] == "+"):
					insertionResult[indel] += 1
				elif(indel[0] == "-"):					
					deletionResult[indel] += 1

		indelResult = {'ins': insertionResult , 'del': deletionResult}


		return indelResult

	
	def getAllConversion(self):
		'''
		Use the nucleotide characters and identify all the nucleotide
		conversion in cases which does not match the reference.
		This is case sensitive (keeps strand info) and returns empty
		list if all bases matches reference.
		'''
		ntConversions = []
		for ntChar in self.ntchars:
			if(ntChar != "." and ntChar != ","):
				ntConversionCurrent = self.refnt + ">" + ntChar
				ntConversions.append(ntConversionCurrent)

		return ntConversions

	def hasMultiNucleotide(self, minNumOfClasses = 3, minCountOfEachClass = 3):
		"""
		Check if the position has nucleotides of the defined number
		of different classes at the specified depth
		"""

		nucClasses = ['A', 'T', 'G', 'C']

		numOfNucTypes = 0
		for nucleotide in nucClasses:
			if self.ntcounts[nucleotide] >= minCountOfEachClass:
				numOfNucTypes += 1

		if numOfNucTypes >= minNumOfClasses:
			return 1
		else:
			return 0



	def meetRDDSpecialFilter(self):
		"""
		Check if the pileup position has nucleotides or indels
		in multiple different classes.
		"""
		self.ntcounts
		indelsInPosn = self.getIndelCount()

		if(indelsInPosn['del'] >= 3 and self.hasMultiNucleotide()):
			return 1
		else:
			return 0



class pileupWithReadPosn(pileup):
	#def __init__(self, refnt, pileupstr, qualstr="NA", mapQstr="NA", readPosnStr="NA"):
	#	super(pileupWithReadPosn, self).__init__(self, refnt, pileupstr)

	def __init__(self, *args, **kwargs):
		super(pileupWithReadPosn, self).__init__(*args, **kwargs)


		if(self.readPosnStr != "NA"):
			#print self.readPosnStr
			self.readPosnVals = self.readPosnStr.strip().split(",")


		#self.ntcountsAtPosn		= self.GetNtcounts()



	def getNucleotidesAtReadPosn(self, readPosition):
		#for readPosnVal in self.readPosnVals.length(:
		#	if()

		ntsAtReadPosnArr = []

		print len(self.readPosnVals)
		print len(self.pileupstrconv)
		print self.indels

		for i in range(len(self.readPosnVals)):

			if(int(self.readPosnVals[i]) == int(readPosition)):
				#print i
				ntsAtReadPosnArr.append(self.pileupstrconv[i])

		return ntsAtReadPosnArr



	def getNucleotidesCountAtReadPosn(self, readPosition):
		ntsAtReadPosnArr = self.getNucleotidesAtReadPosn(readPosition)
		counts = {'A': 0, 'T':0, 'G':0, 'C':0, 'N':0, \
		'a':0, 't':0, 'g':0, 'c':0, 'n':0}
		for nt in ntsAtReadPosnArr:
			counts[nt] += 1
		return counts



	def getMismatchNucleotidesProportionAtReadPosnMinus(self, readPosition):
		ntsAtReadPosnCount = self.getNucleotidesCountAtReadPosn(readPosition)

		try:
			return 1 - float(ntsAtReadPosnCount[self.refnt_lower]) / (ntsAtReadPosnCount['a'] + ntsAtReadPosnCount['t'] + ntsAtReadPosnCount['g'] + ntsAtReadPosnCount['c'])
		except:
			return 0



	def getMismatchNucleotidesProportionAtReadPosnPlus(self, readPosition):
		ntsAtReadPosnCount = self.getNucleotidesCountAtReadPosn(readPosition)

		try:
			return 1 - float(ntsAtReadPosnCount[self.refnt]) / (ntsAtReadPosnCount['A'] + ntsAtReadPosnCount['T'] + ntsAtReadPosnCount['G'] + ntsAtReadPosnCount['C'])
		except:
			return 0



class pileupLine():
	'''
	Class representing a line generated from samtools
	mpileup
	'''
	def __init__(self, line):
		self.line 		= line

		# Check if line is empty?
		# Need deal with cases where there is zero depth but still a line,
		# which causes error


		self.lineArr 	= self.line.strip().split("\t")
		self.chrom 		= self.lineArr[0]
		self.posn 		= self.lineArr[1]
		self.refNt 		= self.lineArr[2]
		self.depth 		= self.lineArr[3]

		# Deal with pileup cases with no depth
		if int(self.depth) != 0:
			self.pileupstr 	= self.lineArr[4]
			self.qualstr 	= self.lineArr[5]
		else:
			self.pileupstr 	= ""
			self.qualstr 	= ""


		



#nts, indels = InterpretPileupString("AT^#GCNGAG+23AGATGCTYUODUKFNHMFMJIDIKMGKKUKJGKLJHDDLJHKD$HJ")
#pileup("A", "AT^#GCNGAG+23AGATGCTYUODUKFNHMFMJIDIKMGKKUKJGKLJHDDLJHKD$HJ", "34234", "asdf")


#apple = pileup("A", "......,,,,,.C.", "34234", "asdf")
#apple = pileup("A", "......,,,,,.C.")

#print apple.VarNT()
#print apple.VAF()


#print apple.ATGCNcount()
#print apple.depth()


