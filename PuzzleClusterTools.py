#tools needed in main software
#include in same file

from math import log,exp,sqrt,pi
from UserDict import UserDict
from random import randint #recall that randint uses closed intervals as input
from time import time

#This class extends the Python dictionary class (UserDict) and stores read information such as the underlying string, the length, and the kmer class probability dictionary:
class Read(UserDict): 
	def __init__(self,readstringinput):
		UserDict.__init__(self)
		self['readstring'] = readstringinput.upper() #Note that we convert to all caps just to be safe
		self['readlength'] = len(self['readstring'])

	def GetKmerClassProbDict(self,k):
		#Since we'll use the kmerclassprobdict of reads repeatedly, we compute it the first time and then store it in memory:
		if(self.has_key(str(k)+'merclassprobdict')):
			return self[str(k)+'merclassprobdict']
		else:			
			kmerprobdict = {}
			totalweight = self['readlength']-k+1
			for kmer in GetKmerList(k):
				kmerprobdict[kmer] = 0
			for startingindex in range(self['readlength']-k+1):
				try:						
					kmerprobdict[self['readstring'][startingindex:startingindex+k]] += float(1)/float(totalweight)
				except KeyError: #This is here to handle letters other than A,T,C, or G (like N) in the genome
					pass
			self[str(k)+'merclassprobdict'] = ConvertKmerProbDictToKmerClassProbDict(kmerprobdict)		
			return self[str(k)+'merclassprobdict']

def GetKmerList(k):
	nucleotides = ['A','T','C','G']
	if(k == 1):
		return nucleotides
	else:		
		previouslist = GetKmerList(k-1)
		newlist = []
		for element in previouslist:
			for nucleotide in nucleotides:
				newlist.append(element + nucleotide)
		return newlist

def GetAverage(listofnumbers):
	totalsum = 0
	for number in listofnumbers:
		totalsum += number
	return float(totalsum)/float(len(listofnumbers))

def GetAverageKmerClassProbDict(k,readorkmerclassprobdictlist):
	if(len(readorkmerclassprobdictlist) == 0):
		raise NameError('Error! Cannot compute the average of a list with no elements!')
	#We allow the input to be either a list of reads or a list of kmer class probability dictionaries. We first try the former:
	try:
		kmerclassprobdictlist = []	
		for read in readorkmerclassprobdictlist:	
			kmerclassprobdictlist.append(read.GetKmerClassProbDict(k))
	except AttributeError:
		kmerclassprobdictlist = readorkmerclassprobdictlist

	averagekmerclassprobdict = {}
	numreads = len(readorkmerclassprobdictlist)	
	
	for kmerclass in GetKmerClassList(k):
		totalkeysum = 0
		for kmerclassprobdict in kmerclassprobdictlist:
			totalkeysum += kmerclassprobdict[kmerclass]
		totalkeyaverage = float(totalkeysum) / float(numreads)
		averagekmerclassprobdict[kmerclass] = totalkeyaverage		

	return averagekmerclassprobdict

def GetDistancesList(k,distancefunctioninput,readsdict,numiterations):
	distancelist = []
	for iterationcounter in range(numiterations):
		print 'Sampling read distances, iteration ',iterationcounter,' of ',numiterations-1		
		randomindex1 = randint(0,len(readsdict)-1) #Recall that the interval following random_integers in inclusive on both ends
		randomindex2 = None
		while(randomindex2 == None):
			randomindex2 = randint(0,len(readsdict)-1)
			if(randomindex2 == randomindex1):
				randomindex2 = None
		newdistance = distancefunctioninput(k,readsdict[randomindex1],readsdict[randomindex2])
		distancelist.append(newdistance)
	return distancelist

def GetInterclusterDistance(k,distancefunctioninput,readsdictinput,clusterset1,clusterset2,numiterations): ###
	distancelist = []
	clusterlist1 = list(clusterset1)
	clusterlist2 = list(clusterset2)
	for iterationcounter in range(numiterations):
		randomindex1 = randint(0,len(clusterlist1)-1)
		randomread1 = readsdictinput[clusterlist1[randomindex1]]
		randomindex2 = randint(0,len(clusterlist2)-1)
		randomread2 = readsdictinput[clusterlist2[randomindex2]]
		distancelist.append(distancefunctioninput(k,randomread1,randomread2))
	return GetAverage(distancelist)

def GetReverseComplement(sequence):
	ComplementDictionary = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	ReverseComplement = ""
	for i in range(-1,-len(sequence)-1,-1):
		if(ComplementDictionary.has_key(sequence[i])):
			ReverseComplement += ComplementDictionary[sequence[i]]
		else:
			ReverseComplement += sequence[i]
	return ReverseComplement

#This method inputs a dictionary whose keys are kmers and whose values are sets, and it outputs a dictionary whose keys are kmer reverse complement equivalence classes (in the form of tuples) and whose values are the unions of the sets over the (one or two) constituent classes
def ConvertSetsDictToClassesDict(setsdict):
	unusedkeysset = set(setsdict.keys())
	classesdict = {}	
	while(len(unusedkeysset) > 0):
		nextkey = unusedkeysset.pop() #We remove the next kmer in the list of unused kmers
		nextkeycomplement = GetReverseComplement(nextkey)
		if(nextkey != nextkeycomplement and setsdict.__contains__(nextkeycomplement)):
			unusedkeysset.discard(nextkeycomplement) #We should be sure to also remove the complement from the list of unused kmers				
			classesdict[(nextkey,nextkeycomplement)] = setsdict[nextkey].union(setsdict[nextkeycomplement])	
			setsdict.__delitem__(nextkey)
			setsdict.__delitem__(nextkeycomplement)	
		elif(nextkey == nextkeycomplement): #Recall that some strings are equal to their own reverse complement, and therefore lie in a single element equivalence class
			classesdict[nextkey] = setsdict[nextkey]
			setsdict.__delitem__(nextkey)
	return classesdict

def GetKmerClassList(k): #This gives a canonical ordering of the kmer classes
	unusedkmerslist = GetKmerList(k)
	kmerclasslist = []
	while(len(unusedkmerslist) > 0):
		nextkey = unusedkmerslist.pop()	
		nextkeycomplement = GetReverseComplement(nextkey)
		if(nextkey != nextkeycomplement):
			unusedkmerslist.remove(nextkeycomplement)
			kmerclasslist.append((nextkey,nextkeycomplement))
		elif(nextkey == nextkeycomplement):
			kmerclasslist.append(nextkey)		
	return kmerclasslist

#This method takes a kmer probability dictionary as input and adds together any two keys in the same reverse complements equivalence class
def ConvertKmerProbDictToKmerClassProbDict(kmerprobdict):
	k = len(kmerprobdict.keys()[0]) #We extract the value of k automatically instead of requiring it as input
	unusedkmerslist = GetKmerList(k)
	kmerclassprobdict = {}	
	while(len(unusedkmerslist) > 0):
		nextkey = unusedkmerslist.pop() #We remove the next kmer in the list of unused kmers
		nextkeycomplement = GetReverseComplement(nextkey)
		if(nextkey != nextkeycomplement):
			unusedkmerslist.remove(nextkeycomplement) #We should be sure to also remove the complement from the list of unused kmers		
			kmerclassprobdict[(nextkey,nextkeycomplement)] = kmerprobdict[nextkey] + kmerprobdict[nextkeycomplement]	
		elif(nextkey == nextkeycomplement): #Recall that some strings are equal to their own reverse complement, and therefore lie in a single element equivalence class
			kmerclassprobdict[nextkey] = kmerprobdict[nextkey]
	return kmerclassprobdict

#This distance is the Jensen-Shannon divergence. Here we are taking into account the occurence of reverse complements.
def KJensenShannonDistance(k,readorkmerclassprobdict1,readorkmerclassprobdict2):
	#We allow either reads or probability dictionaries as input. If the inputs are reads:
	try:
		kmerclassprobdict1 = readorkmerclassprobdict1.GetKmerClassProbDict(k)
		kmerclassprobdict2 = readorkmerclassprobdict2.GetKmerClassProbDict(k)
		averagekmerclassprobdict = GetAverageKmerClassProbDict(k,[kmerclassprobdict1,kmerclassprobdict2])	
	#If the inputs are kmer class probability dictionaries:
	except AttributeError:	
		kmerclassprobdict1 = readorkmerclassprobdict1
		kmerclassprobdict2 = readorkmerclassprobdict2
		averagekmerclassprobdict = GetAverageKmerClassProbDict(k,[kmerclassprobdict1,kmerclassprobdict2])

	distance1 = 0
	distance2 = 0

	#We ignore any events with probability zero from the set of events of the average distribution. 
	for key in averagekmerclassprobdict.keys():	
		if(averagekmerclassprobdict[key] != 0):
			#If either distribution has a zero term, we ignore that term as well (since we are interpretting 0*log(0) as 0):
			if(kmerclassprobdict1[key] != 0):		
				distance1 += kmerclassprobdict1[key]*log(float(kmerclassprobdict1[key])/float(averagekmerclassprobdict[key]))
			if(kmerclassprobdict2[key] != 0):
				distance2 += kmerclassprobdict2[key]*log(float(kmerclassprobdict2[key])/float(averagekmerclassprobdict[key]))
	
	#Finally, we return the average of the two ways of computing the sum:
	return .5*(distance1 + distance2)			
			
def GetMeanAndVariance(datapoints):
	totalsum = 0
	for datapoint in datapoints:
		totalsum += datapoint
	mean = float(totalsum)/float(len(datapoints))
	totalsum2 = 0
	for datapoint in datapoints:
		totalsum2 += pow((datapoint - mean),2)
	variance = float(totalsum2)/float(len(datapoints))
	return [mean,variance]

def gaussianpdf(x,m,sqsd):
	numerator = exp(-float(pow((x-m),2))/float(2*sqsd))	
	denominator = float(sqrt(2*pi*sqsd))
	return (numerator/denominator)

#This method performs the expectation maximization algorithm for fit a set of data points to two mixed Gaussian distributions
def ExpectationMaximizationTwoMixedGaussians(datapoints,maxnumiterations,logsavefilename=None):
	n = len(datapoints)
		
	#We put all the datapoints in a dictionary (with integer keys) called x:
	xdict = {}
	for i in range(1,n+1):
		xdict[i] = datapoints[i-1]		
	
	totalmean,totalvariance = GetMeanAndVariance(xdict.values())

	#We use a simple heuristic for computing initial guesses (this could definitely be improved):
	a1guess = .5
	a2guess = 1 - a1guess
	m1guess = totalmean - sqrt(totalvariance)
	m2guess = totalmean + sqrt(totalvariance)
	sqsd1guess = sqrt(totalvariance)
	sqsd2guess = sqrt(totalvariance)

	#Enter the inital guesses into the dictionary:
	akt = {(1,1):a1guess, (2,1):a2guess}
	mkt = {(1,1):m1guess, (2,1):m2guess}
	sqsdkt = {(1,1): sqsd1guess, (2,1): sqsd2guess}

	t = 1
	print 'Running expectation-maximization algorithm for two mixed normal distributions.'
	print 'Initial guesses:'
	print 'akt[1,1]: ',akt[1,1]
	print 'akt[2,1]: ',akt[2,1]
	print 'mkt[1,1]: ',mkt[1,1]
	print 'mkt[2,1]: ',mkt[2,1]
	print 'sdkt[1,1]: ',sqrt(sqsdkt[1,1])
	print 'sdkt[2,1]: ',sqrt(sqsdkt[2,1])
	print

	if(logsavefilename != None):
		logsavefile = open(logsavefilename,'w')
		logsavefile.write('Initial guesses:\n')
		logsavefile.write('akt[1,1]: '+str(akt[1,1])+'\n')
		logsavefile.write('akt[2,1]: '+str(akt[2,1])+'\n')
		logsavefile.write('mkt[1,1]: '+str(mkt[1,1])+'\n')
		logsavefile.write('mkt[2,1]: '+str(mkt[2,1])+'\n')
		logsavefile.write('sqsdkt[1,1]: '+str(sqsdkt[1,1])+'\n')
		logsavefile.write('sqsdkt[2,1]: '+str(sqsdkt[2,1])+'\n\n')

	while(t <= maxnumiterations):
 		prob_yi_k_cond_xi_theta_t = {}
		
		for i in range(1,n+1):
			denominator = gaussianpdf(xdict[i],mkt[(1,t)],sqsdkt[(1,t)])*akt[(1,t)] + gaussianpdf(xdict[i],mkt[(2,t)],sqsdkt[(2,t)])*akt[(2,t)]
			for k in [1,2]:
				prob_yi_k_cond_xi_theta_t[(k,i,t)] = gaussianpdf(xdict[i],mkt[(k,t)],sqsdkt[(k,t)])*akt[(k,t)]/denominator	
	
		for k in [1,2]:
			akt[k,t+1] = 0
			for i in range(1,n+1):
				akt[k,t+1] += prob_yi_k_cond_xi_theta_t[(k,i,t)]
			akt[k,t+1] = akt[k,t+1]/float(n)
	
		for k in [1,2]:
			mkt[k,t+1] = 0
			for i in range(1,n+1):
				mkt[k,t+1] += float(xdict[i]*prob_yi_k_cond_xi_theta_t[(k,i,t)])
			mkt[k,t+1] = mkt[k,t+1]/(akt[k,t+1]*float(n))		
	
		for k in [1,2]:
			sqsdkt[k,t+1] = 0
			for i in range(1,n+1):
				sqsdkt[k,t+1] += pow((xdict[i] - mkt[k,t+1]),2)*prob_yi_k_cond_xi_theta_t[(k,i,t)]
			sqsdkt[k,t+1] = sqsdkt[k,t+1]/(akt[k,t+1]*float(n))
			
		t += 1
	
	
		print 'After iteration ',t-1,':'	
		print 'akt[1,',t,']: ',akt[1,t]
		print 'akt[2,',t,']: ',akt[2,t]
		print 'mkt[1,',t,']: ',mkt[1,t]
		print 'mkt[2,',t,']: ',mkt[2,t]
		print 'sdkt[1,',t,']: ',sqrt(sqsdkt[1,t]) #It's probably more useful to output the standard deviation
		print 'sdkt[2,',t,']: ',sqrt(sqsdkt[2,t])
		print	

		if(logsavefilename != None):	
			logsavefile.write('After iteration '+str(t-1)+':\n')	
			logsavefile.write('akt[1,'+str(t)+']: '+str(akt[1,t])+'\n')
			logsavefile.write('akt[2,'+str(t)+']: '+str(akt[2,t])+'\n')
			logsavefile.write('mkt[1,'+str(t)+']: '+str(mkt[1,t])+'\n')
			logsavefile.write('mkt[2,'+str(t)+']: '+str(mkt[2,t])+'\n')
			logsavefile.write('sqsdkt[1,'+str(t)+']: '+str(sqsdkt[1,t])+'\n')
			logsavefile.write('sqsdkt[2,'+str(t)+']: '+str(sqsdkt[2,t])+'\n\n')
		
		#We check for convergence:
		if(round(akt[1,t],12) == round(akt[1,t-1],12) and round(akt[2,t],12) == round(akt[2,t-1],12) and round(mkt[1,t],12) == round(mkt[1,t-1],12) and round(mkt[2,t],12) == round(mkt[2,t-1],12) and round(sqsdkt[1,t],12) == round(sqsdkt[1,t-1],12) and round(sqsdkt[2,t],12) == round(sqsdkt[2,t-1],12)):
			print 'Convergence achieved!'
			if(logsavefilename != None):
				logsavefile.write('Convergence achieved!')
			break							

	if(logsavefilename != None):
		logsavefile.close()
	return [akt[1,t],akt[2,t],mkt[1,t],mkt[2,t],sqrt(sqsdkt[1,t]),sqrt(sqsdkt[2,t])]


