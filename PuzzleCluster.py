#main program of clustering software

from PuzzleClusterTools import *
import sys
from time import time

def ConvertReadsfileToDict(readsfilename):
	readsfile = open(readsfilename,'r')
	readsdict = {}
	proceedprocessingfile = True
	readsfile.readline()
	readindex = 0
	while(proceedprocessingfile):
		newreadstring = ''
		proceedprocessingread = True
		while(proceedprocessingread):
			newline = ''.join(readsfile.readline().split())
			if(newline == ''):
				proceedprocessingread = False
				proceedprocessingfile = False
			elif(newline[0] == '>'):
				proceedprocessingread = False
			elif(newline[0] == 'A' or newline[0] == 'G' or newline[0] == 'T' or newline[0] == 'C' or newline[0] == 'N'):
				newreadstring += newline
			else:
				print 'Error encountered while processing',readsfilename,'!'
		readsdict[readindex] = Read(''.join(newreadstring.split()).upper())
		readindex += 1
	readsfile.close()
	return readsdict

#This method inputs an index-to-read dictionary readsdict and a segment length, and it returns a list of sets of read (indices) which agree on a word of length at least segmentlength. Here we do not distinguish between a word and its reverse complement.
def GetWordClassAgreementList(readsdict,segmentlength):
	segmentdict = {}	
	for readindex in readsdict.keys():
		print 'Looking for word agreements of length '+str(segmentlength)+', or higher - read index: ',readindex,' of ',len(readsdict)		
		readlength = readsdict[readindex]['readlength'] #Note that the reads can vary in length
		for startindex in range(readlength-segmentlength+1):
			newsegment = readsdict[readindex]['readstring'][startindex:startindex + segmentlength]			
			if(segmentdict.has_key(newsegment)):
				segmentdict[newsegment].add(readindex)
			else: #For the first time a word appears in the set, we add a new entry to segmentdict
				segmentdict[newsegment] = set([readindex])

	#Next, we collapse any two keys of segmentdict which are reverse complements:
	print 'Collapsing word agreement dictionary to equivalence classes...'
	wordclassagreementdict = ConvertSetsDictToClassesDict(segmentdict)
	print 'Finished.'	

	#We now restrict to those words occuring at least twice:
	repeatedwordclassagreementdict = {} #repeatedwordclassagreementdict will be the restriction of wordclassagreementdict to those keys whose value is a set of length at least 2
	for wordclass in wordclassagreementdict.keys():
		if(len(wordclassagreementdict[wordclass]) >= 2):
			repeatedwordclassagreementdict[wordclass] = wordclassagreementdict[wordclass]	
	
	#Finally, we return the resultant sets of agreeing (indices of) reads:
	return repeatedwordclassagreementdict.values()

#This method uses the agreementlist of the previous method to group together all reads into word-agreeing equivalence classes:
def ClusterReadsByAgreementList(readsdict,agreementlist):
	readindextoclusterindexdict = {}	
	clusterdict = {}
	currentclusternumber = 0
	readsetindex = 0

	for readset in agreementlist:
		print 'Using agreement set ',readsetindex,' of ',len(agreementlist)
		newsupercluster = set([])
		for readindex in readset:
			try:		
				newsupercluster.update(clusterdict[readindextoclusterindexdict[readindex]])
				oldclusterindex = readindextoclusterindexdict[readindex]		
				for partnerreadindex in clusterdict[readindextoclusterindexdict[readindex]]:
					readindextoclusterindexdict[partnerreadindex] = currentclusternumber
	
				clusterdict.pop(oldclusterindex)
						
			except KeyError:
				newsupercluster.add(readindex)
				readindextoclusterindexdict[readindex] = currentclusternumber	
			
		readsetindex += 1

		clusterdict[currentclusternumber] = newsupercluster
		currentclusternumber += 1				
			
	#Finally, we add the remaining read indices to singleton clusters
	for readindex in range(len(readsdict)):
		print readindex,' of ',len(readsdict)
		if(readindextoclusterindexdict.has_key(readindex)):
			pass
		else:
			clusterdict[currentclusternumber] = set([readindex])
			currentclusternumber += 1

	reorderedclusterdict = {}	
	newclusterindex = 0
	for oldclusterindex in clusterdict.keys():
		reorderedclusterdict[newclusterindex] = clusterdict[oldclusterindex]
		newclusterindex += 1		

	return reorderedclusterdict

#This method partitions the feature vector space into a grid, and then picks a read from the highest density grid
def GetHighDensityProbDict(kuniv,readsdict,numbins):
	maxdensity = 0
	maxdensitykeyslist = []
	clusterdict = {}
	
	#This is analogous to a choice of projection: 
	roundingchoicedict = {}
	for kmerclass in GetKmerClassList(kuniv):	
		roundingchoicedict[kmerclass] = randint(0,1)

	for readindex in readsdict.keys():
		kmerclassprobdict = readsdict[readindex].GetKmerClassProbDict(kuniv)
		roundedkmerclassprobdict = {}		
		for kmerclass in kmerclassprobdict.keys():		
			if(roundingchoicedict[kmerclass] == 0):	
				roundedkmerclassprobdict[kmerclass] = round(numbins*kmerclassprobdict[kmerclass])/float(numbins) #We round to the nearest gridpoint
			elif(roundingchoicedict[kmerclass] == 1):
				roundedkmerclassprobdict[kmerclass] = (round(numbins*kmerclassprobdict[kmerclass] + .5)-.5)/float(numbins) #We round to the nearest half-gridpoint
		#Now we convert roundedkmerclassprobdict to a tuple so that it's hashable:
		roundedkmerclassproblist = []
		for kmerclass in GetKmerClassList(kuniv):
			roundedkmerclassproblist.append(roundedkmerclassprobdict[kmerclass])
		roundedkmerclassprobtuple = tuple(roundedkmerclassproblist)
		if(clusterdict.has_key(roundedkmerclassprobtuple)):
			clusterdict[roundedkmerclassprobtuple].add(readindex)
		else:
			clusterdict[roundedkmerclassprobtuple] = set([readindex])

	#We would like to keep track of all clusters sharing the maximum density:		
	for clusterkey in clusterdict.keys():
		if(len(clusterdict[clusterkey]) > maxdensity):
			maxdensity = len(clusterdict[clusterkey])
			maxdensitykeyslist = [clusterkey]
		elif(len(clusterdict[clusterkey]) == maxdensity):
			maxdensitykeyslist.append(clusterkey)
			
	print 'Found a read in a region with density: ',maxdensity
	
	randomindex = randint(0,len(maxdensitykeyslist)-1)	
	randommaxdensitykey = maxdensitykeyslist[randomindex]

	return GetAverageKmerClassProbDict(kuniv,[readsdict[readindex] for readindex in clusterdict[randommaxdensitykey]])

#This is the main method of the program:
def PuzzleCluster(kinput,distancefunctioninput,mso,readsdictinput,groupsdictinput,minsphereradius,typicalsphereradius,radiusincrement,maxsphereradius,maxnumspecies,stoppingratio,singletondistancethreshold,numbins,numinitialcenterguesses,numpointstosample,logsavefilename,commentary):
	logsavefile = open(logsavefilename,'w')
	logsavefile.write('Using parameters: \n')
	logsavefile.write('Distance function: '+distancefunctioninput.__name__+'\n')	
	logsavefile.write('k: '+str(kinput)+'\n')
	logsavefile.write('mso: '+str(mso)+'\n')
	logsavefile.write('minsphereradius: '+str(minsphereradius)+'\n')	
	logsavefile.write('radiusincrement: '+str(radiusincrement)+'\n')
	logsavefile.write('typicalsphereradius: '+str(typicalsphereradius)+'\n')
	logsavefile.write('maxsphereradius: '+str(maxsphereradius)+'\n')		
	logsavefile.write('stoppingratio: '+str(stoppingratio)+'\n')
	logsavefile.write('singletondistancethreshold: '+str(singletondistancethreshold)+'\n')
	logsavefile.write('numbins: '+str(numbins)+'\n')
	if(commentary):
		print'Using parameters: \n'
		print'Distance function: '+distancefunctioninput.__name__+'\n'
		print'k: '+str(kinput)+'\n'
		print'mso: '+str(mso)+'\n'
		print'minsphereradius: '+str(minsphereradius)+'\n'
		print'radiusincrement: '+str(radiusincrement)+'\n'
		print 'typicalsphereradius: '+str(typicalsphereradius)+'\n'
		print'maxsphereradius: '+str(maxsphereradius)+'\n'		
		print'stoppingratio: '+str(stoppingratio)+'\n'
		print'singletondistancethreshold: '+str(singletondistancethreshold)+'\n'
		print'numbins: '+str(numbins)+'\n'

	kuniv = kinput
	distancefunctionuniv = distancefunctioninput
	groupsdict = groupsdictinput.copy()
	readsdict = readsdictinput.copy()	
	stoppingquantity = int(float(len(readsdict))*float(stoppingratio))

	maxpointsinsphere = 0
	clustersearchindex = 0
	clusterslist = []

	while(len(readsdict) > stoppingquantity and len(clusterslist) < maxnumspecies):
		logsavefile.write('Looking for an approximate center of a cluster by finding a maximal location for a center of radius '+str(minsphereradius)+'. We proceed by picking random high density locations (found by a grid method) '+str(numinitialcenterguesses)+' times, each time performing the centroid method for further maximization, and then we take the maximum of these trials.\n')
		if(commentary):
			print'Looking for an approximate center of a cluster by finding a maximal location for a center of radius '+str(minsphereradius)+'. We proceed by picking random high density locations (found by a grid method) '+str(numinitialcenterguesses)+' times, each time performing the centroid method for further maximization, and then we take the maximum of these trials.\n'

		optimalcluster = set([])
		optimalcenterprobdict = None
		for initialguessindex in range(numinitialcenterguesses):	
			centerprobdict = GetHighDensityProbDict(kuniv,readsdict,numbins)
			#We guess a new location for the sphere and optimize that location using the centroid method:			
			potentialclusterindexset = set([]) #In the beginning we should compare new spheres to the empty sphere					
			continuecentroidmethod = True			
			while(continuecentroidmethod):
				proxypotentialclusterindexset = set([])
				for readindex in readsdict.keys():
					readprobdict = readsdictinput[readindex].GetKmerClassProbDict(kuniv)
					if(distancefunctionuniv(kuniv,readprobdict,centerprobdict) <= typicalsphereradius):
						proxypotentialclusterindexset.add(readindex)
				if(len(proxypotentialclusterindexset) <= len(potentialclusterindexset)):
					continuecentroidmethod = False
				elif(len(proxypotentialclusterindexset) > len(potentialclusterindexset)):				
					potentialclusterindexset = proxypotentialclusterindexset
					centerprobdict = GetAverageKmerClassProbDict(kuniv,[readsdictinput[readindex] for readindex in potentialclusterindexset])
			#Now we compare the result to our current best location:
			if(len(potentialclusterindexset) > len(optimalcluster)):	
				optimalcluster = potentialclusterindexset
				optimalcenterprobdict = centerprobdict
		
		logsavefile.write('Finished finding approximate center of cluster. Found a location for a sphere containing '+str(len(optimalcluster))+' reads.\n')
		if(commentary):
			print 'Finished finding approximate center of cluster. Found a location for a sphere containing '+str(len(optimalcluster))+' reads.\n'

		currentsphereradius = minsphereradius	
	
		try:
			currentcenterprobdict = optimalcenterprobdict
		except UnboundLocalError:
			raise NameError('Error! Failed to find an initial sphere containing any reads. Perhaps try a larger typical sphere radius.')

		levelslist = [optimalcluster]
		levelsradiilist = [currentsphereradius]		
		normaveragedisplacementslist = []
		increasesthreshold = 2
		numincreasesnormaveragedisplacements = 0

		logsavefile.write('Looking for appropriate radius of cluster...\n')
		if(commentary):
			print'Looking for appropriate radius of cluster...\n'
	
		continueincreasingradius = True
		while(continueincreasingradius):	
			currentsphereradius += radiusincrement			
			levelsradiilist.append(currentsphereradius)

			logsavefile.write('Trying radius '+str(currentsphereradius)+'\n')
			if(commentary):			
				print'Trying radius '+str(currentsphereradius)+'\n'

			#We fill up the sphere of the next biggest radius:
			levelslist.append(set([]))			
			for readindex in readsdict.keys():
				readprobdict = readsdictinput[readindex].GetKmerClassProbDict(kuniv)
				if(distancefunctionuniv(kuniv,readprobdict,currentcenterprobdict) <= currentsphereradius):
					levelslist[-1].add(readindex)	

			logsavefile.write('Contains '+str(len(levelslist[-1]))+' reads.\n')
			if(commentary):
				print'Contains '+str(len(levelslist[-1]))+' reads.\n'
			
			if(currentsphereradius >= maxsphereradius):
				clustercoreset = levelslist[-1]
				continueincreasingradius = False

				logsavefile.write('Reached maximum allowable radius.\n')
				if(commentary):
					print'Reached maximum allowable radius.\n'

				if(currentsphereradius >= maxsphereradius):
					clustercoreset = levelslist[-1]			
					continueincreasingradius = False
	
					print 'Reached maximum allowable radius.'
					logsavefile.write('Reached maximum allowable radius.')

			elif(len(levelslist[-1]) > 0):
				levelcentroidprobdict = GetAverageKmerClassProbDict(kuniv,[readsdictinput[readindex] for readindex in levelslist[-1]])
				normaveragedisplacementslist.append(distancefunctionuniv(kuniv,currentcenterprobdict,levelcentroidprobdict))

				logsavefile.write('Latest norm average displacement: '+str(normaveragedisplacementslist[-1])+'\n')
				if(commentary):
					print'Latest norm average displacement: '+str(normaveragedisplacementslist[-1])+'\n'
				
				if(len(normaveragedisplacementslist) >= increasesthreshold):
					if(normaveragedisplacementslist[-1] >= normaveragedisplacementslist[-2]):
						numincreasesnormaveragedisplacements += 1			
					elif(normaveragedisplacementslist[-1] < normaveragedisplacementslist[-2]):
						numincreasesnormaveragedisplacements = 0
				if(numincreasesnormaveragedisplacements == increasesthreshold):
					clustercoreset = levelslist[-increasesthreshold - 1]

					logsavefile.write('Finished finding radius of cluster core.\n')
					if(commentary):		
						print'Finished finding radius of cluster core.\n'
	
					continueincreasingradius = False
				elif(numincreasesnormaveragedisplacements != increasesthreshold):
					currentcenterprobdict = GetAverageKmerClassProbDict(kuniv,[readsdictinput[readindex] for readindex in levelslist[-1]])	
					
		logsavefile.write('Done looking for radius. Cluster core contains '+str(len(clustercoreset))+' reads.\n')
		if(commentary):
			print'Done looking for radius. Cluster core contains '+str(len(clustercoreset))+' reads.\n'

		clusterslist.append(clustercoreset)	
		
		#Now we add to the sphere all the reads which are connected by word agreement to reads in the sphere:
		for groupindex in groupsdict.keys():
			for readindex in groupsdict[groupindex]:
				if(readindex in clusterslist[-1]): #If a single read from a group is in the sphere, we add the entire group to the sphere		
					clusterslist[-1].update(groupsdict[groupindex])	
					groupsdict.pop(groupindex) #We remove that group from groupsdict to make future searches quicker
					break											
		for readindex in clusterslist[-1]: 
			readsdict.pop(readindex) #We also remove every read in this cluster from readsdict to make future searches quicker	

		logsavefile.write('Added a cluster with '+str(len(clusterslist[-1]))+' elements!\n')	
		if(commentary):		
			print'Added a cluster with ',len(clusterslist[-1]), ' elements!'
				
		logsavefile.write('Number of clusters: '+str(len(clusterslist))+'\n')
		logsavefile.write('Number of unclustered reads: '+str(len(readsdict))+'\n')
		if(commentary):
			print'Number of clusters: '+str(len(clusterslist))+'\n'
			print'Number of unclustered reads: '+str(len(readsdict))+'\n'	
	
	logsavefile.write('Finished discovering primary clusters.\n')
	if(commentary):
		print 'Finished discovering primary clusters.'		
	
	#We now add any remaining read to its closest cluster, if the distance to that cluster is under singletondistancethreshold:
	centroidslist = []
	for cluster in clusterslist:
		clusterreadlist = [readsdictinput[readindex] for readindex in cluster]	
		centroidslist.append(GetAverageKmerClassProbDict(kuniv,clusterreadlist)) 

	for readindex in readsdict.keys():
		closestclusterindex = None
		distancetoclosestcluster = None
		readprobdict = readsdictinput[readindex].GetKmerClassProbDict(kuniv)
		for potentialclosestclusterindex in range(len(clusterslist)):
			potentialdistancetoclosestcluster = distancefunctionuniv(kuniv,readprobdict,centroidslist[potentialclosestclusterindex])		
			if(closestclusterindex == None):#The first cluster is automatically the current closest
				distancetoclosestcluster = potentialdistancetoclosestcluster
				closestclusterindex = potentialclosestclusterindex
			else: 
				if(potentialdistancetoclosestcluster < distancetoclosestcluster):
					distancetoclosestcluster = potentialdistancetoclosestcluster
					closestclusterindex = potentialclosestclusterindex
		if(distancetoclosestcluster <= singletondistancethreshold):
			clusterslist[closestclusterindex].add(readindex)

			logsavefile.write('Added read '+str(readindex)+' to cluster '+str(closestclusterindex)+' since its distance was below the specified threshold.\n')	
			if(commentary):
				print 'Added read ',readindex,' to cluster ',closestclusterindex, ' since its distance was below the specified threshold.'			
			readsdict.pop(readindex)
		
	#Then compile a list of the centroids of each cluster:		
	logsavefile.write('Computing centroids and average size of each cluster...\n')
	if(commentary):
		print 'Computing centroids and average size of each cluster...\n'

	centroidslist = []
	averagedistancetocenterlist = []
	for cluster in clusterslist:
		centroidslist.append(GetAverageKmerClassProbDict(kuniv,[readsdictinput[readindex] for readindex in cluster])) 
		distancestocenterlist = [distancefunctionuniv(kuniv,readsdictinput[readindex].GetKmerClassProbDict(kuniv),centroidslist[-1]) for readindex in cluster]
		averagedistancetocenterlist.append(GetAverage(distancestocenterlist))	

	#Next we perform hierarchical clustering on the clusters we have so far, using clusterdistancethreshold as the stopping threshold:
	#Find the minimum intercluster distance. Since we're assuming there are only a few large clusters, we can ignore the case of ties.		
	continuehierarchicalclustering = True		
	while(continuehierarchicalclustering and len(clusterslist) >= 2):	
		#We next find the minimum intercluster distance and the pair that achieves this minimum:
		mininterclusterdistance = None
		mininterclusterdistancepair = None
		for clusterindex1 in range(len(clusterslist)):
			for clusterindex2 in range(clusterindex1+1,len(clusterslist)):
				potentialmininterclusterdistance = GetInterclusterDistance(kuniv,distancefunctionuniv,readsdictinput,clusterslist[clusterindex1],clusterslist[clusterindex2],numpointstosample)
				potentialmininterclusterdistancepair = (clusterindex1,clusterindex2)
				if(mininterclusterdistance == None): 
					mininterclusterdistance = potentialmininterclusterdistance
					mininterclusterdistancepair = potentialmininterclusterdistancepair	
				elif(potentialmininterclusterdistance < mininterclusterdistance):
					mininterclusterdistance = potentialmininterclusterdistance
		if(mininterclusterdistance <= max(averagedistancetocenterlist[mininterclusterdistancepair[0]],averagedistancetocenterlist[mininterclusterdistancepair[1]])):	

			logsavefile.write('Added cluster '+str(mininterclusterdistancepair[1])+' ('+str(len(clusterslist[mininterclusterdistancepair[1]]))+' elements) to cluster '+str(mininterclusterdistancepair[0])+' ('+str(len(clusterslist[mininterclusterdistancepair[0]]))+' elements).\n')
			logsavefile.write('Cluster '+str(mininterclusterdistancepair[1])+':\n')
			logsavefile.write(str(clusterslist[mininterclusterdistancepair[1]])+'\n')	
			logsavefile.write('Cluster '+str(mininterclusterdistancepair[0])+':\n')
			logsavefile.write(str(clusterslist[mininterclusterdistancepair[0]])+'\n')
			if(commentary):
				print 'Added cluster ',mininterclusterdistancepair[1],' ('+str(len(clusterslist[mininterclusterdistancepair[1]])),' elements) to cluster ',mininterclusterdistancepair[0],' ('+str(len(clusterslist[mininterclusterdistancepair[0]])),' elements).'
				print 'Cluster '+str(mininterclusterdistancepair[1])+':'
				print str(clusterslist[mininterclusterdistancepair[1]])
				print 'Cluster '+str(mininterclusterdistancepair[0])+':'
				print str(clusterslist[mininterclusterdistancepair[0]])
	
			clusterslist[mininterclusterdistancepair[0]].update(clusterslist[mininterclusterdistancepair[1]])
			centroidslist[mininterclusterdistancepair[0]] = GetAverageKmerClassProbDict(kuniv,[readsdictinput[readindex] for readindex in clusterslist[mininterclusterdistancepair[0]]])			
			distancestocenterlist = [distancefunctionuniv(kuniv,readsdictinput[readindex].GetKmerClassProbDict(kuniv),centroidslist[mininterclusterdistancepair[0]]) for readindex in clusterslist[mininterclusterdistancepair[0]]]
			averagedistancetocenterlist[mininterclusterdistancepair[0]] = GetAverage(distancestocenterlist)
			clusterslist.pop(mininterclusterdistancepair[1])
			centroidslist.pop(mininterclusterdistancepair[1])
			averagedistancetocenterlist.pop(mininterclusterdistancepair[1])		
		else:
			continuehierarchicalclustering = False

	singletonsset = set(readsdict.keys())
	clusterslistproxy = clusterslist[:]	

	#Now we classify the clusters smaller than the desired threshold as unclustered:
	clusterindexcounter = 0
	for cluster in clusterslist:
		if(len(cluster) <  stoppingquantity - float(len(readsdict))/float(len(clusterslist))):
			logsavefile.write('Removed cluster '+str(clusterindexcounter)+' since it only contained '+str(len(cluster))+' reads.\n')
			if(commentary):
				print'Removed cluster '+str(clusterindexcounter)+' since it only contained '+str(len(cluster))+' reads.\n'

			singletonsset.update(cluster)
			clusterslistproxy.remove(cluster)
		clusterindexcounter += 1
	clusterslist = clusterslistproxy

	for clusterindex in range(len(clusterslist)):
		logsavefile.write('Cluster '+str(clusterindex)+': '+str(len(clusterslist[clusterindex]))+' elements.\n')
		if(commentary):
			print 'Cluster '+str(clusterindex)+': '+str(len(clusterslist[clusterindex]))+' elements.\n'

	logsavefile.write('Number of unclustered reads: '+str(len(singletonsset))+'\n')
	if(commentary):		
		print 'Number of unclustered reads: '+str(len(singletonsset))+'\n'

	logsavefile.write(str(clusterslist))

	logsavefile.close()

	return clusterslist

def main(args):
	kuniv = 4 #This is the universal value for the length of kmers, which will be used throughout the program.
	mso = 25 #This stands for the minimum significant overlap, a number to be determined theoretically.	
	distancefunction = KJensenShannonDistance

	starttime = time()
	readsfilename = args[1]	
	readsdict = ConvertReadsfileToDict(readsfilename)

	wordagreementlist = GetWordClassAgreementList(readsdict,mso)	
	wordagreementclusters = ClusterReadsByAgreementList(readsdict,wordagreementlist)

	maxnumemiterations = 500
	numemdatapoints = 10000	

	datapoints = GetDistancesList(kuniv,distancefunction,readsdict,numemdatapoints)
	curvefittingdata = ExpectationMaximizationTwoMixedGaussians(datapoints,maxnumemiterations)
	mean1 = curvefittingdata[2]
	mean2 = curvefittingdata[3]
	sd1 = curvefittingdata[4]
	sd2 = curvefittingdata[5]

	stoppingratio = .05
	singletondistancethreshold = mean1 - sd1
	numbins = 30
	logsavefilename = readsfilename+'_log.txt'	
	commentary = True
	minsphereradius = float(mean1 - 0.50*sd1)/float(2)
	radiusincrement = 0.10*sd1
	typicalsphereradius = float(mean1)/float(2)
	maxsphereradius = float(mean2 + sd2)/float(2)
	maxnumspecies = 10
	numinitialcenterguesses = 5
	numpointstosample = 1000

	puzzle = PuzzleCluster(kuniv,distancefunction,mso,readsdict,wordagreementclusters,minsphereradius,typicalsphereradius,radiusincrement,maxsphereradius,maxnumspecies,stoppingratio,singletondistancethreshold,numbins,numinitialcenterguesses,numpointstosample,logsavefilename,commentary)
	print puzzle

	SaveClusters(puzzle,readsfilename)

	print 'Total elapsed time ',time()-starttime,' seconds.'

#This method saves the clusters in a format convenient to evaluating the clustering results in Matlab program:
def SaveClusters(clusterslist,savefilename):
	for clusterindex in range(len(clusterslist)):
		clustersavefilename = savefilename+'_cluster_'+str(clusterindex)
		clustersavefile = open(clustersavefilename,'w')
		clustersavefile.write(str(list(clusterslist[clusterindex]))[1:-1])
		print 'Wrote cluster '+str(clusterindex)+' data to file '+clustersavefilename+'.'
		clustersavefile.close()	

if __name__ == "__main__":
	sys.exit(main(sys.argv))



