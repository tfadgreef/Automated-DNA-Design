# critonsCheck.py 
# This program is part of the TUPack software package and makes use of NuPack 3.0. NuPack can be found on <www.nupack.org>>
# This class creates an object performing criton checks and scoring the results in order to determine whether the sequence is optimum or not.
# written by Zandra Felix Garza
# Eindhoven University of Technology
# 2013
import sys
sys.path.append('/home/s120315/lib64/python2.5/site-packages')

import Strand

import math
import os
import sys
import Tkinter
import tkFileDialog
import copy
import collections 
from collections import defaultdict
        
class CritonsCheck:
    def __init__(self,runNum,nRuns,lambdaFactor,dataFileName,infoFileName,lc=3,strandsObject=[]):
    #def __init__(self,runNum,nRuns,lambdaFactor,lc=3,strandsObject=[],excTemp1=['GACTC','GACTC'],excTemp2=['CTCGA','GACT','GACTC'],excTemp3=['CTCGA','GACTC']):

        self.dataFileName = dataFileName
        self.infoFileName = infoFileName
        

        #Lists holding the Strands
        self.Strands = []

        #List holding repeatedBaseCritons for comparisons
        self.standardCritons=[]

        #Lists holding exceptions for critons
        self.exceptionsTemplate1=[]
        self.exceptionsCountTemplate1=[]        
        self.exceptionsTemplate2=[]
        self.exceptionsCountTemplate2=[]        
        self.exceptionsTemplate3=[]
        self.exceptionsCountTemplate3=[]       
        
        #Lists holding the scores from each criton check
        self.checkCritonsNumber=[]
        self.checkRepeatedCritonsStrand=[]
        self.checkRepeatedBaseCritonsStrand=[]
        self.checkRepeatedAdjacentBaseCritonsStrand=[]
        self.checkRepeatedCritonsNumberStrand=[]
        self.checkRepeatedCritonsSystems=0
        self.GCcritons=0

        #Define strand ID names
        self.SeqId=['T1','T2','T3','F','u','a','Y','inh','cT1','cT2','cT3']

        #Define criton size, number of Runs and the factor for computing the critons cost function
        self.critonSize=lc
        self.runNum = runNum
        self.nRuns = nRuns
        self.lambdaFactor=lambdaFactor
        self.critonWeight=0

        #Define exception sequences
        self.exceptionT1=['GACTC','GACTC']
        #self.exceptionT1=excTemp1
        self.exceptionT2=['CTCGA','GACT','GACTC']
        #self.exceptionT2=excTemp2
        self.exceptionT3=['CTCGA','GACTC']
        #self.exceptionT3=excTemp3

        #Define exceptions for critons
        [self.exceptionsTemplate1,self.exceptionsCountTemplate1] = self._setExceptions(self.exceptionT1)
        [self.exceptionsTemplate2,self.exceptionsCountTemplate2] = self._setExceptions(self.exceptionT2)
        [self.exceptionsTemplate3,self.exceptionsCountTemplate3] = self._setExceptions(self.exceptionT3)

        #Define comparisonCritons
        self.standardCritons=self._setStandardCritons(self.critonSize)


    	 #Defining the local strands object
        try:
            self.Strands=strandsObject
        except(TypeError,IOError):
	     self._readSeqFile(self,'/home/s120315/DNAStrandDesign/TUPack_Upgrade/model/sequences.txt')
      
    # Function for reading the sequence input
    def _readSeqFile(self, seqFile):

        try:
            ID = 1
            for line in seqFile:
                line = line.strip("\n").strip("\r").split("\t")
                newStrand = Strand.Strand(line[1].strip("\r"), line[0], ID)
                self.Strands += [newStrand]
                newStrand.defineN()

                # define the strand type
                if line[2] == "T":
                    newStrand.Type = "template"
                elif line[2] == "B":
                    newStrand.Type = "binding"
                elif line[2] == "I":
                    newStrand.Type = "inhibitor"
                ID += 1  
        except IndexError:
	     sys.exit("The sequence file has a wrong format. Check the user manual for the right file formats.")
            

    #Function that manages the execution of the criton checks. (Local main loop)
    def _checkCritons(self):

        for strand in self.Strands:
            if strand.ID<4:  
                indexNumber=strand.ID-1
                #print "-----------------------------------------------"  
                #print "Starting criton checks for " + str(strand.Name)"
                #print "-----------------------------------------------"   
                currentSequence=strand.Sequence
                critons=self._getCritons(currentSequence)
                #print "critons: "+ str(critons).replace('[','').replace(',','').replace('\'','').replace(']','\n')"
                self.checkCritonsNumber.append(self._getNumberCritons(critons))
                #print "Number of critons for " + str(strand.Name)+": "+ str(self.checkCritonsNumber[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n')"
                #print "Checking the internal uniqueness of the strand""""
                self._checkUniqueCriton(critons,indexNumber)
                #print "Repeated critons in " + str(strand.Name)+": "+ str( self.checkRepeatedCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']',' ')+str(self.checkRepeatedCritonsNumberStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n')"
                #print "Repeated critons of duplicated base for " + str(strand.Name)+": "+ str(self.checkRepeatedBaseCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n')"""
                #print "Adjacent critons of type duplicated base for " + str(strand.Name)+": "+str(self.checkRepeatedAdjacentBaseCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n')"
                #print "Internal Uniqueness check is finished" """         
        #print "-----------------------------------------------" """
        #print "Checking the uniqueness between strands"
        self._checkStrandUniqueness()
        #print "Checks for all strands have been completed.\n"
        #print "Number of repeated critons between strands: " + str( self.checkRepeatedCritonsSystems )+"\n"
        self.critonWeight=self._getCritonWeight(self.lambdaFactor)
        if self.runNum == 1 or (self.runNum%100000) == 0 or self.runNum == self.nRuns:            
            self._writeDataFile(self.dataFileName+str(self.runNum)+".txt")
            self._writeInfoFile(self.infoFileName+str(self.runNum)+".txt")
        return self.critonWeight

    # function for setting the duplicated base critons for comparison base of criton's length
    def _setStandardCritons(self,critonSize):
        critons=[]
        bases=['A','T','C','G']
        temp=[]
        for base in bases:
            temp=[]
            for i in range (0,critonSize):
                temp.extend(base)
            temp = ''.join(temp)
            critons.append(temp)

        return critons

    # function for setting exception critons 
    def _setExceptions(self,exceptionSequences):
        critons=[]
        exceptions=[]
        exceptionsCount=[]
        tally = defaultdict(list)
        
        for sequence in exceptionSequences:
            temp=self._getCritons(sequence)
            critons.extend(temp)

        for i,item in enumerate(critons):
            tally[item].append(i)
            
        for key,locs in tally.items():
            if len(locs)>0:
                exceptions.append(key)
		exceptionsCount.append(len(locs))

        return exceptions,exceptionsCount
    
    
    # function for getting critons for current sequence
    def _getCritons(self,sequence):
        critonNumber = len(sequence)-(self.critonSize-1)
        critons = []

        for i in range (0,critonNumber):
            critons.append(str(sequence[i:i+self.critonSize]))	 
        return critons 
            
    #Inner checks: within each strand         
    def _getNumberCritons(self,critons):

        critonsNumber = len(critons)
        return critonsNumber

    #Get information of duplicates by using a defaultdict to keep a list of all seen locations for any item, and returning those items that were seen more than once.
    def _list_duplicates(self,critonsCurrent):
        tally = defaultdict(list)
        for i,item in enumerate(critonsCurrent):
            tally[item].append(i)
        return ((key,locs) for key,locs in tally.items() if len(locs)>1)

    #Strand check: within the strand
    #Perform check 1: In case of no exceptions critons are not repeated, if there are any exceptions then critons do not appear in the strands more than the predefined number of times
    def _checkUniqueCriton(self,critonsCurrent,strandId):
	 #Initialize local arrays for storage of criton related data
        criton=[]
        location=[]
        repeated=[]
        repeatedAdjacentBase=[]
        repeatedBase=[]
        repeatedCritons=[]
        flag=0
        countGC=0

        #Define the set exceptions according to the template ID
        if strandId == 0:
            exceptionCritons = self.exceptionsTemplate1
            exceptionCount = self.exceptionsCountTemplate1
        elif strandId == 1:
            exceptionCritons = self.exceptionsTemplate2
            exceptionCount = self.exceptionsCountTemplate2           
        elif strandId == 2:
            exceptionCritons = self.exceptionsTemplate3
            exceptionCount = self.exceptionsCountTemplate3

        #Get all the duplicates for the strand
        for duplicate in sorted(self._list_duplicates(critonsCurrent)):
            criton.append(duplicate[0])
            location.append(duplicate[1])

 	 #Check conditions and store those that comply
        for current in criton:
            #If criton is of a repeated base type:
            if str(current) == self.standardCritons[0] or str(current) == self.standardCritons[1] or str(current) == self.standardCritons[2] or str(current) == self.standardCritons[3]:
                #Add the criton and its repetition count to the correspondant variables if it is not already there
                if str(current) not in repeatedCritons:
                    repeatedBase.append(str(current))
                    repeatedBase.append(len(location[criton.index(current)])-1)

                #If the repeated base criton is 'GGG' or 'CCC' consider it for the weight function
                if str(current)== self.standardCritons[2] or str(current)== self.standardCritons[3]:
                    countGC+=len(location[criton.index(current)])

                #If the repeated base criton is adjacent to another criton of the same type, save 
		  #the criton and repetition count information
                for i in range(len(location[criton.index(current)])):
                    if i+1<len(location[criton.index(current)]):
                        if location[criton.index(current)][i+1]==location[criton.index(current)][i]+1:
                            flag+=1
                if flag>0:
                    repeatedAdjacentBase.append(str(current))
                    repeatedAdjacentBase.append(flag)

            #If the criton is not of the repeated base type then consider the existing exceptions 
	     #and save the criton and  its repetition count
            if str(current)!= self.standardCritons[0] and str(current) != self.standardCritons[1] and str(current) != self.standardCritons[2] and str(current) != self.standardCritons[3]:
                #If criton is not found within the exception critons for this strand add it to the repeated critons
		  # if the criton is not already there
                if str(current) not in exceptionCritons:
                    if str(current) not in repeatedCritons:
                        repeatedCritons.append(str(current))
                        repeated.append(len(location[criton.index(current)])-1)
                #If the criton is found within the exception critons for this strand then check if it is repeated
		  #more than the allowed number of times
                if str(current) in exceptionCritons:
                    if len(location[criton.index(current)])>exceptionCount[exceptionCritons.index(current)]:
                        if str(current) not in repeatedCritons:
                            repeatedCritons.append(str(current))
                            repeated.append(len(location[criton.index(current)])-exceptionCount[exceptionCritons.index(current)])
                    
        #Add results to correspondant global variable
        self.checkRepeatedBaseCritonsStrand.append(repeatedBase)
        self.checkRepeatedAdjacentBaseCritonsStrand.append(repeatedAdjacentBase)
        self.checkRepeatedCritonsStrand.append(repeatedCritons)
        self.checkRepeatedCritonsNumberStrand.append(repeated)  
        self.GCcritons += countGC                      
      
    #System checks: between strands
    #Perform check 2: Strand sequence is unique, aside from exceptions
    def _checkStrandUniqueness(self):
	 #Initialize local variables for storing sequences information
        sequences=[]
        critonsTO=[]
        critonsTT=[]
        critonsTTh=[]
	
	 #Get the sequences and their critons
        for strand in self.Strands:
            sequences.append(strand.Sequence)
        critonsTO = self._getCritons(sequences[0])
        critonsTT = self._getCritons(sequences[1])
        critonsTTh = self._getCritons(sequences[2])

        #Obtain the number of repeated critons in the system
        self.checkRepeatedCritonsSystems = self._compareCritonsInStrands(critonsTO,critonsTT,sequences[0])
        self.checkRepeatedCritonsSystems += self._compareCritonsInStrands(critonsTO,critonsTTh,sequences[0])
        self.checkRepeatedCritonsSystems += self._compareCritonsInStrands(critonsTT,critonsTTh,sequences[1])
        
    #Compare the critons from a certain strand with those of the other strands
    def _compareCritonsInStrands(self,strand1,strand2,sequence):
        repeatedCritons = 0
        criton = self.removeDuplicates(strand1)
        originalCritonsS1=self._getCritons(sequence)

        for currentCriton in criton:           
            critonCount=originalCritonsS1.count(currentCriton)+strand2.count(currentCriton)

            if critonCount>4: #Number of maximum accepted repetitions per criton in the system
                repeatedCritons = repeatedCritons+1

        return repeatedCritons

    #Remove duplicated critons for getting a list of existing critons in a certain strand
    def removeDuplicates(self,critons):
        for criton in critons:
            if critons.count(criton)>1:
                critons.remove(criton)
        return critons

    #Calculate the error releated to the critons based on the checks results
    def _getCritonWeight(self,lambdaFactor):
        RB=0 # Initialize the duplicated base critons element
        RC=0 # Initialize the duplicated critons element

	 #Get values for each component of the cost function
        for strand in self.checkRepeatedBaseCritonsStrand:
	    if len(strand)>0:
		    if len(strand) == 2:
			    RB += strand[1]
		    else:
			    RB += sum(strand[1::2])
        RB= (RB +self.GCcritons)*(2)
        for strand in self.checkRepeatedCritonsNumberStrand:
	    if len(strand)>0:
		    RC += sum(strand)
        RS=self.checkRepeatedCritonsSystems #Duplicated critons in the system

        #Compute cost function
        weight=float((RB+RC+RS)*(RB+RC+RS))*lambdaFactor        
        
        return weight

    #Write data file containing the number of critons, duplicated critons per strand and their counts, as well as an indication of which ones are of duplicated base type
    #and wether they are adjacent to each other. Also, the total number of repeated critons is included.
    def _writeDataFile(self,fileName):
        f=open(fileName, "w")
        f.write('#Number of critons\n')
        f.write(str(self.checkCritonsNumber).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
        f.write('#Repeated critons per strand\tCritons count for repeated critons\n')
        for item in self.checkRepeatedCritonsStrand:
            f.writelines(str(item).replace('[','').replace(',','').replace('\'','').replace(']',' ')+
                         str(self.checkRepeatedCritonsNumberStrand[self.checkRepeatedCritonsStrand.index(item)]).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
        f.write('#Repeated base critons per strand\tCritons count for repeated critons\n')
	for item in self.checkRepeatedBaseCritonsStrand:
            f.writelines(str(item).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
        f.write('#Repeated adjacent base critons per strand\tCritons count for adjacent repeated critons cases\n')
        for item in self.checkRepeatedAdjacentBaseCritonsStrand:
            f.writelines(str(item).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
        f.write('#Number of repeated critons between strands\n')
        f.write(str(self.checkRepeatedCritonsSystems)+'\n')
        f.write('#Critons contribution to global error function\n'+str(self.critonWeight)+'\n')
        f.close()

    #Write info file containing the number of critons, duplicated critons per strand and their counts, as well as an indication of which ones are of duplicated base type
    #and wether they are adjacent to each other. Also, the total number of repeated critons is included. The information is organized by 
    #strand
    def _writeInfoFile(self,fileName):
        f=open(fileName,"w")
        f.write('Critons check info file\n')
        f.write('------------------------------------------------------------------------------\n')
        for strand in self.Strands:
            if strand.ID<4:
                indexNumber=strand.ID-1
                f.write('Starting criton checks for\t'+ str(strand.Name)+'\n'+'Number of critons for '+str(strand.Name)+": "+str(self.checkCritonsNumber[indexNumber])+'\n'+'Critons for '+str(strand.Name)+": "+str(self._getCritons(strand.Sequence)).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
                f.write('Checking the internal uniqueness of the strand\nRepeated critons in '+str(strand.Name)+": "+ str( self.checkRepeatedCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']',' ')
                        +str(self.checkRepeatedCritonsNumberStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
                f.write('Repeated critons of duplicated base for ' +str(strand.Name)+": "+ str(self.checkRepeatedBaseCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
                f.write('Adjacent critons of duplicated base for '+str(strand.Name)+": "+str(self.checkRepeatedAdjacentBaseCritonsStrand[indexNumber]).replace('[','').replace(',','').replace('\'','').replace(']','\n'))
                f.write('Internal Uniqueness check is finished\n')
                f.write('------------------------------------------------------------------------------\n')
        f.write('Checking the uniqueness between strands\n'+'Number of repeated critons between strands: ' + str( self.checkRepeatedCritonsSystems )+'\n')
        f.write('Checks for all strands have been completed.\n')
        f.write('Critons contribution to global error function: '+ str(self.critonWeight)+'\n')
        f.close()
