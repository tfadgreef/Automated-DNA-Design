# seqGen class
# This program is part of the TUPack software package and makes use of NuPack 3.0. NuPack can be found on <www.nupack.org>>
# This class creates an object performing simulated annealing and random sequence mutations to find the global optimum to the target free energies.
# For more info about the program, please read the user manual.
#
# written by Sander Rodenburg
# modified by Zandra Felix Garza
# Eindhoven University of Technology
# 2012
# modified in 2013

import sys
sys.path.append('/home/hroekel/lib64/python2.5/site-packages')

# importing the NUPACK object, for running complexes and for reading the NuPack output files
import nupack

# importing the Strand library, containing the objects Strand and StrandRelation for implementing the DNA network
import Strand

# importing the CritonCheck object, for running the critons checks.
import critonsCheck

# importing the ppserver library, for running processes in parallel
import pp

import randomizer

# importing other libraries
import random
import math
import time
import os
import sys
import Tkinter
import tkFileDialog
import copy

class calcConc:
  
  def __init__(self, nRuns, T_start, T_end, lambdaFactor, maxCxSize=1, seqFile=None, bindingFile=None, saveRuns=[], mode="SA"):

    # defining seed for random
    random.seed(time.time())
    #random.seed(1)

    # lists holding the Strands, duplexes and Inhibitors
    self.Strands = []
    self.Duplexes = []
    self.Inhibitors = []

    if mode not in ["SA", "random"]:
      raise sys.exit("'mode' parameter should be either 'SA', or 'random'")

    # list for keeping the order in which the sequences should be adjusted
    self.adjustOrder = []

    # NuPack maximum complex size
    self.maxCxSize = maxCxSize

    # runs in which the program should make a distinct save
    self.saveRuns = saveRuns

    #Section included for the initialization of parallel processing related parameter
    ######################################################
    # tuple of all parallel python servers to connect with
    self.ppservers = ()

    if len(sys.argv) > 1:
      self.ncpus = int(sys.argv[1])
      # Creates jobserver with ncpus workers
      self.job_server = pp.Server(ncpus, ppservers=self.ppservers)
    else:
      # Creates jobserver with automatically detected number of workers
      self.job_server = pp.Server(ppservers=self.ppservers)
    #######################################################
    
    # critons check object
    self.critonsChck=()

    # criton weight for error function
    self.critonWeight=0

    # total delta error and number of total positive delta errors
    self.totalDError = 0
    self.totalPosError = 0

    self.gibbsError = 0

    # if the files are not defined in the class parameters, a filechooser pops up
    try:
      self.seqFile = open(seqFile, "r")
    except (TypeError, IOError):
      self.seqFile = self._chooseFile("Select sequence file")
      
    try:
      self.bindingFile = open(bindingFile, "r")
    except (TypeError, IOError):
      self.bindingFile = self._chooseFile("Select binding file")

    # if there are still no files selected, raise an error
    if self.seqFile == None or self.bindingFile == None:
      sys.exit("One or more input files are missing.")

    #defining file path
    wd = os.path.dirname(self.seqFile.name)

    # changing directory to input file path
    ioDir = wd+"/NuPackIO"
    if not os.path.exists(ioDir):
      os.mkdir(ioDir)
    os.chdir(wd+"/NuPackIO")

    # initializing the NUPACK object in silent mode
    self.nupack = nupack.NUPACK(prefix=os.getcwd()+"/cxFile",
                                paramFile=os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))+"/parameters.txt",
                                silent=True, degenerate=False, material="dna", dangles="all", ordered=True, timeonly=False, mfe=True, quiet=True, sort=1)

    # reading input files for initial sequence generations   
    self._readSeqFile(self.seqFile)
    self._readBindingFile(self.bindingFile)

    # close the files
    self.seqFile.close()
    self.bindingFile.close()

    # adjust initial sequences
    self._adjustSequenceSet()
    self._adjustInhibitors()
    """
    if mode == "SA":
      print "Mode:\tSimulated Annealing (SA)"
      
      # run the mainloop for simulated annealing and sequence mutation
      self._mainloop(nRuns, T_start, T_end, lambdaFactor, screenOut=False)
      print "Done."    
      print "Files stored in " +os.getcwd()

      if mode == "SA":
        if nRuns > 100:
          avgDError = self.totalDError/self.totalPosError
          print "Average delta error: "+str(avgDError)
        
        print "Lowest error: "+str(self.lowestError)
      
    elif mode ==  "random":
      print "Mode:\tRandom"
      #prog = randomizer.randomizer(nRuns, 1, self.seqFile.name, self.bindingFile.name, saveRuns)
    """
    # initializing the NUPACK object in silent mode for the calculation of the concentrations in the first steady state
    self.concentrations1 = nupack.NUPACK(prefix=os.getcwd()+"/Concentrations/LowestError",
                                paramFile=os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))+"/parameters.txt",
                                silent=True, degenerate=False, material="dna", dangles="all", ordered=True, timeonly=False, mfe=True, quiet=True, sort=1)
    # initializing the NUPACK object in silent mode for the calculation of the concentrations in the second steady state
    self.concentrations2 = nupack.NUPACK(prefix=os.getcwd()+"/Concentrations/LowestErro",
                                paramFile=os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))+"/parameters.txt",
                                silent=True, degenerate=False, material="dna", dangles="all", ordered=True, timeonly=False, mfe=True, quiet=True, sort=1)
    # Calculate all complexes in solution at first steady state and calculate concentrations
    self.concentrations1.runComplexes()
    self.concentrations1.runConcentrations()
    complexIndices1, mult1, conc1, ener1 = self.concentrations1.readConcOutput()
    
    # Save information in file
    self._writeConcentrationFile(complexIndices1, mult1, conc1, ener1, 'ConcNuPack1')

    # Calculate all complexes in solution at second steady state and calculate concentrations
    self.concentrations2.runComplexes()
    self.concentrations2.runConcentrations()
    complexIndices2, mult2, conc2, ener2 = self.concentrations2.readConcOutput()

    # Save information in file
    self._writeConcentrationFile(complexIndices2, mult2, conc2, ener2, 'ConcNuPack2')
  
  # function for the SA loop
  def _mainloop(self, nRuns, T_start, T_end, lambdaFactor, screenOut=False):

    # opening output files
    gcFile = open("GC_Content15_seed1_Temp15000_1e6.txt", "w")
    errFile = open("Errors15_seed1_Temp15000_1e6.txt", "w")
    feFile = open("Free Energies15_seed1_Temp15000_1e6.txt", "w")

    if screenOut:
      self.printSequences()
    T = T_start
    # T = 10000000

    # initializing temperature decrement factor
    factor = (float(T_end)/float(T))**(float(1)/float(nRuns))
    # factor = 1
    
    run = 0
    
    # while the end temperature is not reached
    while T > T_end:
      run += 1


      # make random pointmutation and adjust this in the other sequences
      self._processMutation()	
      if lambdaFactor>0:
        # This section corresponds to the execution of the criton check class as a parallel process to the error calculation
        #####################################################################
        # initialize parallel process that performs the criton check. Criton check returns the weight value that should be added to error function
        # self.critonsChck = critonsCheck.CritonsCheck(run,nRuns,lambdaFactor,"/home/hroekel/DNAStrandDesign_MM/TUPack_Upgrade_Zandra/model/CritonsInfoFile" ,lc=3,strandsObject=self.Strands)     
        # self.critonsChck = critonsCheck.CritonsCheck(run,nRuns,lambdaFactor,"/home/hroekel/DNAStrandDesign_MM/TUPack_Upgrade_Zandra/model/CritonsDataFile_Lambda15_seed1_Temp15000_1e6_","/home/hroekel/DNAStrandDesign_MM/TUPack_Upgrade_Zandra/model/CritonsInfoFile_Lambda15_seed1_Temp15000_1e6_",lc=3,strandsObject=self.Strands,excTemp1=['GACTC','GACTC'], excTemp2=['CTCGA','GACT','GACTC'], excTemp3=['CTCGA','GACTC'])
        self.critonsChck = critonsCheck.CritonsCheck(run,nRuns,lambdaFactor,"/home/hroekel/DNAstranddesign_MM/TUPack_Upgrade_Zandra15_1/model/CritonsDataFile_Lambda15_seed1_Temp15000_1e6_","/home/hroekel/DNAstranddesign_MM/TUPack_Upgrade_Zandra15_1/model/CritonsInfoFile_Lambda15_seed1_Temp15000_1e6_",lc=3,strandsObject=self.Strands)
        print "Starting pp with ", self.job_server.get_ncpus(), " workers"
        job1 = self.job_server.submit(self.critonsChck._checkCritons,())
        # Retrieves the result calculated by job1
        # The value of job1() is the same as critonsChck._checkCritons() 
        # If the job has not been finished yet, execution will
        # wait here until result is available
        self.critonWeight = job1()	
        #####################################################################

      # run complexes, read the output
      NpOut = self._runNuPack()
      
      """structpenalty=0
      for elem in mfestruct:
        b=float(elem)
        structpenalty+=b**2"""
        
      # calculate error
      #print 'NpOut: ',NpOut
      self.gibbsError = self._fillErrorMatrix(NpOut)
      newError = self.gibbsError + self.critonWeight

      try:
        dError = newError - Error
      except NameError:
        Error = newError + 1
        self.lowestError = Error
        dError = newError - Error
        
      if dError > 0:
        self.totalDError += dError
        self.totalPosError += 1

      # if the new error is lower than the current error, accept mutation
      if newError <= Error:
        Error = newError
        fe = NpOut

        # save the lowest error
        if self.lowestError > Error:
          self.lowestError = Error
          self._saveRun(os.getcwd()+"/Concentrations/LowestError3",1)
          self._saveRun(os.getcwd()+"/Concentrations/LowestError4",1)

      # if the new error is not lower, accept the mutation with a small chance
      elif random.random() < math.exp(-dError/T):
        Error = newError
        fe = NpOut

      # reject and rollback mutation
      else:
        self._rollback()

      err = [Error, dError]

      gc = []
      for strand in self.Strands:
        gc += [strand.getGC()]

      # write the output files
      if run == 1:
        self._writeFeFile(fe, feFile, True)
        self._writeGcFile(gc, gcFile, True)
        self._writeErrFile(newError, err, errFile, True)
      elif (run%100) == 0:
        self._writeFeFile(fe, feFile, False)
      elif (run%10) == 0:
        self._writeGcFile(gc, gcFile, False)
        self._writeErrFile(newError, err, errFile, False)

      # if the run is in the save list, run an extra time on another file name to prevent overwriting the files
      if run in self.saveRuns:
        self._saveRun("run"+str(run),0)

      if run % 1 == 0:
        print "Done "+str(run)+" of "+str(nRuns)+"\tCurrent Error: "+str(Error)+"\tTemperature: "+str(T)

      # decrease the temperature
      T = T * factor

      if screenOut:
        self.printSequences(True)

    self._saveRun(os.getcwd()+"/Concentrations/Final1",1)
    self._saveRun(os.getcwd()+"/Concentrations/Final2",1)

    # close the output files
    gcFile.close()
    feFile.close()
    errFile.close()


  # function for saving NuPack output files
  def _saveRun(self, title, condition):
    if condition == 0:
      os.system("cp cxFile.in " +title+".in")
      os.system("cp cxFile.cx " +title+ ".cx")
      os.system("cp cxFile.ocx " +title+ ".ocx")
      os.system("cp cxFile.ocx-mfe " +title+".ocx-mfe")
    else:
      self._writeNpInput(title, self.maxCxSize, 1)

  # extra function for printing the sequences, mostly for debugging
  def printSequences(self, comp=False):
    for seq in self.Strands:
      if seq.Type == "binding" or seq.Type == "inhibitor":
        if comp:
          seq.Sequence.reverse()
          print seq.Sequence + " (REVERSED)"
          seq.Sequence.reverse()
        else:
          print seq.Sequence
      else:
        print seq.Sequence
    print "\n=================\n"

  # function for writing the output file containing the GC content
  def _writeGcFile(self, gcList, gcFile, printHeader=True):

    # if the header should be written
    if printHeader:
      header = ""
      for strand in self.Strands:
        header += strand.Name + "\t"

      gcFile.write(header.rstrip("\t") + "\n")
      
    gcData = ""
    for strandGC in gcList:
      gcData += str(strandGC) + "\t"
    gcFile.write(gcData.rstrip("\t") + "\n")

  # function for writing the output file containing the errors
  def _writeErrFile(self, newError, errList, errFile, printHeader=True):
    if printHeader:
      errFile.write("#Error:\tdError:\n")
      errFile.write("#Final error\t\delta Final error\t\Gibbs contribution to error function\t\Critons contribution to error function\t\Total error\t\Critons to Total error ratio\n")
    
    errFile.write(str(errList[0])+"\t"+str(errList[1])+"\t"+str(self.gibbsError)+"\t"+str(self.critonWeight)+"\t"+str(newError)+"\t"+str(self.critonWeight/newError)+"\n")

  # function for reading and writing NuPack complexes
  def _runNuPack(self):
    
    # write the nupack input files '.in' and '.list'
    self._writeNpInput("cxFile", self.maxCxSize, 0)
        
    # run the nupack binary complexes
    self.nupack.runComplexes()      

    # add the output to the free energy list
    NpOut = self.nupack.readLastOutput()

    return NpOut

  # makes a random mutation in a sequence, and adjusts this in dependant sequences
  def _processMutation(self):

    # make a backup of the current sequences
    self._backup()

    # make the random mutation
    mutated = self._makeRandomMutation()

    # make a list with all sequences
    seqList = []
    for strand in self.Strands:
      seqList.append(strand.Name)

    # list containing the adjust order
    self.adjustOrder = []

    # make the order in which the sequences should be made complementary
    self._makeAdjustOrder(mutated, seqList)

    # make the adjustments in the sequences, and an extra function for the inhibitors
    self._adjustSequenceSet()
    self._adjustInhibitors()

    #self._checkComplement() 

  # rolls back the point mutation and restores the sequences to the former state
  def _rollback(self):
    self.Strands = []
    self.Duplexes = []
    self.Inhibitors = []

    # appending new lists with deepcopies of the original strands, duplexes and inhibitors
    for strand in self.StrandBackup:
      original = copy.deepcopy(strand)      
      self.Strands += [original]

    for duplex in self.DuplexBackup:      
      s1 = duplex.Strand1
      s2 = duplex.Strand2

      for strand in self.Strands:
        if strand.Name == s1.Name:
          strand1 = strand
        if strand.Name == s2.Name:
          strand2 = strand

      original = copy.deepcopy(duplex)
      original.updateStrands(strand1, strand2)
      
      self.Duplexes += [original]
      
    for inh in self.InhBackup:
      s1 = inh.Strand1
      s2 = inh.Strand2

      for strand in self.Strands:
        if strand.Name == s1.Name:
          strand1 = strand
        if strand.Name == s2.Name:
          strand2 = strand

      original = copy.deepcopy(inh)
      original.updateStrands(strand1, strand2)

      self.Inhibitors += [original]

  # makes a backup of the current sequences
  def _backup(self):
    self.StrandBackup = []
    self.DuplexBackup = []
    self.InhBackup = []

    for strand in self.Strands:
      backup = copy.deepcopy(strand)
      self.StrandBackup += [backup]

    for duplex in self.Duplexes:      
      s1 = duplex.Strand1
      s2 = duplex.Strand2

      for strand in self.StrandBackup:
        if strand.Name == s1.Name:
          strand1 = strand
        if strand.Name == s2.Name:
          strand2 = strand

      backup = copy.deepcopy(duplex)
      backup.updateStrands(strand1, strand2)
      
      self.DuplexBackup += [backup]

    for inh in self.Inhibitors:
      s1 = inh.Strand1
      s2 = inh.Strand2

      for strand in self.StrandBackup:
        if strand.Name == s1.Name:
          strand1 = strand
        if strand.Name == s2.Name:
          strand2 = strand

      backup = copy.deepcopy(inh)
      backup.updateStrands(strand1, strand2)
      
      self.InhBackup.append(backup)

  # makes an error matrix containing all individual errors of the duplexes
  # secondary structure penalty is not taken into account
  def _fillErrorMatrix(self, output, screenOut=False):
    selfError = {}
    errorMatrix = {}
    
    for duplex in self.Duplexes + self.Inhibitors:
      if duplex.Strand1 not in errorMatrix.keys():
        errorMatrix[duplex.Strand1.Name] = {}
      if duplex.Strand2 not in errorMatrix.keys():
        errorMatrix[duplex.Strand2.Name] = {}
    
    for row in output:
      #print row
      # if the permutation is done with one strand
      if len(row) == 2:
        for strand in self.Strands:
          if strand.ID == int(row[0]):

            # define individual errors
            # print "Single species energy: "+str(float(row[1]))
            selfError[strand.Name] = self._calcError(10*float(row[1]), None)

      # if the permutation is done with two strands
      if len(row) == 3:
        for strand in self.Strands:
          if int(row[0]) == strand.ID:
            s1 = strand.Name

          if int(row[1]) == strand.ID:
            s2 = strand.Name

        # print 'check: ',s1, s2, output[row[0]-1][1], output[row[1]-1][1]
        energy = float(row[2])-float(output[row[0]-1][1])-float(output[row[1]-1][1])
        # print 'Energy: ',energy, float(row[2])
        
        duplex = None
        for dup in self.Duplexes + self.Inhibitors:
          if (dup.Strand1.Name == s1 and dup.Strand2.Name == s2) or (dup.Strand1.Name == s2 and dup.Strand2.Name == s1):
            duplex = dup           

        # define individual errors
        errorMatrix[s1][s2] = self._calcError(energy, duplex)
        
        # print errorMatrix

    # calculate total error
    totalScore = float()
    for key in errorMatrix.keys():
      row = errorMatrix[key]
      for cell in row:
        totalScore += row[cell]

    # if screen output is allowed
    if screenOut:
      header = "\t|\t"
      for strand in self.Strands:
        header += strand.Name + "\t\t"
      print header
      print "="*(16*len(self.Strands)+16)

      count = 0
      for strand1 in self.Strands:
        for row in errorMatrix.keys():
          if strand1.Name == row:
            out = ""
            for strand2 in self.Strands:
              for col in errorMatrix[row].keys():
                if strand2.Name == col:
                  out += str(errorMatrix[row][col]) +"\t"
            print row+"\t|\t"+("\t\t"*count)+out
        count += 1

      print "="*len(header)*5
      print "Total error:\t"+str(totalScore)+"\n"

      for key in selfError.keys():
        for strand in self.Strands:
          if key == strand.Name:
            print strand.Name + "\t" + str(selfError[key])
    for key in selfError.keys():
        totalScore+=selfError[key]
    return totalScore                  

  # defines individual errors, this function is A.K.A. the cost function
  def _calcError(self, predictedE, duplex):
    if duplex != None:
      target = duplex.TargetEnergy
      rob = duplex.Robustness
    else:
      target = 0.0
      rob = 1.0
    # print "Prediction: "+str(float(predictedE))+" Target: "+str(float(target))+" Difference: "+str(float(predictedE)-float(target))
    score = (1/rob)*(float(predictedE)-float(target))
    """print "Old error = "+str(float(predictedE)-float(target))
    print "New error = "+str(score)"""
    score = score**2

    return score

  # makes the strands and inhibitors complementary in the right order
  def _adjustInhibitors(self):
    for duplex in self.Inhibitors:
      if duplex.RelationType == "inhibiting":
        duplex.adjustSequences(duplex.Strand1.Name)

    for duplex in self.Inhibitors:
      if duplex.RelationType == "binding":
        duplex.adjustSequences(duplex.Strand2.Name)

  # makes all strands complementary    
  def _adjustSequenceSet(self):

    # if there is a specified order
    if self.adjustOrder != []:
      for order in self.adjustOrder:
        dupIndex = 0
        for dup in self.Duplexes:
          if (dup.Strand1.Name == order[0] and dup.Strand2.Name == order[1]):
            dup.adjustSequences(dup.Strand1.Name)
          if (dup.Strand1.Name == order[1] and dup.Strand2.Name == order[0]):
            dup.adjustSequences(dup.Strand2.Name)

    # if the order is not important
    else:
      for duplex in self.Duplexes:
        if duplex.Strand1.Mutated and duplex.Strand2.Mutated:
          pass
          
        if duplex.Strand1.Mutated and not duplex.Strand2.Mutated:
          duplex.adjustSequences(duplex.Strand1.Name)

        if not duplex.Strand1.Mutated and duplex.Strand2.Mutated:
          duplex.adjustSequences(duplex.Strand2.Name)

        if not duplex.Strand1.Mutated and not duplex.Strand2.Mutated:
          duplex.adjustSequences(duplex.Strand1.Name)

      for strand in self.Strands:
        strand.Mutated = False

  # makes a random mutation on a random template strand
  def _makeRandomMutation(self):
    templates = []

    index = []
    for strand in self.Strands:

      # only templates can be mutated
      if strand.Type == "template":
        templates += [strand]

    # pick a random template strand
    template = random.choice(templates)

    # make 1 random mutation
    template.randomlyMutate()

    return template

  # makes the order in which the sequences should be adjusted
  def _makeAdjustOrder(self, baseStrand, seqList):
    if baseStrand.Name in seqList:
      seqList.remove(baseStrand.Name)

    for dup in self.Duplexes:
      Continue = False
      if dup.Strand1.Name == baseStrand.Name:
        compStrand = dup.Strand2
        Continue = True
      if dup.Strand2.Name == baseStrand.Name:
        compStrand = dup.Strand1
        Continue = True

      if Continue:
        if compStrand.Name in seqList:
          self.adjustOrder += [[baseStrand.Name, compStrand.Name]]

          # call itself to check if sequences depend on this one
          self._makeAdjustOrder(compStrand, seqList)

  # function for reading the sequence input
  def _readSeqFile(self, seqFile):

    try:
      
      # reading the sequences file
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

  # reads duplex input
  def _readBindingFile(self, bindingFile):
    try:
      # for each line in the binding file
      for line in bindingFile:
        line = line.strip("\n").strip("\r").split("\t")

        strandFound = False
        for strand in self.Strands:
          if strand.Name == line[0]:
            Strand_1 = strand
            strandFound = True
          if strand.Name == line[1]:
            Strand_2 = strand
            strandFound = True

        if strandFound == False:
          sys.exit(line[0]+" or "+line[1]+" is not defined in sequence file.")

        # if type is inhibitor, define in the duplex object whether it is a normal binding or an inhibiting duplex
        if Strand_2.Type == "inhibitor":
          if line[4] == "B":
            duplex = Strand.StrandRelation(Strand_1, Strand_2, "binding")
            duplex.setImmutable()
          elif line[4] == "I":
            duplex = Strand.StrandRelation(Strand_1, Strand_2, "inhibiting")
          else:
            sys.exit("Inhibitor binding types must be specified.")

          # define target energy, binding structure
          duplex.TargetEnergy = float(line[3])
          duplex.Robustness = float(line[5])
          duplex.defineBindingStructure(line[2], ["(", ")", "#"])
          self.Inhibitors += [duplex]

        else:            
          duplex = Strand.StrandRelation(Strand_1, Strand_2)
          duplex.TargetEnergy = float(line[3])
          duplex.Robustness = float(line[5])
          duplex.defineBindingStructure(line[2], ["(", ")"])

          self.Duplexes += [duplex]
        
    except IndexError:
      sys.exit("The binding file has a wrong format. Check the user manual for the right file formats.")

  # function for writing the free energies to a file.
  def _writeFeFile(self, feList, feFile, printHeader=True):

    header = ""
    freeE = ""
    for j in range(len(feList)):
      
      # if the permutation is done with one strand
      if len(feList[j]) == 2:
        for strand in self.Strands:
          if strand.ID == feList[j][0]:
            if printHeader:
              header += strand.Name+"\t" 

            # add free energies to list
            freeE += str(feList[j][1]) +"\t"

      # if the permutation is done with two strands
      if len(feList[j]) == 3:
        for duplex in self.Duplexes + self.Inhibitors:
          if (duplex.Strand1.ID == feList[j][0] and duplex.Strand2.ID == feList[j][1]) or (duplex.Strand1.ID == feList[j][1] and duplex.Strand2.ID == feList[j][0]):
            if printHeader:
              header += duplex.Name+"\t"
              
            freeE += str(feList[j][2]) +"\t"

    if printHeader:
      feFile.write(header.rstrip("\t") + "\n")
    feFile.write(freeE.rstrip("\t") + "\n")

  # function for writing the NuPack input files
  def _writeNpInput(self, fileName, maxCxSize, condition):
    # open the input files for NuPack input
    cxFile = open(fileName+".in", "w")

    # define the number of sequences and the maximum complex size
    nSeqs = len(self.Strands)

    # write the '.in' file
    if condition == 0:
      cxFile.write(str(nSeqs)+"\n")
      for strand in self.Strands:      
        cxFile.write(str(strand.Sequence)+"\n")
      cxFile.write(str(maxCxSize))
    else:
      cxFile.write(str(nSeqs+3)+"\n")
      for strand in self.Strands:      
        cxFile.write(str(strand.Sequence)+"\n")
      cxFile.write(str(self.Strands[4].Sequence)+str(self.Strands[5].Sequence+"\n"))
      cxFile.write(str(self.Strands[5].Sequence)+str(self.Strands[6].Sequence+"\n"))
      cxFile.write(str(self.Strands[5].Sequence)+str(self.Strands[7].Sequence+"\n"))
      cxFile.write(str(3))

    # close the '.in' file
    cxFile.close()

    if maxCxSize == 1 and condition == 0:
      listFile = open(fileName+".list", "w")
      output = ""
      
      # write the '.list' file
      # write normal duplexes
      for duplex in self.Duplexes:
        id1 = duplex.Strand1.ID
        id2 = duplex.Strand2.ID
        
        output += str(id1)+" "+str(id2)+"\n"
        
      if self.Inhibitors != []:
        for duplex in self.Inhibitors:
          id1 = duplex.Strand1.ID
          id2 = duplex.Strand2.ID
          
          output += str(id1)+" "+str(id2)+"\n"
          
      output = output.rstrip("\n")
      listFile.write(output)
      # close the '.list' file
      listFile.close()

  # Relevant information will be put in file, i.e. species names and concentrations.
  def _writeConcentrationFile(self, ind, mult, conc, ener, fileName):
    ID = {1:'T1', 2:'T2', 3:'T3', 4: 'F', 5:'U', 6:'a', 7:'Y', 8:'inh', 9:'Ua', 10:'aY', 11:'a Inh'}
    concFile = open(fileName+".txt", "w")
    concFile.write('Species\t\tConcentration [nM]\t\tFree Energy [kcal/mol]\n')
    for i in range(len(ind)):
      for j in range(len(ind[i])):
        if j == 0:
          namestring=ID[ind[i][j]]+(mult[i][j]-1)*('.'+ID[ind[i][j]])
        else:
          namestring=namestring+mult[i][j]*('.'+ID[ind[i][j]])
      concFile.write(namestring+'\t\t'+str(conc[i]*10**9)+'\t\t'+str(ener[i])+'\n')
    concFile.close()
        
      
  # function for popping up a file chooser
  def _chooseFile(self, Title):

    # initialize the Tk object
    root = Tkinter.Tk()

    # withdraw the main window
    root.withdraw()

    # open file chooser on the current directory
    File = tkFileDialog.askopenfile(parent=root, mode='r', title=Title, initialdir=os.getcwd())

    # exit the windows
    root.quit()

    # return files
    return File
# saveRuns=[1,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000]
prog = calcConc(100,15000,10,0,maxCxSize=2,
              seqFile='/home/hroekel/DNAstranddesign_MM/TUPack_Upgrade_Zandra15_1/model/sequences.txt',
              bindingFile='/home/hroekel/DNAstranddesign_MM/TUPack_Upgrade_Zandra15_1/model/binding.txt',
              saveRuns=range(1,200001,1000),
              mode='SA')
