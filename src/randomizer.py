# randomizer class
# This program is part of the TUPack software package and makes use of NuPack 3.0. NuPack can be found on <www.nupack.org>g>
# This class creates an object performing random sequence mutations to outline the spread of the free energy predictions.
# For more info about the program, please read the user manual.
#
# written by Sander Rodenburg
# Eindhoven University of Technology
# 2012

# importing the NUPACK object, for running complexes and for reading the NuPack output files
import nupack

# importing the Strand library, containing the objects Strand and StrandRelation for implementing the DNA network
import Strand

# importing other libraries
import random
import math
import time
import os
import sys
import Tkinter
import tkFileDialog
import copy

class randomizer:
  
  def __init__(self, nRuns, maxCxSize=1, seqFile=None, bindingFile=None, saveRuns=[]):

    # defining seed for random
    random.seed(time.time())
    #random.seed(1)

    # lists holding the Strands, duplexes and Inhibitors
    self.Strands = []
    self.Duplexes = []
    self.Inhibitors = []

    # list for keeping the order in which the sequences should be adjusted
    self.adjustOrder = []

    # NuPack maximum complex size
    self.maxCxSize = maxCxSize

    # runs in which the program should make a distinct save
    self.saveRuns = saveRuns

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
                                silent=True)

    # reading input files for initial sequence generations   
    self._readSeqFile(self.seqFile)
    self._readBindingFile(self.bindingFile)

    # close the files
    self.seqFile.close()
    self.bindingFile.close()

    # adjust initial sequences
    self._adjustSequenceSet()
    self._adjustInhibitors()

    # run the mainloop for simulated annealing and sequence mutation
    self._mainloop(nRuns, False)

    print "Done."    
    print "Files stored in " +os.getcwd()

  def _mainloop(self, nRuns, screenOut=False):
    gcFile = open("GC_Content.txt", "w")
    feFile = open("Free Energies.txt", "w")

    if screenOut:
      self.printSequences()
    
    run = 0
    
    # for each run
    for Run in range(nRuns):
      run += 1
      
      self._processMutation()

      NpOut = self._runNuPack()
      fe = NpOut

      gc = []
      for strand in self.Strands:
        gc += [strand.getGC()]

      if run == 1:
        self._writeFeFile(fe, feFile, True)
        self._writeGcFile(gc, gcFile, True)
      else:
        self._writeFeFile(fe, feFile, False)
        self._writeGcFile(gc, gcFile, False)      

      # if the run is in the save list, run an extra time on another file name to prevent overwriting the files
      if run in self.saveRuns:
        self._saveRun("run"+str(run))

      if run % 1000 == 0:
        print "Done "+str(run)+" of "+str(nRuns)

      if screenOut:
        self.printSequences(True)

    #self.printSequences(True)   

    gcFile.close()
    feFile.close()

  # copy the output files, to prevent overwriting
  def _saveRun(self, title):
    os.system("cp cxFile.in " +title+".in")
    os.system("cp cxFile.cx " +title+ ".cx")
    os.system("cp cxFile.ocx-mfe " +title+".ocx-mfe") 

  # print the sequences, reversed or normal
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

  # write the GC file
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
  def _writeErrFile(self, errList, errFile, printHeader=True):
    if printHeader:
      errFile.write("Error:\tdError:\n")
    
    errFile.write(str(errList[0])+"\t"+str(errList[1])+"\n")

  # function for reading and writing NuPack complexes
  def _runNuPack(self):
    # write the nupack input files '.in' and '.list'
    self._writeNpInput("cxFile", self.maxCxSize)
        
    # run the nupack binary complexes
    self.nupack.runComplexes()      

    # add the output to the free energy list
    NpOut = self.nupack.readLastOutput()

    return NpOut

  # makes a random mutation in a sequence, and adjusts this in dependant sequences
  def _processMutation(self):

    # make a random mutation
    mutated = self._makeRandomMutation()
    
    seqList = []
    for strand in self.Strands:
      seqList.append(strand.Name)

    # make the order in which the sequences should be made complementary
    self.adjustOrder = []
    self._makeAdjustOrder(mutated, seqList)

    # make the sequences and the inhibitors complementary
    self._adjustSequenceSet()
    self._adjustInhibitors()

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

  # make a random mutation
  def _makeRandomMutation(self):
    templates = []

    index = []
    for strand in self.Strands:

      # only templates can be mutated
      if strand.Type == "template":
        templates += [strand]

    # pick a random template strand
    template = random.choice(templates)

    # randomize the whole sequence
    template.randomize()

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

        if line[2] == "T":
          newStrand.Type = "template"
        elif line[2] == "B":
          newStrand.Type = "binding"
        elif line[2] == "I":
          newStrand.Type = "inhibitor"
        ID += 1

        #print newStrand.Name, newStrand.Type

    except IndexError:
      sys.exit("The sequence file has a wrong format. Check the user manual for the right file formats.")

  # reads the binding file
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

          duplex.TargetEnergy = float(line[3])
          duplex.defineBindingStructure(line[2], ["(", ")", "#"])
          self.Inhibitors += [duplex]

        else:            
          duplex = Strand.StrandRelation(Strand_1, Strand_2)
          duplex.TargetEnergy = float(line[3])
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
      if len(feList[j]) == 3:
        for strand in self.Strands:
          if strand.ID == feList[j][0]:
            if printHeader:
              header += strand.Name+"\t"

            # add free energies to list
            freeE += str(feList[j][1]) +"\t"

      # if the permutation is done with two strands
      if len(feList[j]) == 4:
        for duplex in self.Duplexes + self.Inhibitors:
          if (duplex.Strand1.ID == feList[j][0] and duplex.Strand2.ID == feList[j][1]) or (duplex.Strand1.ID == feList[j][1] and duplex.Strand2.ID == feList[j][0]):
            if printHeader:
              header += duplex.Name+"\t"
              
            freeE += str(feList[j][2]) +"\t"

    if printHeader:
      feFile.write(header.rstrip("\t") + "\n")
    feFile.write(freeE.rstrip("\t") + "\n")

  # function for writing the NuPack input files
  def _writeNpInput(self, fileName, maxCxSize):
    # open the input files for NuPack input
    cxFile = open(fileName+".in", "w")

    if maxCxSize == 1:
      listFile = open(fileName+".list", "w")

    # define the number of sequences and the maximum complex size
    nSeqs = len(self.Strands)

    # write the '.in' file
    cxFile.write(str(nSeqs)+"\n")
    
    for strand in self.Strands:      
      cxFile.write(str(strand.Sequence)+"\n")      
    cxFile.write(str(maxCxSize))

    # close the '.in' file
    cxFile.close()

    if maxCxSize == 1:
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
