# Strand class and StrandRelation class
# This program is part of the TUPack software package and makes use of NuPack 3.0. NuPack can be found on <www.nupack.org>>
# This file specifies an object representing a DNA strand, and an object which represents a relationship between strands.
# For more info about the program, please read the user manual.
#
# written by Sander Rodenburg
# Eindhoven University of Technology
# 2012

# importing libraries
from Bio.Seq import MutableSeq
from Bio.Seq import IUPAC
from Bio.SeqUtils import GC

import random
import time

# Strand object
class Strand:
  
  def __init__(self, seq, name, ID, strandType=None):
    random.seed(time.time())
    #random.seed(1)
    
    self.Name = name              # sequence name
    self.ID = ID                  # the sequence ID. Used for identification
    self.definedN = False         # are all undefined bases (N) defined?
    self.Mutated = False          # is this sequence already mutated?
    self.Immutables = []          # immutable base locations

    # sequence cannot contain any other symbols than ATGCN
    for base in seq:
      if base not in "ATGCN":
        raise ValueError("All bases in the sequence must be A/T/G/C or N.")
    self.Sequence = MutableSeq(seq, IUPAC.unambiguous_dna)  

    # strand type must be declared from the three choices below, or None
    if strandType not in ["template", "binding", "inhibitor"] and strandType != None:
      raise ValueError("'Type' parameter must be 'template', 'binding', 'inhibitor'.")
    self.Type = strandType

  # defines undefined bases
  def defineN(self):
    i = 0
    for base in self.Sequence:

      # randomly define bases
      if base == "N":
        self.Sequence[i] = random.choice("ATGC")

      # if bases are already defined, set them immutable
      if base != "N":
        self.Immutables += [i]
      i += 1

    self.definedN = True

  # randomizes the whole sequence, except for the immutable bases
  def randomize(self):
    for baseIndex in range(len(self.Sequence)):
      if baseIndex not in self.Immutables:
        self.Sequence[baseIndex] = random.choice("ATGC")

  # makes one random mutation
  def randomlyMutate(self):
    if self.definedN == False:
      self.defineN()

    # choose a random base location for the mutation
    baseIndex = random.randrange(0, len(self.Sequence)-1)

    # while this location is immutable, choose another one
    while baseIndex in self.Immutables:
      baseIndex = random.randrange(0, len(self.Sequence)-1)

    originalBase = self.Sequence[baseIndex]
    self.Sequence[baseIndex] = random.choice("ATGC")

    # go on mutating until the base has another value
    while self.Sequence[baseIndex] == originalBase:
      self.Sequence[baseIndex] = random.choice("ATGC")

    # now the sequence is mutated
    self.Mutated = True

  # returns the GC content of itself
  def getGC(self):
    gc = GC(self.Sequence)
    return gc

# StrandRelation class
class StrandRelation:
  
  def __init__(self, strand_1, strand_2, relationType="binding"):
    
    self.Name = str(strand_1.Name)+"-"+str(strand_2.Name) # the name is equal to the two strand names separated by "-"
    self.Strand1 = strand_1                               # define the participating strands
    self.Strand2 = strand_2                         
    self.Bindings = {}                                    # dictionary which contains the corresponding base pairs of the two strands 
    self.Incomplements = {}                               # contains the incomplementary base pair locations of the two strands
    self.TargetEnergy = None                              # the target Gibbs free energy of the duplex

    # define the type
    if relationType not in ["binding", "inhibiting"]:
      raise ValueError("Parameter 'relationType' must be 'binding' or 'inhibiting'.")
    self.RelationType = relationType

    # complementary bases
    self.complement = {"A":"T", "T":"A", "G":"C", "C":"G"}

  # makes the two Strand objects complementary to each other on the basis of one Strand
  def adjustSequences(self, base_strand_name):
    immutables = []
    incomp = {}

    # if not all undefined bases are defined, define them
    if self.Strand1.definedN == False:
      self.Strand1.defineN()
    if self.Strand2.definedN == False:
      self.Strand2.defineN()

    if self.Bindings == None:
      raise RuntimeError("Binding structure is not defined.")

    # checks which strand is the base strand, and from there get the binding locations etc.
    if base_strand_name == self.Strand1.Name:
      baseStrand = self.Strand1.Sequence
      compStrand = self.Strand2.Sequence
      bindings = self.Bindings
      incomp = self.Incomplements
      immutables = self.Strand1.Immutables
        
    elif base_strand_name == self.Strand2.Name:
      baseStrand = self.Strand2.Sequence
      compStrand = self.Strand1.Sequence

      # swap dictionary keys for values, because the keys are based on strand 1
      bindings = self._swapDict(self.Bindings)
      incomp = self._swapDict(self.Incomplements)
      immutables = self.Strand2.Immutables

    else:
      raise ValueError("Base strand name not specified!")

    # make the mutations
    baseIndex = 0
    for base in baseStrand:
      if baseIndex in bindings.keys():
        compStrand[bindings[baseIndex]] = self.complement[baseStrand[baseIndex]]
      if baseIndex in incomp.keys():
        compStrand[incomp[baseIndex]] = baseStrand[baseIndex]
      baseIndex += 1

    self.Strand1.Mutated = True
    self.Strand2.Mutated = True

  # defines the binding locations of complementary bases, with the input of a DPP notation
  def defineBindingStructure(self, dppNotation, symbols, oppositeDir=True):
    self.Structure = dppNotation
    
    structures = dppNotation.split("+")
    str1 = structures[0]
    str2 = structures[1]
    
    # if one symbol is defined, take that symbol twice for reading the structure index    
    if len(symbols) == 1:
      sym1 = symbols[0]
      sym2 = symbols[0]
      incompSym = None
    elif len(symbols) == 2:
      sym1 = symbols[0]
      sym2 = symbols[1]
      incompSym = None
    elif len(symbols) == 3:
      sym1 = symbols[0]
      sym2 = symbols[1]
      incompSym = symbols[2]
      
    else:
      raise IOError("Symbols must have one or two elements.")

    # check which opening symbols are with which closing symbols in the DPP notation
    if sym1 != None and sym2 != None:
      indices = {}
      opencount = 0
      str1index = 0
      for s1 in str1:
        if s1 == sym1:
          opencount += 1
          closecount = 0
          if oppositeDir:
            str2index = len(str2)-1
            for s2 in reversed(str2):
              if s2 == sym2:
                closecount += 1
                if opencount == closecount:
                  indices[str1index]=str2index
              str2index -= 1
          else:
            str2index = 0
            for s2 in str2:
              if s2 == sym2:
                closecount += 1
                if opencount == closecount:
                  indices[str1index]=str2index
              str2index += i
        str1index += 1

      self.Bindings = indices

    # if there are some incomplementary bases, do the same as before, but with other symbols as opening/closing symbols
    if incompSym != None:      
      indices = {}
      opencount = 0
      str1index = 0
      for s1 in str1:
        if s1 == incompSym:
          opencount += 1
          closecount = 0
          if oppositeDir:
            str2index = len(str2)-1
            for s2 in reversed(str2):
              if s2 == incompSym:
                closecount += 1
                if opencount == closecount:
                  indices[str1index]=str2index
              str2index -= 1
          else:
            str2index = 0
            for s2 in str2:
              if s2 == incompSym:
                closecount += 1
                if opencount == closecount:
                  indices[str1index]=str2index
              str2index += 1
        str1index += 1
      self.Incomplements = indices

  # swap dictionary keys for values
  def _swapDict(self, dictionary):
    swap = {}
    for key, val in dictionary.iteritems():
      swap[val] = key

    return swap

  # update the strands in the relation
  def updateStrands(self, strand_1, strand_2):
    self.Strand1 = strand_1
    self.Strand2 = strand_2

  # set a strand to immutable
  def setImmutable(self):
    for key in self.Bindings.keys():
      Strand1.Immutables += [key]
