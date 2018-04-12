#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cas_object.py hold Cas info as part of ExcisionFinder. Written in Python v 3.6.1.
Kathleen Keough and Michael Olvera 2018.

cas_object.py will do two things:
	 + hold the Cas9 object info (one object per Cas enzyme) and 
	 will handle input error checking.
	 + contain a function to create the objects from an external list
	 (I was thinking just a text file, so if users want to add more
	 enzymes in the future they can)
"""
import os, sys
CAS_PATH=os.path.dirname(os.path.realpath(sys.argv[0])) + '/CAS_LIST.txt'
IUPAC={'Y':'[CT]','R':'[AG]','W':'[AT]',
	'S':'[GC]','K':'[TG]','M':'[CA]','D':'[AGT]',
	'V':'[ACG]','H':'[ACT]','B':'[CGT]','N':'[ACGT]',
	'A':'A','C':'C','G':'G','T':'T'}


class Cas(object):
	"""
	This is the actual Cas object, holding the enzyme name, primeness, and PAM.
	The regex version of the PAM can be outputed. The input PAM is assumed to be
	in the 5'->3' orientation, and can be converted to its reverse orientation.

	"""
	def __init__(self, name, forwardPam, primeness):
		self.name = name
		self.forwardPam = forwardPam
		self.primeness = primeness

	def __str__(self):
		return f"<Cas Object {self.name}, with {self.primeness} PAM:{self.forwardPam}>"


	def forwardPam_regex(self):
		"""
		Returns the regex of the forward PAM. 
		"""
		return ''.join([IUPAC[c] for c in self.forwardPam])
	def reversePam_regex(self):
		"""
		Returns the regex of the reverse PAM. 
		"""
		pam = self.getReversePam()
		return ''.join([IUPAC[c] for c in pam])
	def getReversePam(self):
		"""
		Returns the reverse PAM of the input.
		"""
		watson ='ACGTYRSWKMBDHVN'
		crick = 'TGCARYSWMKVHDBN'
		return self.forwardPam[::-1].translate(self.forwardPam[::-1].maketrans(watson, crick))

	@property
	def name(self):
		return self._name
	@name.setter
	def name(self, value):
		self._name = value

	@property
	def forwardPam(self):
		return self._forwardPam
	@forwardPam.setter
	def forwardPam(self, value):
		value = value.upper()
		if not all([c in list(IUPAC.keys()) for c in value]):
			raise ValueError(f"Reported PAM, {value}, contains non-IUPAC characters.")
		self._forwardPam = value

	@property
	def reversePam(self):
		return self.getReversePam()

	@property
	def primeness(self):
		return self._primeness
	@primeness.setter
	def primeness(self, value):
		if value not in ["3'","5'"]:
			raise ValueError("must specify 3' or 5'.")
		self._primeness = value



def get_cas_enzyme(name, cas_file=CAS_PATH):
	"""
	Function to return a specific Cas protein from a master list, in CAS_PATH
	by default. Since CAS_PATH is a text file, it is easily updated when new Cas
	enzymes are discovered.
	"""
	for line in open(cas_file):
		if not line.startswith("#"):
			cas = line.rstrip().split('\t')
			if cas[0] == name:
				return Cas(*cas)
	raise ValueError(f"Cas not found in {cas_file}: {name}") 

def get_cas_list(cas_file=CAS_PATH):
	"""
	Return list of all Cas9 names.
	"""
	cas_list = []
	for line in open(cas_file):
		if not line.startswith("#"):
			cas_list.append(line.split('\t')[0])
	return cas_list 

## Below: for testing, will be removed in next pull or by 4/15/2018

# test1 = get_cas_enzyme('SpCas9')
# print(test1)
# print(test1.forwardPam)
# print(test1.forwardPam_regex())
# print(test1.reversePam)
# print(test1.reversePam_regex())
# print(test1.primeness)

# print(len(test1.forwardPam))

# print(get_cas_list())