#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cas_object.py hold Cas info as part of ExcisionFinder. Written in Python v 3.6.1.
Kathleen Keough et al 2018.

cas_object.py will do two things:
	 + hold the Cas object info (one object per Cas enzyme) and 
	 will handle input error checking.
	 + contain a function to create the objects from an external list
	 CAS_LIST.txt
"""
import os, sys

CAS_PATH = f"{os.path.dirname(__file__)}/CAS_LIST.txt"

IUPAC = {
    "Y": "[CT]",
    "R": "[AG]",
    "W": "[AT]",
    "S": "[GC]",
    "K": "[TG]",
    "M": "[CA]",
    "D": "[AGT]",
    "V": "[ACG]",
    "H": "[ACT]",
    "B": "[CGT]",
    "N": "[ACGT]",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
}


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
        return "".join([IUPAC[c] for c in self.forwardPam])

    def reversePam_regex(self):
        """
		Returns the regex of the reverse PAM. 
		"""
        pam = self.getReversePam()
        return "".join([IUPAC[c] for c in pam])

    def getReversePam(self):
        """
		Returns the reverse PAM of the input.
		"""
        watson = "ACGTYRSWKMBDHVN"
        crick = "TGCARYSWMKVHDBN"
        return self.forwardPam[::-1].translate(
            self.forwardPam[::-1].maketrans(watson, crick)
        )

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
        if value not in ["3'", "5'"]:
            raise ValueError("must specify 3' or 5'.")
        self._primeness = value


def print_cas_types(cas_file=CAS_PATH):
    """
	Prints out all Cas9 names from get_cas_list() in a formated way.
	"""
    print("Available cas enzymes: ")
    print("\t", end="")
    for c in get_cas_list(cas_file):
        print(c, end=",")
    print()


def get_cas_enzyme(name, cas_file=CAS_PATH):
    """
	Function to return a specific Cas protein from a master list, in CAS_PATH
	by default. Since CAS_PATH is a text file, it is easily updated when new Cas
	enzymes are discovered.
	"""
    for line in open(cas_file):
        if not line.startswith("#"):
            cas = line.rstrip().split("\t")
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
            cas_list.append(line.split("\t")[0])
    return cas_list


def validate_cas_list(in_cas_list):
    """
	Takes in a list of cas enzymes, and validates if they are in CAS_LIST.txt. Returns two list, the first
	being the list of cas enzymes in both in_cas_list and CAS_LIST.txt, and the other being the list of
	enzymes not in CAS_LIST.txt. The second list is usefull for printing error messages in the main program.
	"""
    master = get_cas_list()
    in_both = list(set(master).intersection(in_cas_list))

    return in_both, list(set(in_cas_list).difference(in_both))
