#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DNA Analysis"""

__author__ = "shiva senthilkumar"


class Codon:
    """Amino-acid encoding with three nucleotides"""

    def __init__(self, bases):
        """accepting the correct DNA string; accepting both upper and lower
        case letters"""
        self.bases = bases
        if type(bases) is not str:
            raise ValueError("input is not a string")
        if len(bases) != 3:
            raise ValueError("input string does not contain exactly three " +
                             "characters")
        for element in range(0, len(bases)):
            index = bases.upper()[element]
            if (index == "A") or (index == "G") or (index == "T") or \
               (index == "C"):
                pass
            else:
                raise ValueError("input string contain invalid characters")

    def __str__(self):
        """Convert to string enclosed in square brackets.
        E.g. [GCT]
        """
        return ("[" + str(self.bases.upper()) + "]")

    def __eq__(self, other):
        """Compare if two Codon objects are equal (same sequence of bases)"""
        self.other = other
        if (self.bases.upper() == other):
            return True

    def transcribe(self):
        """Return a string of transcribed bases enclosed in angle brackets.

        E.g. <GCU>
        """
        bases = self.bases
        if "T" in bases.upper():
            transcribed_bases = bases.upper().replace("T", "U")
            return ("<" + transcribed_bases + ">")
        else:
            return ("<" + bases.upper() + ">")


class Gene:
    """Protein encoding with a sequence of codons"""

    def __init__(self, seq):

        self.seq = seq
        if type(seq) is not str:
            raise ValueError("input is not a string")
        for element in range(0, len(seq)):
            index = seq.upper()[element]
            if (index == "A") or (index == "G") or (index == "T") or \
               (index == "C"):
                pass
            else:
                raise ValueError("input string contain invalid characters")
        if len(seq) % 3 == 0:
            pass
        elif len(seq) % 3 == 1:
            seq = seq[0:-1]
        else:
            seq = seq[0:-2]
        self.seq_codons = [seq[x:x+3] for x in range(0, len(seq), 3)]

    def __str__(self):
        """Convert to string using a sequence of codons (in square brackets).

        E.g. [GCT][GGC]...
        """
        new_string = ""
        new_list = self.seq_codons
        i = 0
        while i < len(new_list):
            new_string = new_string + "[" + new_list[i].upper() + "]"
            i = i + 1
        return(str(new_string))

    def transcribe(self):
        """Return a string of transcribed codons (in anle brackets).

        E.g. <GCU><GGC>...
        """
        transcribed_string = ""
        transcribed_list = self.seq_codons
        i = 0
        while i < len(transcribed_list):
            transcribed_string = (transcribed_string + "<" +
                                  transcribed_list[i].upper() + ">")
            i = i + 1
        transcribed_string = transcribed_string.replace("T", "U")
        return(str(transcribed_string))

    def __contains__(self, codon):
        """Check if the gene sequence contains the given codon"""
        x = self.__str__()
        if str(codon).upper() in x:
            return True

    def gc_content(self):
        """Return the fraction of G and C bases relative to all bases"""
        seq = str(self.seq).upper()
        i = 0
        y = 0
        if len(seq) % 3 == 0:
            pass
        elif len(seq) % 3 == 1:
            seq = seq[0:-1]
        else:
            seq = seq[0:-2]
        while i < len(seq):
            if ('G' == seq[i]) or ('C' == seq[i]):
                y += 1
            i += 1
        z = float(y) / len(seq)
        return (z)


if __name__ == "__main__":

    # Read DNA sample file (ignore comments and newlines)
    dna_lines = []
    with open("dna_sample.txt") as datafile:
        for line in datafile:
            line = line.strip()
            if not line.startswith("#"):
                dna_lines.append(line)
    dna_sample = "".join(dna_lines)

    gene = Gene(dna_sample)

    # Feel free to change, delete or add to this testing code below
    # This part of your code and the printed outputs below will not be graded
    print(Codon("act") in gene)
    print(gene.gc_content())
