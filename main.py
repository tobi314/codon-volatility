#!/usr/bin/python3
# -*- coding: utf-8 -*-

WRONG_SEQ_ERROR_STR = "Sequence must consist only of Letters A, C, G and U and be divisible by 3!"

CODONS = {('A', 'A', 'A'): "K", ('A', 'A', 'C'): "N", ('A', 'A', 'G'): "K", ('A', 'A', 'U'): "N", 
		  ('A', 'C', 'A'): "T", ('A', 'C', 'C'): "T", ('A', 'C', 'G'): "T", ('A', 'C', 'U'): "T", 
		  ('A', 'G', 'A'): "R", ('A', 'G', 'C'): "S", ('A', 'G', 'G'): "R", ('A', 'G', 'U'): "S", 
		  ('A', 'U', 'A'): "I", ('A', 'U', 'C'): "I", ('A', 'U', 'G'): "M", ('A', 'U', 'U'): "I", 
		  ('C', 'A', 'A'): "Q", ('C', 'A', 'C'): "H", ('C', 'A', 'G'): "Q", ('C', 'A', 'U'): "H", 
		  ('C', 'C', 'A'): "P", ('C', 'C', 'C'): "P", ('C', 'C', 'G'): "P", ('C', 'C', 'U'): "P", 
		  ('C', 'G', 'A'): "R", ('C', 'G', 'C'): "R", ('C', 'G', 'G'): "R", ('C', 'G', 'U'): "R", 
		  ('C', 'U', 'A'): "L", ('C', 'U', 'C'): "L", ('C', 'U', 'G'): "L", ('C', 'U', 'U'): "L", 
		  ('G', 'A', 'A'): "E", ('G', 'A', 'C'): "D", ('G', 'A', 'G'): "E", ('G', 'A', 'U'): "D", 
		  ('G', 'C', 'A'): "A", ('G', 'C', 'C'): "A", ('G', 'C', 'G'): "A", ('G', 'C', 'U'): "A", 
		  ('G', 'G', 'A'): "G", ('G', 'G', 'C'): "G", ('G', 'G', 'G'): "G", ('G', 'G', 'U'): "G", 
		  ('G', 'U', 'A'): "V", ('G', 'U', 'C'): "V", ('G', 'U', 'G'): "V", ('G', 'U', 'U'): "V", 
		  ('U', 'A', 'A'): "Z", ('U', 'A', 'C'): "Y", ('U', 'A', 'G'): "Z", ('U', 'A', 'U'): "Y", 
		  ('U', 'C', 'A'): "S", ('U', 'C', 'C'): "S", ('U', 'C', 'G'): "S", ('U', 'C', 'U'): "S", 
		  ('U', 'G', 'A'): "Z", ('U', 'G', 'C'): "C", ('U', 'G', 'G'): "W", ('U', 'G', 'U'): "C", 
		  ('U', 'U', 'A'): "L", ('U', 'U', 'C'): "F", ('U', 'U', 'G'): "L", ('U', 'U', 'U'): "F"}

def hamming_d(c1: tuple, c2: tuple) -> int:
	d = 0
	for i in range(len(c1)):
		if c1[i] != c2[i]:
			d += 1
	return d

class Codon:
	def __init__(self, bases: str):
		self.bases = tuple(bases)
		self.amino_letter = CODONS[self.bases]

	def __repr__(self):
		return "".join(self.bases)+f" ({self.amino_letter})"

	def get_volatility(self) -> int:
		neighbors = 0; different = 0
		for bases, amino_letter in CODONS.items():
			if amino_letter == "Z" or hamming_d(self.bases, bases) != 1:
				continue
			neighbors += 1
			if amino_letter != self.amino_letter:
				different += 1
		return different/neighbors

class Sequence:
	def __init__(self, seq: str):
		if not (all(x in ("A", "C", "G", "U") for x in seq) and len(seq) % 3 == 0):
			raise TypeError(WRONG_SEQ_ERROR_STR)
		else:
			self.codons = []
			for chars in [seq[i:i+3] for i in range(0, len(seq), 3)]:
				self.codons.append(Codon(chars))

	def plot_volatility(self, chunk_len: int) -> list:
		data =  []
		for chunk in [self.codons[i:i+chunk_len] for i in range(0, len(self.codons), chunk_len)]:
			volatility = []
			for codon in chunk:
				volatility.append(codon.get_volatility())
			data.append(sum(volatility)/len(volatility))
		print(data)

if __name__ == '__main__':
	s = Sequence("CGGGCUUGGAGCCCUUGCCGU")
	s.plot_volatility(2)