import unittest
from assembly import *


class TestAssembly(unittest.TestCase):

	def test_best_overlap_standard(self):
		self.assertEqual(find_best_overlap("ATTAGACCTG", "AGACCTGCCG"), (7,3))
		self.assertEqual(find_best_overlap("ATGC", "ATGC"), (4,0))

	def test_best_overlap_nooverlap(self):
		self.assertEqual(find_best_overlap("ATTT", "GCCCCCCC")[0], -1)
		self.assertEqual(find_best_overlap("GCCCCCCC", "ATTT")[0], -1)

	def test_best_match(self):
		self.assertEqual(best_match("ATGC", ["ATATG", "GC", "ATGC"]), ((4,0), 'ATGC'))
		self.assertEqual(best_match("AAAA", ["TTTT", "GGGG", "CCCC"]), None)

	def test_merge_seq(self):
		self.assertEqual(merge_seq("ATGCA","TGCATG", 1), "ATGCATG")
		self.assertEqual(merge_seq("TGCATG","ATGCA", -1), "ATGCATG")
		self.assertEqual(merge_seq("AAAA","AAAA", 0), "AAAA")

	def test_assemble_standard(self):
		self.assertEqual(assemble(["ATTAGACCTG", "CCTGCCGGAA", "AGACCTGCCG", "GCCGGAATAC"]), "ATTAGACCTGCCGGAATAC")
		self.assertEqual(assemble(["AAAAC", "CGGGG"]), "AAAACGGGG")

	def test_assemble_complete_overlap(self):
		self.assertEqual(assemble(["AA", "AAAA","AAAAAA"]), "AAAAAA")
		self.assertEqual(assemble(["ATGC", "ATGC"]), "ATGC")
		

if __name__ == '__main__':
	unittest.main()

