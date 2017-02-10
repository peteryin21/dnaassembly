#!/bin/sh 
from Bio import SeqIO
import sys 
import argparse


def find_best_overlap(s1, s2):
    """Take in 2 strings, and return score and offset for best overlap"""
    min_offset = 1-len(s2)
    max_offset = len(s1)
    offset_scores = []
    for offset in range(min_offset, max_offset):
        # For each offset, loop through overlap region
        score = 0
        for pos in range(max(offset, 0), min(len(s1),len(s2)+offset)):
            if s1[pos] == s2[pos-offset]:
                score += 1
            else: 
                offset_scores.append((-1, offset))
                break
        offset_scores.append((score, offset))
    return sorted(offset_scores, reverse=True)[0]


def best_match(seq, list_seq):
    """Take in a sequence and a list of sequences, return the best matching sequence and the offset"""
    results=[]
    for s in list_seq:
        if seq != s:
            score = find_best_overlap(seq, s)
            results.append((score, s))
    return sorted(results, key=lambda x: x[0], reverse=True)[0]


def merge_seq(s1, s2, offset):
    """Return the merge of 2 sequences based on overlap offset"""
    if offset > 0:
        return s1[:offset] + s2
    else:
        return s2[:-offset] + s1


def assemble(seqs):
    """Take in a list of strings and return fully assembled sequence"""
    seq_list = sorted(seqs)
    s1 = seq_list[0]
    seq_list.remove(s1)
    while seq_list:
        best = best_match(s1, seq_list)
        offset = best[0][1]
        s2 = best[1]
        new_seq = merge_seq(s1, s2, offset)
        seq_list.remove(s2)
        s1 = new_seq
    return new_seq
        

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Assemble sequences from a FASTA file, writes to result.txt")
	parser.add_argument('fasta', help='Provide the FASTA file')
	args = parser.parse_args()
	fasta_seq = SeqIO.parse(args.fasta, 'fasta')
	seq_list = [str(fasta.seq) for fasta in fasta_seq]
	with open("result.txt", 'w') as out:
		out.write(assemble(seq_list))
	print("Assembled sequence written to result.txt")

