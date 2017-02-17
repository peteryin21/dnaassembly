# Python De Novo DNA Assembly

## Usage
Get the requirements
```
pip install -r requirements.txt
```
Assemble away!
```
python assembly.py fasta.txt
```

Run unit tests
```
python test.py
```

This DNA sequence assembler works under the assumption that there are no read errors, and that there is overlap for each of the sequences in the file. 

I took an iterative approach. First, I take one sequence and find the sequence in the list that has the best overlap score. Merge the sequences and repeat with the remainder of the list until all sequences have been combined.

If n is the average length of sequence, and m is the number of sequences, this has O(mn^2) complexity. Finding best overlap of two strings is O(n^2) in order to explore all the offsets and check overlap matching. We have to repeat this m times for each sequence. In this example with 50 sequences of ~1000 length, this naive implementation suffices, however if we want to scale to many sequences of longer length, a more efficient algorithm (i.e De Bruijn graph) should be considered.

TODO:

- Improve efficiency 
- Handle sequences with no overlap and imperfect overlap