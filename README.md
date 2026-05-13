# SeededAlignment.jl (now with functional_filter_and_align)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/SeededAlignment.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml?query=branch%3Amain)

**`SeededAlignment.jl`** aims to provide a frameshift robust framework for codon alignments. It can support
downstream analyses such as selection detection, phylogenetic inference, and evolutionary modeling. The package
is design to be flexible and performant enough for large-scale codon alignments.

## Installation

You can install this package using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/MurrellGroup/SeededAlignment.jl")
```

---

## Core Methods

-**`seed_chain_align`** - frameshift robust seed-and-chain strategy for codon alignment. Significantly faster than `nw_align` on biologically meaningful data. 

-**`msa_codon_align`** - frameshift robust codon-level MSA (visual) constructed by scaffolding pairwise codon alignments against a trusted reference.

-**`clean_frameshifts`** - cleans pairwise (or multiple) sequence alignments of frameshift errors if they contain a trusted reference.

-**`nw_align`** - frameshift robust codon-level Needleman-Wunsch algorithm that detects and cleans frameshift errors.

-**`filter_and_align`** - filter sequences with frameshift errors and codon align functional sequences. 

**Note:** **`nw_align`** and **`seed_chain_align`** also support De-Novo nucleotide alignments. 

**FASTA Support** – **`read_fasta`** / **`write_fasta`**

--- 

## Examples

### seed_chain_align

---

This is the main method for computing pairwise codon alignments. 

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("example.fasta")
# this sequence has no frameshift errors - reference sequence
ref_CDS = dna_seqs[1]
# this sequence might have frameshift errors - query sequence
query_CDS = dna_seqs[2]
# compute codon alignment by cleaning frameshift errors
codon_alignment = seed_chain_align(ref=ref_seq, query=query_seq)
```
If we want to view what was cleaned up in prinouts we can run the last line again with the verbose kwarg
```julia
codon_alignment = seed_chain_align(ref=ref_seq, query=query_seq, verbose=true)
```

---

### clean_frameshifts

---

In the situation where we have a pairwise alignment that has a broken reading frame due to frameshift errors. We can clean up 
the frameshift errors if one of the sequences is a trusted reference sequence. 

```julia
using SeededAlignment

# read in raw_alignment that contains frameshift errors
seq_names, raw_alignment = read_fasta("alignment_frameshift_error.fasta")
# aligned reference sequence that is frameshift free
aligned_ref = raw_alignment[1]
# aligned query sequence that contains frameshift errors
aligned_query = raw_alignment[2]
# cleans up frameshift errors
cleaned_alignment = clean_frameshifts(aligned_ref, aligned_query)
# write the clean result to fasta file
write_fasta("cleaned_codon_alignment.fasta", cleaned_alignment, seq_names=seq_names)
```
By default the changes are visibile in printouts just like `seed_chain_align` with verbose kwarg.

---

### msa_codon_align

---

This method produces a reference-guided codon multiple sequence alignment that is frameshift robust. This might be helpful 
because traditionally if a single sequence in the alignment contains a frameshift error it breaks the entire codon alignment and 
one would have to currate a correction manually.

```julia
using SeededAlignment

# read in dna_sequences
seq_names, dna_seqs = read_fasta("example.fasta")
# choose a reference sequence with intact reading frame
ref_seq = dna_seqs[1]
# extract the query sequences
query_seqs = dna_seqs[2:end]
# produce a codon alignment based on the reference sequence reading frame
msa_codon_alignment = msa_codon_align(ref_seq, query_seqs)
# write alignment to fasta file
write_fasta("example_output.fasta", msa_codon_alignment, seq_names = seq_names)
```
Note that the resulting alignment will be biased by the reference sequence and might need to be supplemented if the
intended alignment is complex.  

---

### filter_and_align

---

This method produces a reference-guided codon multiple sequence alignment in fasta format for those sequences whose 
longest open reading frame is free of frameshift errors and matches the reference to a user specified theshold.  
Sequences failing the filter are annotated and written to a rejection fasta file.
```julia
using SeededAlignment

# extract the longest open reading frame from each query sequence and codon align to a reference.
# keep those ORFs that have no reading frame errors and that match the reference to a user supplied threshold
# otherwise reject and annotate the sequence name with the rejection reason
# On completion, perform a codon aware alignment of passing queries and write aligned passing functionals and
# unaligned rejects to fasta files. Also populate a dataframe with housekeeping entries.
match=0.7
hk = filter_and_align(ref_filename, query_filename, functionals_filename, rejects_filename, match_thresh=match)
@show hk
```


