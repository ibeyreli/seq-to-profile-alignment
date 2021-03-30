# Sequence To Profile Alignment
Aim: Given one alignment and one new sequence, find the sequence to profile
alignment by utilizing naive Needleman-Wunsch method where gap penalty, 
match score, and mismatch penalty will be given as parameters.

# alignSeqToProfile.py

alignSeqToProfile is a program that performs sequence to profile alignment.

## Installation

alignSeqToProfile uses Python 3 standard libraries.
If you have installed Python 3 before, no additional package is needed.
If not, you may install Python 3 from "https://www.python.org/downloads/".
Make sure you have selected the appropriate distribution for your operating system.

## Usage

You can run the program using terminal.

The syntax:

"python alignSeqToProfile.py --fasta seq.fa --aln aligned_sequences.aln --out seq.aln
 --gap ${gap_penalty} --match ${match_score} --mismatch ${mismatch_penalty}"
 
--fasta					: Indicator that the next input is the input file which
						contains the sequence to be aligned
seq.fa					: Sequence fasta file to be aligned to the given profile
--aln 					: Indicator that the next input is the input file which 
						contains the given alignments
aligned_sequences.aln	: Given alignments which contains n DNA sequences line-by-line 
--gap 					: Indicator that the next input is gap penalty score
--match 				: Indicator that the next input is matching score
--mismatch				: Indicator that the next input is mismatch penalty score
--out 					: Indicator that the next input is the name of the output 
						file which will contain the final alignment
seq.aln 				: Output alignment file

## Examples:

aligned_sequences.aln :
sequence1 ATAC---CTAATTGGAGATGATCAAATTTATAAT
sequence2 TTAT---CTAATTGGCGACGATCAAATTTATAAT
sequence3 ATAT---CTAATTGGTGATGATCAAATTTATAAT
sequence4 ATCA---TTAATTGGAGATGATCAATCCTAATGA
sequence5 CTGTACTTTTATTGGTGATAGTCAAATCTATAAT

seq.fa :
> sequence
CTAGATAATTGGAGATGATCAAATTTATAT

Input >>> python alignSeqToProfile.py --fasta seq.fa --aln aligned_sequences.aln --out seq.aln  --gap -2 --match 1 --mismatch -1

Output >>>
out.aln :
sequence1 ATAC---CTAATTGG-AGATGATCAAATTTATAAT-----
sequence2 TTAT---CTAATTGG-CGACGATCAAATTTATAAT-----
sequence3 ATAT---CTAATTGG-TGATGATCAAATTTATAAT-----
sequence4 ATCA---TTAATTGG-AGATGATCAATCCTAATGA-----
sequence5 CTGTACTTTTATTGG-TGATAGTCAAATCTATAAT-----
sequence  C-T---AGA-TAATTGGAGA-TGA-TCAAAT--T-TATAT
