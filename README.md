ModTect
==========

Overview
-------------
ModTect is a computational tool for detecting RNA modifications which disrupts base-pairing, from typical RNA-sequencing datasets. The method relies upon two features induced by base-pair-disrupting RNA modifications, (1) the multi-nucleotide-mismatch signal, and (2) the deletion signal. ModTect automatically extracts both signals from typical RNA-sequencing libraries, using a specialized statistical model, to determine possible RNA modification sites.

ModTect is written in Python and should not take >1-2 minutes to set-up on a system pre-installed with both Python2.7 and SAMtools.


System requirements and dependencies
-------------

ModTect was developed in a Linux environment. The following packages are also required for running:
- Python 2.7 (ModTect can also be executed using pypy for a ~2x speedup)
- SAMtools (version 0.1.19, or 1.3.1 and later)


Installation
-------------

1) Download or clone the repository.

```bash
git clone https://github.com/ktan8/ModTect.git
```

2) Ensure that both Python 2.7 and SAMtools are installed and accessible in your $PATH. Note that python and samtools has to be callable via the following command.

Check if python 2.7 is accessible directly in your path:
```bash
python -V
```

Check if samtools is accessible directly in your path:
```bash
samtools
```

NOTE: It is critical for samtools to be directly callable with the "samtools" command, as ModTect makes a direct call to samtools in your environment

Running ModTect
-------------
The easiest way to run ModTect is as follows

```bash
python ModTect.py <bamfile> <Reference Genome> <chromosome> <position_start> <position_end>
```

A full example is:
```bash
python ModTect_1_7_5.py ./sample/CML.MALAT1.bam Homo_sapiens_assembly19.fasta 11 65273620 65273640
```
NOTE: The hg19 reference genome used for the sample dataset is the hg19 Broad variant (ftp://ftp.broadinstitute.org/pub/svtoolkit/reference_metadata_bundles/Homo_sapiens_assembly19_25Jan2015.tar.gz). In this reference genome, chromosomes are labelled "1,2,3,4,..." instead of "chr1,chr2,chr3,chr4,..."

A few optional arguments are also provided if you like to make modifications to it.

```
[--basesToTrimFromEdge BASESTOTRIMFROMEDGE]
[--scoreCutoff SCORECUTOFF] 
[--minDepth MINDEPTH]
[--readlength READLENGTH]
[--mismatchPropFilter MISMATCHPROPFILTER]
[--deletionCountCutoff DELETIONCOUNTCUTOFF]
[--threads THREADS]
[--regionFile REGIONFILE]
[--label LABEL]
```


Sample Test Run
-------------

We have included a small sample dataset which you can use to test out ModTect. To run it, use the following command:
```bash
python ModTect_1_7_5.py ./sample/CML.MALAT1.bam Homo_sapiens_assembly19.fasta 11 65273620 65273640
```
Please remember to point the "Homo_sapiens_assembly19.fasta" reference genome file to the path in your local machine.

You should see the following output:
```
chrom	position	reference_nt	mono-alleleLogP	bi-alleleLogP	tri-alleleLogP	tetra-alleleLogP	variant_proportion	depth	deletionCount	deletion_proportion	refUpperCount	refLowerCount	altUpperCount	altLowerCount	refPlusProp	altPlusProp	A_count	T_count	G_count	C_count	a_count	t_count	g_count	c_count	has_reference_nt	types_of_nt	ModTect_score	variant_readposn_median_fwd	variant_readposn_median_rev	median_abs_dev_fwd	median_abs_dev_rev
11	65273630	A	-183.627994199	-82.6673603595	-29.4593841182	-29.4593841182	0.734375	69	0.0724637681159	9	8	26	21	0.529411764706	0.553191489362	9	15	11	0	8	15	6	0	1	3	53.2079762413	44.5	38	21.5	22.0
```

We see the positive sites identified from this region, which has a modification score of "53.2079762413". This indicates that it is 10 ^ 53.2 more likely to be a tri-nucleotide site than a bi-nucleotide site.


Description of output
-------------
These are the columns that you might find most useful:
- Modification score
- Deletion rate
- Variant proportion
- Depth



References
-------------
If you use ModTect for your work, please cite our study:

Kar-Tong Tan, Ling-Wen Ding, Chan-Shuo Wu, Daniel G. Tenen, Henry Yang. Repurposing RNA-Sequencing for de novo Discovery of Base-pair-disrupting mRNA modifications with ModTect.



Contact
-------------
If you have any queries or spot any bugs, please contact me at Kar-Tong Tan (ktan@broadinstitute.org)

