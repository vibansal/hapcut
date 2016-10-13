
This directory contains code for to extract fragment sequences from the bam files for each pool containing sequence data from long DNA fragments. The code was first started in early 2014 and has been updated sporadically. The current code is able to reconstruct the locations of the fragments but needs some work to handle overlapping fragments. 


HISTORY:

code first started in late 2013, working version developed in Dec 2014 


How to phase fosmid pooled sequencing data (Kitzman et al. Nat. Biotech. 2011)

    Align reads for each pool to reference genome using your favorite aligner.
    Use 'samtools targetcut' command to generate consensus sequences for each long fragment/fosmid
    run 'extracthairs' to generate variant calls at heterozygous sites in each fragment (this should be done for each pool separately)
    combine the fragment files from different pools into a single fragment file
    run HAPCUT using this fragment file and VCF variant file using options "--fosmids 1" and "--sf 1"


