

This code is primarily designed to process BAM files (along with VCF files) to generate haplotype fragment files that can then be phased using HapCUT. Currently, it has the following capabilities:

1. identify fragments from dilution pool sequence data (see Kitzman et al. 2011 and Kaper et al. PNAS 2012) using a DP based segmentation algorithm that accounts for mappability and fragment lengths as well as overlapping fragments. BAM file for each pool is processed separately and outputs a HapCUT compatible fragment file

2. identify fragments using an input bed file that specifies the location of the fragment intervals. Such a bed file can be generated using samtools targetcut 

3. Process 10X barcoded BAM file (each read has a barcode) to generate haplotype fragments. This requires an input bed file with the intervals and corresponding barcodes. 




HISTORY:

code first started in late 2013, working version developed in Dec 2014.  

How to phase fosmid pooled sequencing data (Kitzman et al. Nat. Biotech. 2011)

    Align reads for each pool to reference genome using your favorite aligner.
    Use 'samtools targetcut' command to generate consensus sequences for each long fragment/fosmid
    run 'extracthairs' to generate variant calls at heterozygous sites in each fragment (this should be done for each pool separately)
    combine the fragment files from different pools into a single fragment file
    run HAPCUT using this fragment file and VCF variant file using options "--fosmids 1" and "--sf 1"


