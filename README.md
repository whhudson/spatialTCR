# spatialTCR

The folder spatial_images contains microscopy images of H&E stained human brain metastases used in spatial transcriptomics experiments. Files are in .jpg format and are numbered by patient sample.

For TCR analysis with the Visium platform, downsampled data (100,000 of ~2.1 million reads) of TCRβ library sequencing from patient 16 are located in the folder downsampled_reads. Read 2 of these sequences can be given to MiXCR for antigen receptor calling with the sample code below:

mixcr analyze shotgun -s hsa --starting-material rna --align "-OsaveOriginalReads=true" --assemble "--write-alignments" path/to/pt16_downsample_L001_R2_001.fastq.gz mixcr_output/clones


Reads can be exported with:

mixcr exportReadsForClones -s mixcr_output/clones.clns mixcr_output/reads/reads.fastq.gz 

This will result in the MiXCR output found in the mixcr_ouput folder; MiXCR version 3.0.13 was used for this analysis. The mixcr_output/reads folder will contain reads mapping to each clone in the format reads_clnN.fastq.gz, where N is the clone number. Clone information will be given in the mixcr_output/clones.clonotypes.ALL.txt file. The reads given in the mixcr_output/reads folder can be mapped to a spatial location by using the spatial barcode and UMI sequences found in the paired read 1 file (downsampled_reads /pt16_downsample_L001_R1_001.fastq.gz here). 

An example Python script implementation is given in the code_samples folder. combine_counts.py takes the mixcr reads folder (mixcr_output/reads here), the read 1 sequencing files (downsampled_reads /pt16_downsample_L001_R1_001.fastq.gz), and a spatial barcode file from 10X Genomics’s cellranger output. cellranger output for patient 16 is given here in the spaceranger_output folder, spaceranger_output/spatial/tissue_positions_list.csv is the file needed by combine_counts.py. This python script will output a matrix where spatial position barcodes are given in rows and UMI counts are given in columns. This output file is given in code_samples/output_counts.txt

Finally, an example visualization of a clone with the Seurat package in R is shown in the script plot_clone.R. Note that this shows only clones found in the downsampled reads and is not indicative of all clones found in the complete sequencing dataset.
