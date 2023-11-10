#a) from the help, the pipeline can be run with: nextflow run IARCbioinfo/fastqc-nf -r v1.1 -profile singularity --input_folder input --output_folder results . We just need to adapt the "input" to match the proper data test folder with fastq files.
#b)
nextflow run IARCbioinfo/fastqc-nf -r v1.1 -profile singularity --input_folder data_test/FASTQ --output_folder results
