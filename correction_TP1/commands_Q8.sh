wget http://refgenomes.databio.org/v2/asset/hg38/salmon_partial_sa_index/archive?tag=default salmon_sa_default.tgz
tar -xzvf salmon_sa_default.tgz

nextflow run SalmonSRA.nf -c nextflow.config -profile conda \
    --ref /home/nalcala/Medical_Genomics_TP1/data_test/REF/17.fasta

