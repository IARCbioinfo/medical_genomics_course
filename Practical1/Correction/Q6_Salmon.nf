nextflow.enable.dsl=2

params.ref = "ref.fa"
params.input = "fastq/"

refch = file(params.ref)
fastqch = channel.fromFilePairs("${params.input}/*_{1,2}.fastq.gz")

process buildIndex(){
	input:
		path ref
	output:
		path "transcripts_index"
	publishDir "index", mode: "copy"
	shell:
	'''
	salmon index -t !{ref} -i transcripts_index
	'''
}

process quant(){
	input:
		path index
		tuple val(ID), path(fastq)
	output:
		path "${ID}"
	publishDir "expression", mode: "copy"
	shell:
        '''
        salmon quant -i !{index} -l A -1 !{fastq[0]} -2 !{fastq[1]} -o !{ID}
        '''
}

workflow {
   buildIndex(refch)
   quant( buildIndex.out, fastqch)
}

