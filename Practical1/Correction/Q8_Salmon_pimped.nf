nextflow.enable.dsl=2

params.index = null
params.ref = "ref.fa"
params.input = "fastq/"
params.cpus = 2
params.memory = 5

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
	cpus params.cpus
	memory params.memory+"GB"

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
   if(!params.index){
     ind = buildIndex(refch)
   }else{
     ind = channel.from(params.ref)
   }
   quant( ind, fastqch)
}

