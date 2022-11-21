params.index = null
params.ref = "ref.fa"
params.input = "fastq/"
params.cpus = 2
params.memory = 5
//d)
params.output = "results"

refch = file(params.ref)
fastqch = channel.fromFilePairs("${params.input}/*_{1,2}.fastq.gz")

process buildIndex(){
	//b)
	cpus params.cpus
        memory params.memory+"GB"

	input:
		path ref
	output:
		path "transcripts_index"
	publishDir "$params.output/index", mode: "copy"
	shell:
	'''
	salmon index -t !{ref} -i transcripts_index
	'''
}

process quant(){
	//b)
	cpus params.cpus
	memory params.memory+"GB"

	input:
		path index
		tuple val(ID), path(fastq)
	output:
		path "${ID}"
	publishDir "$params.output/expression", mode: "copy"
	shell:
        '''
        salmon quant -i !{index} -l A -1 !{fastq[0]} -2 !{fastq[1]} -o !{ID}
        '''
}

workflow {
   //a)
   if(!params.index){
     ind = buildIndex(refch)
   }else{
     ind = channel.from(params.ref)
   }
   quant( ind, fastqch)
}

