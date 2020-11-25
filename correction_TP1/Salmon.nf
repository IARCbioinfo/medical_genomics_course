fastq = Channel.fromFilePairs(params.input_folder+"/*_{1,2}.fastq.gz")
		.view()

ref = file(params.ref)

process index {
	input:
	path ref

	output:
	path("INDEX") into index
	publishDir "results/", mode: "copy"

	shell:
	'''
	salmon index -t !{ref} -i INDEX
	'''
}


process quantif {
	cpus 2
 
	input:
	path index
	tuple val(name), path(reads) from fastq

	output:
	path("quantif_${name}") into quantif_files
	publishDir "results", mode: 'copy'

        shell:
        '''
        salmon quant -i !{index} -l A -1 !{reads[0]} -2 !{reads[1]} -o quantif_!{name} -p 2
	'''
}
