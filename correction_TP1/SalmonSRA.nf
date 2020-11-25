fastq = Channel.fromSRA("SRP156394",max:1)
		.view()

index_ch = file(params.index_folder)

process quantif {
	cpus 2
 
	input:
	path index_ch
	tuple val(name), file(reads) from fastq

	output:
	path("quantif_${name}") into quantif_files
	publishDir "results", mode: 'copy'

        shell:
        '''
        salmon quant -i !{index_ch} -l A -1 !{reads[0]} -2 !{reads[1]} -o quantif_!{name} -p 2
	'''
}
