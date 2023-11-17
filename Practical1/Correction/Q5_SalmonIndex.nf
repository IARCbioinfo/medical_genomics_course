nextflow.enable.dsl=2

params.ref = "ref.fa"
refch = file(params.ref)

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

workflow {
   buildIndex(refch)
}

