//a)
params.ref = "ref.fa"
refch = file(params.ref)

//b)
process buildIndex(){
	input:
		path ref
	output:
		path "transcripts_index"
	//c)
	publishDir "index", mode: "copy"

	//d)
	shell:
	'''
	salmon index -t !{ref} -i transcripts_index
	'''
}

workflow {
   buildIndex(refch)
}

