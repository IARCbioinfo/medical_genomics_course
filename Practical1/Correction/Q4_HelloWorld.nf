nextflow.enable.dsl=2

params.greeting = "Hello World!"
greetch = channel.from(params.greeting)

process sayhello(){
	input:
		val x
	output:
		stdout
	shell:
	'''
	echo !{x}
	'''
}

workflow {
   sayhello(greetch).view()
}
