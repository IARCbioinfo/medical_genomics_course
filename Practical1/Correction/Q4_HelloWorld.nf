//a) 
params.greeting = "Hello World!"
greetch = channel.from(params.greeting)

//b)
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

//c)
workflow {
   sayhello(greetch).view()
}
