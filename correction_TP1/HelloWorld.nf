params.greeting = "Hello World!"
greetch = Channel.from(params.greeting)

process sayhello {
	input:
	val greet from greetch

	output:
	stdout into out

	shell:
	'''
	echo "!{greet}"
	'''
}

out.view()
