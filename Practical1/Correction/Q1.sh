#a)
nextflow run HelloWorld
#b) Nextflow starts running (version, etc), then fetches from github the pipeline HelloWorld [unique tag for run] mentions github revision number and [branch], mentions where the script is executed (e.g., local, amazon cloud, cluster), then executes the process and monitor the number of successful executions (%), and finally prints the results.
#c) 
nextflow run HelloWorld
nextflow run HelloWorld
# the order of the outputs change, because Nextflow submits them simultaneously and displays the results in the order in which they finish executing, which depends on many factors on each machine (which CPU is available, etc) and varies in each execution.
#d)
nextflow run HelloWorld -resume
# nextflow mentions "cached", which means that processes are not re-executed but the output is just read from the previous work directories.
