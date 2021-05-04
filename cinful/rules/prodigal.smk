SAMPLES, = glob_wildcards("{sample}.fna")

rule final:
	input:
		expand("{sample}.faa", sample = SAMPLES)

rule prodigal:
	input:
		"{sample}.fna"
	output:
		gff3="{sample}.gff3",
		aa="{sample}.faa"
	shell:
		"prodigal -i {input} -o {output.gff3} -a {output.aa}"