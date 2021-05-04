SAMPLES, = glob_wildcards("{sample}.fna")

rule final:
	input:
		expand("{sample}_cinfulOut/{sample}.faa", sample = SAMPLES)

rule prodigal:
	input:
		"{sample}.fna"
	output:
		gff3="{sample}_cinfulOut/{sample}.gff3",
		aa="{sample}_cinfulOut/{sample}.faa"
	shell:
		"prodigal -i {input} -o {output.gff3} -a {output.aa}"