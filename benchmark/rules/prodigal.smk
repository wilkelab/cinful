SAMPLES, = glob_wildcards("{sample}.fna")

# rule final:
# 	input:
# 		expand("cinfulOut/01_orf_homology/prodigal_out/{sample}.faa", sample = SAMPLES)

rule prodigal:
	input:
		"{sample}.fna"
	output:
		gff3="cinfulOut/01_orf_homology/prodigal_out/{sample}.gff3",
		cds="cinfulOut/01_orf_homology/prodigal_out/{sample}.cds",
		aa="cinfulOut/01_orf_homology/prodigal_out/{sample}.faa"
	benchmark:
		"cinfulOut/benchmark/prodigal/{sample}.txt"
	shell:
		"prodigal -p meta -i {input} -o {output.gff3} -a {output.aa} -d {output.cds}"