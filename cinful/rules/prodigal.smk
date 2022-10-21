SAMPLES, = glob_wildcards("{sample}.fna")

rule prodigal:
	input:
		"{sample}.fna"
	output:
		gff3=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.gff3",
		cds=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.cds",
		aa=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.faa"
	shell:
		"prodigal -p single -i {input} -o {output.gff3} -a {output.aa} -d {output.cds}"
