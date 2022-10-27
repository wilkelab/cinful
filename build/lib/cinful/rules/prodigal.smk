SAMPLES, = glob_wildcards("{sample}.fna")

rule prodigal:
	input:
		"{sample}.fna"
	output:
#		gff3=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.gff3",
#		cds=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.cds",
		aa=config["outdir"] + "/01_orf_homology/prodigal_out/{sample}.faa"
	priority: 50
	shell:
		"prodigal -p meta -i {input} -a {output.aa} -q"

		# """
		# let a=`wc -m < {input}`
		# let b=30000
		# if [ $a -gt $b ]
		# then
		# 	prodigal -p single -i {input} -a {output.aa} -q
		# else
		# 	prodigal -p meta -i {input} -a {output.aa} -q
		# fi
		# """ 
		
		# -o {output.gff3} -d {output.cds}
