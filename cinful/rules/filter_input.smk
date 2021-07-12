
rule filter_microcin:
	input:
		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/microcins/{sample}.filtered.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"

rule filter_immunity_protein:
	input:
		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/immunity_proteins/{sample}.filtered.fa"
	shell:
		"seqkit seq -m 30 -M 250  {input} | seqkit rmdup -s > {output}"

rule filter_CvaB:
	input:
		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/{sample}.filtered.fa"
	shell:
		"seqkit rmdup -s {input} > {output}"
