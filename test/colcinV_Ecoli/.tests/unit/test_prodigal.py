import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_prodigal():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/prodigal/data")
        expected_path = PurePosixPath(".tests/unit/prodigal/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
