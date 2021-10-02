import os 
import subprocess as sp
import argparse


class readable_dir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

class valid_file(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_file=values
        if not os.path.exists(prospective_file):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_file))












def main():
	# set up command line arguments
	parser = argparse.ArgumentParser(description='cinful')
	parser.add_argument('-d', '--directory', action=readable_dir,required=True)
	parser.add_argument('-o', '--outDir', type=str,default="cinful_out")
	parser.add_argument('-t', '--threads', type=int, default=1)
	args = parser.parse_args()

	print(args)
	threads = args.threads
	workdir = args.directory
	outdir = args.outDir

	# absolute path of this file, needed to find snakefile
	currentAbsPath = os.path.dirname(os.path.abspath(__file__))
	snakefile = f"{currentAbsPath}/Snakefile"
	print(os.path.dirname(os.path.abspath(__file__)))
	print(snakefile)

		# snakemake command to run
	cmd = [		"python",
				"-m",
				"snakemake",  
				"--snakefile",
				snakefile,
				"-j",
				str(threads),
				"--directory",
				workdir,
				"--config",
				f"outdir={outdir}"
			]
	print("Running the following command:")
	print(" ".join(cmd))
	sp.check_output(cmd)



if __name__ == "__main__":
    main()


