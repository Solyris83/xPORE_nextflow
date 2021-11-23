# xPORE_nextflow

Install
1. Nanopolish \
		git clone --recursive https://github.com/jts/nanopolish.git\
		cd nanopolish/ \
		make \
2. Miniconda
3. Install environment with yaml file \
		conda env create --file xPORE.yml

USAGE
1. conda activate xPORE
2. nextflow run process_xpore.nf -resume -with-report logs/report.html
