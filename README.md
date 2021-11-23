# xPORE_nextflow

Install
1. Nanopolish
	a. git clone --recursive https://github.com/jts/nanopolish.git
	b. cd nanopolish/
	c. make
2. Miniconda
3. Install environment with yaml file
	a. conda env create --file xPORE.yml

USAGE
1. conda activate xPORE
2. nextflow run process_xpore.nf -resume -with-report logs/report.html
