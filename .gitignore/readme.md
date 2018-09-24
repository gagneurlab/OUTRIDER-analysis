# OUTRIDER analysis

This is the accompanying analysis repository of the paper: ` OUTRIDER: A statistical method for detecting aberrantly expressed genes in RNA sequencing data` 

The paper can be found [here](https://www.biorxiv.org/content/early/2018/06/14/322149).
A full copy of the resulting files can be found on our webserver [https://i12g-gagneurweb.in.tum.de/public/paper/OUTRIDER/](https://i12g-gagneurweb.in.tum.de/public/paper/OUTRIDER/).

This repository contains the full pipeline and code to reproduce the results published in the paper. 

## Project structure

This project is setup as a [wBuild](https://github.com/gagneurlab/wBuild). This is an automatic build tool for R reports based on [snakemake](https://snakemake.readthedocs.io/en/stable/).

* The `Scripts` folder contains scripts which will be rendered as HTML reports
* The `src` folder contains additional helper functions and scripts
* The `Output` folder will contain all files produced in the analysis pipeline
    * `Output/data` has all raw RDS output files
    * `Output/html` contains the final HTML report
    * `Output/paper_figures` has all paper figures

## Dependencies

This project depends on the python package `wBuild` and the R package `OUTRIDER`.
The data will be downloaded automatically. Since the genotypes are not publicly shareable one has to apply for the data at [dbGaP]( https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v6.p1). We included a fake VCF to run the full pipeline. Please replace the path to the VCF files in the [wbuild.yaml](./wbuild.yaml) file to include the real enrichment data.

## Repository setup

First download the repo and its dependencies:

```
git clone https://github.com/gagneurlab/OUTRIDER OUTRIDER
git clone https://github.com/gagneurlab/OUTRIDER-analysis OUTRIDER-analysis
cd OUTRIDER-analysis
```

and install wbuild using pip by running.

```
pip install wBuild
wBuild init
```

Since `wBuild init` will reset the current `Snakefile`, ` readme.md`, and `wbuild.yaml` we have to revert them again with git.

```
git checkout Snakefile
git checkout wbuild.yaml
git checkout readme.md
```

To make sure all packages needed in the analysis are installed source the following file in R

```
Rscript ./src/r/install_dependencies.R
```

## Run the full pipeline

To run the full pipeline, execute the following command with 10 cores in parallel:

```
snakemake -j 10
```

or to run it on the cluster with SLUM installed: 

```
snakemake -k --restart-times 4 --cluster "sbatch -N 1 -n 10 --mem 30G" --jobs 16
```
