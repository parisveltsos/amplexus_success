# Scripts used for Meitreile field estimates of sexual selection on different Y haplogypes

## Preliminaries

The scripts should work if you put the analysis files in `~/git/amplexus_success/`. If not you will need to change the working directory in the scripts.

The output folder is not synced and you will need to create it

	mkdir `~/git/amplexus_success/output/'

Prerequisite packages (they appear in the beginning of R scripts) can be installed with

	install.packages("package_name", dependencies=TRUE)

When all scripts are run, the output folder is <1 Mb.

## Amplexus success

Run the script `scripts/amplexus.r`. Most of the analysis results are produced by this script.

The data headers are

id - individual ID
sex - individual sex
amplexus_dead - whether found in amplexus with a dead female
amplexus - whether found in amplexus
l1 - snout to vent length (cm)
l2 - back-leg length (cm)
w - body weight (g)
date - date of measurement
year - year of measurement
time - time of measurement
day - day batch of measurement
field - measurement in field or lab (weight in lab is discarded)
pop - population of individual (all are Meitreile)
day_numeric - day of capture since beginning of season
genotype - genotype characterisation. This is converted to haplotype-level of differentiation or Dmrt type, within the script.

## Fathering success

Run the script `scripts/fathering.r`. 

The data are the same as the amplexus success data in a "pre" fathering status, plus information on the genotype of the father from genotyping clutces ("mated" fathering status).

## PCA

Run the script `scripts/pca.r`.

The output illustrates the level of differentiation of the Y corresponds to the naming used.
