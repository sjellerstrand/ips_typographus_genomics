# Weak population genetic structure in Eurasian spruce bark beetle over large regional scales in Sweden
This repository contains lab protocols, RADtag library preperation summaries, and bioinformatic scripts used in:


Ellerstrand, S. J., Choudhury, C., Svensson, K., Andersson, M. N., Kirkeby, C., Powell, D., Schlyter, F., Jönsson, A. M., Brydegaard, M., Hansson, B. and Runemark, A. (Accepted May 2022) Weak population genetic structure in Eurasian spruce bark beetle over large regional scales in Sweden. Ecology and Evolution.

## Bioinformatics pipeline info
Scripts for processing of raw data all the way to a filtered VCF file can be found in bioinformatics/scripts/process_data. Scripts are numbered in the order they need to be performed.

Scripts for quality control of several of the steps leading up to the filtered VCF file can be found in bioinformatics/scripts/quality_control. Some of these are run individually, while some are called wihin the scripts for VCF filtering.

Summary statistics of each of the steps leading up to the filtered VCF can be found in bioinformatics/pipeline_summary_statistics.xlsx.

Script for the analyses and visualisation of the results can be found in bioinformatics/scripts/analyses. These are preformed on the filtered VCF file.

Associated metadata for several of the scripts to run can be found in bioinformatics/metadata.
