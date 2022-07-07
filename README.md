# Weak population genetic structure in Eurasian spruce bark beetle over large regional scales in Sweden

<a href="https://doi.org/10.1002/ece3.9078"><img src="https://img.shields.io/badge/Ecology%20and%20Evolution-doi%3A%2010.1002%2Fece3.9078-%23440154"></a>

This repository contains lab protocols, RADtag library preperation summaries, and bioinformatic scripts used in:

Ellerstrand, S. J., Choudhury, S., Svensson, K., Andersson, M., Kirkeby, C., Powell, D., Schlyter, F., Jönsson, A. M., Brydegaard, M., Hansson, B., & Runemark, A. (2022). Weak population genetic structure in Eurasian spruce bark beetle over large regional scales in Sweden. Ecology and Evolution, 12, e9078. https://doi.org/10.1002/ece3.9078

![Figure 1](https://user-images.githubusercontent.com/81574558/177731550-8c86e88b-d173-4900-b5b3-8becac9574e1.jpg)

## Bioinformatics pipeline info
Scripts for processing of raw data all the way to a filtered VCF file can be found in bioinformatics/scripts/process_data. Scripts are numbered in the order they need to be performed.

Scripts for quality control of several of the steps leading up to the filtered VCF file can be found in bioinformatics/scripts/quality_control. Some of these are run individually, while some are called wihin the scripts for VCF filtering.

Summary statistics of each of the steps leading up to the filtered VCF can be found in bioinformatics/pipeline_summary_statistics.xlsx.

Script for the analyses and visualisation of the results can be found in bioinformatics/scripts/analyses. These are preformed on the filtered VCF file.

Associated metadata for several of the scripts to run can be found in bioinformatics/metadata.
