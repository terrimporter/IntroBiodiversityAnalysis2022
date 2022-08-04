# Introduction to Biodiversity Analysis

This repository contains a powerpoint presentation (pdf) as well as some scripts to explore biodiversity data and create some typical biodiversity plots.

This is not intended to be comprehensive, but provide a gentle start towards exploring biodiversity datasets.

The examples shown in the powerpoint slides come from taxonomically assigned exact sequence variants (ESVs) from MetaWorks (results.csv).  The scripts can be adapted to accomodate different file naming conventions (SampleName) and output from different pipelines (ESV table, taxonomically assigned sequence clusters).

*Before running scripts, create a new R project directory and put scripts/ into it.  In the R project directory,  create a directory called infiles/ and put your results.csv file in it.  Also create a directory called outfiles/ to contain the output.*

```bash
├── R project directory
│   ├── scripts
│   ├── infiles
│   │   ├── results.csv
│   ├── outfiles
```

Last updated: August 4, 2022
