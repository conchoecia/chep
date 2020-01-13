# chep
Calculate Heterozygosity with Pileup



# Contents
- Introduction
- Installation
- Usage
- License
- Citation

# Introduction

This program is useful to:
  - calculate and visualize the heterozygosity of a genome assembly
  - to visualize ploidy
  - and to determine if the genome is sufficiently haplotype-collapsed.
  
It works by scanning a bam file with `mpileup` and asking the questions, "What is the read depth at this position, and how many bases match the reference base?" This information is turned into a 3D histogram, where:
  - Column 0 is the read depth
  - Column 1 is how many reads had the reference allele at that position
  - Column 3 is the number of positions in the genome that match Columns 0 and 1
  
This 3D histogram is then used to generate plots of the marginal histogram of read depths across the whole genome (top panel), the heterozygosity calculated at every specific read depth (middle panel), and a plot of the underlying 3D histogram.


  
# Installation

To install, go to a directory where you want to install the program. Execute:

```
git clone https://github.com/conchoecia/chep.git
cd chep
make
```

Currently, `chep` doesn't automatically install the scripts into a system-wide location. As a result you must add the following line to your `~/.bash_profile` or `~/.bashrc` file.

```
export PATH=$PATH:/home/dschultz/chep/bin
```

# Usage

# License

Chep is free software and is licensed under GPLv3.

# Citation

Currently there is no way to cite this software aside from citing this github page.
