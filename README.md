[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14911032.svg)](https://doi.org/10.5281/zenodo.14911032)

This repository contains the data and codes from the following manuscript:

**Tall stature and small leaves: ecological strategies that enhance tree growth across the subtropical Brazilian Atlantic Forest**

Repository content:

The `scripts` folder contains:

Scripts to perform all the analyses from the manuscript\
*code.R* - packages and functions, data visualisation, PCAs and linear models.

The `data` folder contains:

Functional trait data and relative growth rates for all the 121 species\
species.information.csv - csv. file containing data, see description of each column below

plot.rtf: plot and census information

metadata.xlxs: metadata (readme file) with the same information as described below

-   Species: species name

-   Family: Species family

-   Mean.absolute.growth.rate: mean species growth rate based on diameter (mm) in time 1 minus diameter in time 0, divided by census interval (in years) \<-\> (dbh1-dbh0)/time

-   SD.growth.rate: standard deviation from the mean absolute growth rates

-   SLA: specific leaf area, in cm2/g

-   LA: leaf area, in cm2

-   LDMC: leaf dry-matter content, in mg/g

-   LNC: leaf nitrogen content, in %

-   LPC: leaf phosphorus content, in %

-   WD: wood density, in g/cm3

-   Height: maximum height per species within the sampled trees, in meters

-   Number.of.stems: number of stems per species

-   Median.height: species median height within the sampled trees, in meters

-   Category: canopy or understory classification

-   Relative.Growth.Rate: species relative growth rate, based on the following equation: (dbht1/dbhto)/t

The `results` folder contains:

*principal.components.png*: PCAs diagrams

The code was developed and tested using the R statistical software, R version 4.2.2 (R Development Core Team, 2022).
