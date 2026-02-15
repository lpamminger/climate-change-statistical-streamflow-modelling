# Source code for manuscript: Climate Change Alters the Fraction of Rainfall That Becomes Streamflow

## Overview

This code links climate change to shifts in streamflow using the rainfall-runoff relationship. For details see:

Link to manuscript here (currently unavailable)

## Repo Contents

This repo contains the data and source code to accompany the manuscript.

## System Requirements

### Hardware requirements

The code in repository requires only a standard computer with enough RAM to support the in-memory operations. The code was developed using a computer with the following specifications:

-   RAM: 32 Gb

-   CPU: 11th Gen Intel (R) Core (TM) i7-11850H \@ 2.5 GHz (8 cores, 16 threads)

### OS Requirements

The code was developed and tested in both Windows and Linux operating systems: The versions used are:

-   Windows 10

-   Linux Ubuntu 24.04

### R Requirements

The code is designed for R version 4.5.2.

The code requires the following R packages: tidyverse, cmaesr, smoof, truncnorm , sloop, ggExtra, gridExtra, furrr, parallel, arrow, tictoc, dream, patchwork, ozmaps, sf, patchwork, metR, ggmagnify, ggh4x, trend and qpdf.

## Installation Guide

To run on local machine either download (size \~ 650 Mb) or clone the repo.

Install time dependent on internet speed and existing packages installed. For typical computer install time is less than 10 minutes.

## Examples

Examples of how to use the code can be found in the `vignette.R` file.

The vignette provides an example of how to fit a rainfall-runoff model to a catchment.

The expected run time for the vignette is less than 30 seconds.
