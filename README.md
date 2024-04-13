# Sample-Size Calculations for A/B Testing: Novel Extensions & Practical Guide

This repository contains the implementation of the initial part of the methodology presented in the paper "All about Sample-Size Calculations for A/B Testing: Novel Extensions & Practical Guide" by Jing Zhou, Jiannan Lu, Anas Shallah, published in August 2023. The paper can be accessed [here](https://arxiv.org/abs/2305.16459).

## Overview

The code provided in this repository focuses on the innovative approaches for calculating sample sizes in A/B testing as discussed in the paper. These methods are designed to enhance the precision and efficiency of statistical tests in marketing analytics and other fields. Specifically, it includes comprehensive code for the "Sample Size for Correlated Data" section of the paper, and the simulation results here match those obtained by the authors.

This project evaluates the method suggested by Zhou et al. (2023) for estimating sample size when the unit of analysis is smaller than the unit of randomization. To account for data correlation, Zhou et al. recommend using the Delta Method, arguing that it improves statistical significance in assessments of Type I error and power compared to the "standard" sample size method used for independent data. This study replicates and extends the original comparative analysis, focusing on absolute lift in eleven different scenarios. Our findings support those of Zhou et al.

You can view the original paper from the PDF file named "Paper" in this repository. A detailed code explanation and discussion can be found in the "Discussion" document.
