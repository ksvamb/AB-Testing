# Sample Size Calculations for A/B Testing: Novel Extensions & Practical Guide_

This repository hosts the implementation of section 3 - "Sample Size for Correlated Data" from the paper _"All about Sample-Size Calculations for A/B Testing: Novel Extensions & Practical Guide"_ by Jing Zhou, Jiannan Lu, Anas Shallah, published in August 2023. The full paper is available [here](https://arxiv.org/abs/2305.16459).

## Overview

The provided code focuses on the Delta method, which is suggested for calculating sample sizes in A/B testing when the unit of analysis is smaller than the unit of randomization. Zhou et al. (2023) argue that this approach enhances the assessment of statistical significance, specifically Type I error and power, compared to traditional sample size methods used for independent data.

This repository not only replicates the original comparative analysis but also extends it, focusing on the absolute lift across eleven different scenarios. Our findings support the conclusions drawn by Zhou et al.

The original paper can be accessed from the PDF file named "Paper" within this repository. For a detailed explanation of the code and further discussion, please refer to the "Discussion" document.
