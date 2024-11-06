# A map continuity view of assessing and improving the reliability of neighbor embedding methods
We provide two diagnostic scores to assess and improve the reliability of neighbor embedding methods based on a map-continuity view. The perturbation score is to diagnose points around the location of overconfidence-inducing discontinuity and the singularity score is to diagnose points of fracture-inducing discontinuity. Our approach is flexible and works as a wrapper around many neighbor embedding algorithms. The R package is now available for t-SNE.

The method is based on the paper:

Liu, Z., Ma, R., Zhong, Y. (2023) Assessing and improving reliability of neighbor embedding methods: a map-continuity perspective. arXiv:2410.16608.

## Contents
The directory `neMDBD` contains the R package `neMDBD` for implementation of perturbation score and singularity score.

The directory `RtsneWithP` contains the modified R package `Rtsne` with the similarity matrix as an extra output.

The directory `Code and data` contains the R scripts for reproducing the data analysis and examples in the manuscript.

## Get started
For the quick guide to the two diagnostic scores, please check out https://github.com/zhexuandliu/NE-Reliability-MapContinuity/blob/main/neMDBD/README.md.

To find out how to extract the similarity matrix computed in `Rtsne`, please check out https://github.com/zhexuandliu/NE-Reliability-MapContinuity/blob/main/RtsneWithP/README.md.

To check how to draw the LOO loss landscape and other implementation details, please check out the folder `Code and data`.

For further questions and inquiries, please contact Zhexuan Liu (zhexuan.liu2@wisc.edu).
