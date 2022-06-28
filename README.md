# InpaintingNMF

This is the accompanying repository for the article *Algorithms for audio inpainting based on probabilistic nonnegative matrix factorization* authored by Ondřej Mokrý, Paul Magron, Thomas Oberlin and Cédric Févotte, submitted to Elsevier Signal Processing.

> Audio inpainting, i.e., the task of restoring missing or occluded audio signal samples, usually relies on sparse representations or autoregressive modeling. In this paper, we propose to structure the spectrogram with nonnegative matrix factorization (NMF) in a probabilistic framework. First, we treat the missing samples as latent variables, and derive two expectation–maximization algorithms for estimating the parameters of the model, depending on whether we formulate the problem in the time- or time-frequency domain. Then, we treat the missing samples as parameters, and we address this novel problem by deriving an alternating minimization scheme. We assess the potential of these algorithms for the task of restoring short- to middle-length gaps in music signals. Experiments reveal great convergence properties of the proposed methods, as well as competitive performance when compared to state-of-the-art audio inpainting techniques.

The repository includes the MATLAB source codes needed to reproduce the research, as well as some supplementary data and figures. As one of the state-of-the-art methods used for comparison, we used the codes for [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).

Note that the psychoacoustically motivated metrics (PEMO-Q, PEAQ) are not computed inside the main testing file [inpainting_comparison.m](https://github.com/ondrejmokry/InpaintingNMF/blob/main/inpainting_comparison.m). This evaluation could be performed by
1. uncommenting the lines related to signal saving in [inpainting_comparison.m](https://github.com/ondrejmokry/InpaintingNMF/blob/main/inpainting_comparison.m) and running the whole experiment,
2. running the file [inpainting_comparison_add_odgs.m](https://github.com/ondrejmokry/InpaintingNMF/blob/main/inpainting_comparison_add_odgs.m).

The PEAQ package is included in the [utils](https://github.com/ondrejmokry/InpaintingNMF/tree/main/utils) subfolder (acquired from [TSP Lab of McGill University](http://www-mmsp.ece.mcgill.ca/Documents/Software/). However, the PEMO-Q software is no longer publicly available, thus this functionality is limited.
