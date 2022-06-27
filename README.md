# InpaintingNMF

This is the accompanying repository for the article *Algorithms for audio inpainting based on probabilistic nonnegative matrix factorization* authored by Ondřej Mokrý, Paul Magron, Thomas Oberlin and Cédric Févotte, submitted to Elsevier Signal Processing.

> Audio inpainting, i.e., the task of restoring missing or occluded audio signal samples, usually relies on sparse representations or autoregressive modeling. In this paper, we propose to structure the spectrogram with nonnegative matrix factorization (NMF) in a probabilistic framework. First, we treat the missing samples as latent variables, and derive two expectation–maximization algorithms for estimating the parameters of the model, depending on whether we formulate the problem in the time- or time-frequency domain. Then, we treat the missing samples as parameters, and we address this novel problem by deriving an alternating minimization scheme. We assess the potential of these algorithms for the task of restoring short- to middle-length gaps in music signals. Experiments reveal great convergence properties of the proposed methods, as well as competitive performance when compared to state-of-the-art audio inpainting techniques.

The repository includes the MATLAB source codes needed to reproduce the research, as well as some supplementary data and figures. As one of the state-of-the-art methods used for comparison, we used the codes for [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).

**The source codes are currently being reviewed and polished before publishing them, everything should appear online in the next few days.**
