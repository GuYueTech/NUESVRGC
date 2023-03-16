Created at 3/4/2015 by Dimitris Kugiumtzis

The files in this directory are developed for the computation of the Partial Mutual Information on Mixed Embedding (PMIME), which measures the information flow among all variables from a multivariate time series. There are two main Matlab functions for PMIME:
- PMIME.m : the PMIME using a hard threshold for the termination of the algorithm building the mixed embedding vector. 
- PMIMEsig.m : the PMIME using an adjusted threshold for the termination of the algorithm building the mixed embedding vector. The adjusted threshold is implemented as a significance test making use of randomized (surrogate) data.

As from 3/4/2015, two functions corresponding to PMIME.m and PMIMEsig.m are included, namely:
- PMIMEExact1stLag.m
- PMIMEsigExact1stLag.m
The input and output of these functions are as for PMIME.m and PMIMEsig.m, respectively. These functions differ from PMIME.m and PMIMEsig.m, respectively, only in the first embedding cycle, and they include an exact randomization significance test for the first selected lagged variable. This makes the function slower at the first embedding cycle (with a factor K, where K is the number of time series). They should be used only in the case the time series appear to be random (no auto-correlation and cross-correlation), as for example in returns time series in finance. Otherwise, the user should rely on PMIME.m and PMIMEsig.m. 

The Matlab script "rPMIME.m" is an example for running PMIME.m and PMIMEsig.m on a multivariate time series from coupled Henon maps. The function for the generation of the time series from this map is given in the function "causalhenonmaps2.m". The rest of the files are auxiliary to PMIME.m and PMIMEsig.m. 

The Matlab files are distributed according to the GNU General Public License (see http://www.gnu.org/licenses/ ).

Please consider the following two references if the codes are used for reported results:
1. D. Kugiumtzis, "Direct coupling information measure from non-uniform embedding", Physical Review E, Vol 87, 062918, 2013
2. I. Vlachos, D. Kugiumtzis, "Non-uniform state space reconstruction and coupling detection", Physical Review E, Vol 82, 016207, 2010




