# Entanglement-Moments


The objective of this work is to show that the use of PnCP maps coupled to that of Matrices of MoMents, explicitly coded or via an SDP formulation enables the detection of useful and hard to detect Bound Entangled (BE) States.

We load PPT BE found either in the litterature or constructed by us. We check that their entanglement is revealed by Abhishek algorithm and not by PPT criteria.We apply noise models on them, and check the same it again.

Now the objective is to show that PnCP algorithm makes the most of the information gathered by measurements.

We construct different kind of Matrices of Moments and see if the application of PnCP algorithm enables us to discover entanglement better. Two parameters play a part: the noisier the state is, the more difficult it should be to detect entanglement, the more measurements we have the easier. So our technique must be compared with regard to 1) the order of the MM 2) the threshold value of noise.

This construction is particularly handy for big dimensions because it is computationally not demanding. However, it is not systematic enough. 

The other approach is to implement an SDP that maximizes any quantity, with constraints on the measurements of a set of moments, and the fact that the underlied state must be positive, hermitian, trace 1 and PPT. This consists in finding all the possible matrices of moments and see if a separable state can reproduce these results. We replace the PPT assumption by the fact that the state must be positive under PnCPs in general, and use a library of PnCPs obtained through Abhishek algorithm. 

