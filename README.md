# Embedding and trajectories of temporal networks

This repository provides the codes for\
[Chanon Thongprayoon, Lorenzo Livi, Naoki Masuda.
Embedding and trajectories of temporal networks.
IEEE Access, 11, 41426-41443 (2023)](http://doi.org/10.1109/ACCESS.2023.3268030).

## Python package dependencies
Run
```R
`import pandas as pd`
`import numpy as np`
`import numpy.linalg as LA`
`import math`
`import matplotlib.pyplot as plt`
```

## Functions
- `Adjacency`: Create an adjacency matrix at each time step where at least one event occurs between a pair of nodes. The input is a list of time-stamped contact events.
- `ComputeBk`: Use the tie-decay matrix at from [NM: "at from" sounds wrong. At least I am confused. Rephrase, please.] the previous time step to compute the tie-decay matrix of the current time step. The input arguments are the exponential decay rate `alpha` and an array containing adjacency matrices obtained from `Adjacency`.
- `Tie_decay_matrices`: Using `ComputeBk`, this function computes tie-decay matrices of multiple time steps. Note that the index [NM: index of what? To avoid this (and potentially other) confusion, it is better to say explicitly what the inputs are and what the outputs are.] starts at $1$.
- `Find_eigenpair`: Compute eigenvalues and eigenvectors of a matrix.
- `Centered_Distance_matrix`: Apply double mean centering, i.e., apply matrix $$H_M$$ to the left and right, each of the (squared) distance matrix. The input is the square distance matrix [NM: I added this sentence, but I am not sure whether this is correct. Please check. If correct, you can remove this comment. Otherwise, please correct.]
- `Squared_Frobenius_Distance_matrix`: Compute the squared Frobenius distance matrix of tie-decay matrices. [NM: Again, input?]
- `Degree_matrix`: Compute the degree matrix of the input matrix. [NM: Degree matrix, not the degree vector? Does it mean that the output is a diagonal matrix having the degree? Please be explicit about the ouptut.]
- `Laplacian_matrix`: Compute the Laplacian matrix of the input matrix.
- `Laplacian_matrix`: Compute the Laplacian matrix of the input matrix [NM: Why do you repeat this twice? Do you mean a different function, or otherwise this line should be removed?]. 
- `Find_eigenpair`: Compute eigenvalues and the associated corresponding eigenvectors. The input is [NM: Fill.].
- `Squared_Laplacian_distance`: Compute the squared Laplacian distance between two matrices. The inputs are [NM: Fill.].
- `Squared_Laplacian_distance_matrix`: Construct the squared Laplacian distance matrix [NM: I do not understand the difference between this one and the last one.].
- `Centered_Distance_matrix`: Double mean centering a squared distance matrix [NM: We did it above. Is this function different? Or it is duplicated so this line should be removed?].
- `aux_delta`: Compute egenvalues and eigenvectors of each landmark. The input is [NM: Fill.]
- `Classical_MDS`: Run the classical multidimensional scaling. [NM: Input? Also the output?]
- `Col_sum`: Calculate the vector obtained by averaging each column of a matrix. The input is [NM: Fill.]
- `squared_dist_vector`: Compute a vector encoding a squared distance between a tie-decay matrix and a landmark. [NM: input? Why output is a vector? To avoid that confusion, it is also better to say the output explicitly for this one.]
- `LMDS`: Compute the LMDS. [NM: input and output?]
  
  ### One may use the following functions to aid the computation of LMDS with the tie-decay property

- `determine_interval`: Find the subinterval in which a particular time step lies in. [NM: input and output?]
- `convenient_sq_laplacian_distance`: Combine the tie-decay property with LMDS computation. [NM: Vague. What does this function do concretely?}

## Example usage

[NM: Is it possible to show example usage? The above is too difficult for people to use. What people ideally want to see is a working example. So, ideally, there is a mock data (i.e., time-stamped event sequence), which is also uploaded on this repository. Using that, how can one run a sequence of (?) those commands to get some quantities? Just one example suffices.]
