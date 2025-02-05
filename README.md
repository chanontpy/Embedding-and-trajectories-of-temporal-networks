# Embedding-and-trajectories-of-temporal-networks

This repository provides the codes for\
[Chanon Thongprayoon, Lorenzo Livi, Naoki Masuda.\
Embedding and trajectories of temporal networks.\
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
- `Adjacency`: to create an adjacency matrix at each time step, where at least one event occurs between a pair of nodes. The input is a list of time-stamped contact events.
- `ComputeBk`: use the tie-decay matrix at from the previous time step to compute the tie-decay matrix of the current time step. The input arguments are the exponential decay rate `alpha` and the array containing adjacency matrices from `Adjacency`.
- `Tie_decay_matrices`: by using `ComputeBk`, this function computes tie-decay matrices of multiple time steps. Note that the index starts at $1$.
- `Find_eigenpair`: compute eigenvalues and eigenvectors of a matrix.
- `Centered_Distance_matrix`: apply double mean-centering matrices, each to the left and the right side of the (squared) distance matrix.
- `Squared_Frobenius_Distance_matrix`: compute the squared Frobenius distance matrix of tie-decay matrices.
- `Degree_matrix`: compute the degree matrix of the input matrix.
- `Laplacian_matrix`: compute the Laplacian matrix of the input matrix.
- `Laplacian_matrix`: compute the Laplacian matrix of the input matrix.
- `Find_eigenpair`: compute eigenvalues and their corresponding eigenvectors.
- `Squared_Laplacian_distance`: compute the squared Laplacian distance between two matrices.
- `Squared_Laplacian_distance_matrix`: construct the squared Laplacian distance matrix.
- `Centered_Distance_matrix`: double mean centering a squared distance matrix.
- `aux_delta`: compute egenvalues and eigenvectors of each landmark.
- `Classical_MDS`: classical multidimensional scaling
- `Col_sum`: finding the vector obtained from averaging the columns of a matrix.
- `squared_dist_vector`: compute a vector encoding a squared distance between a tie-decay matrix and a landmark.
- `LMDS`: compute LMDS.
  ### One may use the following functions to aid the computation of LMDS with the tie-decay property
- `determine_interval`: finding the subinterval in which a particular time step lies in.
- `convenient_sq_laplacian_distance`: combinging tie-decay property with LMDS computation.
