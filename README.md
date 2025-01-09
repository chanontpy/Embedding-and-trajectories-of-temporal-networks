# Embedding-and-trajectories-of-temporal-networks

This repository provides the codes for the article [Embedding and trajectories of temporal networks](http://doi.org/10.1109/ACCESS.2023.3268030)

## Python packages
- `import pandas as pd`
- `import numpy as np`
- `import numpy.linalg as LA`
- `import math`
- `import matplotlib.pyplot as plt`

## Functions
- `Adjacency`: to create an adjacency matrix at each time step, where at least one event occurs between a pair of nodes. The input is a list of time-stamped contact events.
- `ComputeBk`: use the tie-decay matrix at from the previous time step to compute the tie-decay matrix of the current time step.
- `Tie_decay_matrices`: by using `ComputeBk`, this function computes tie-decay matrices of multiple time steps. Note that the index starts at $1$.
- `Find_eigenpair`: compute eigenvalues and eigenvectors of a matrix.
- `Centered_Distance_matrix`: apply double mean-centering matrices, each to the left and the right side of the (squared) distance matrix.
- `Squared_Frobenius_Distance_matrix`: compute the squared Frobenius distance matrix of tie-decay matrices.
- `Degree_matrix`: compute the degree matrix of the input matrix.
- `Laplacian_matrix`: compute the Laplacian matrix of the input matrix.
