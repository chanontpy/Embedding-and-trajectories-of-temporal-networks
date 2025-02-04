import pandas as pd
import numpy as np
import numpy.linalg as LA
import math
import matplotlib.pyplot as plt

def Adjacency(C):
    Net =[[ [], -1.0 ]]
    Net2 =[[ [], -1.0 ]]
    time_points = C[:,2]
    nodes0 = C[:,0]
    nodes1 = C[:,1]
    num1 = int(max(np.unique(nodes0)))
    num2 = int(max(np.unique(nodes1)))
    num = max([num1,num2])
    N = np.shape(C[:,0])[0]
    time_points = np.unique(C[:,2])
    
    for t in time_points:
        A = np.zeros((num+1,num+1))
        for k in range(N):
            if C[k][2] == t:
                n = int(C[k][0])
                m = int(C[k][1])
                A[n][m] = 1
                A[m][n] = 1
        A2 = A.tolist()
        Net.append([A2,t])
    return Net

def ComputeBk(net_k, t_k, alpha, k ,B_k_1,t_k_1):
    if k == 1 :
        net1  = net_k
        t1 = t_k
        B1 = np.copy(net1)
        return B1, t1

    numrow = len(net_k)
    numcol = len(net_k[0])
    Bk = np.zeros((numrow, numcol))
    for i in range(0, numrow):
        for j in range(0,numcol):
            node_k   = net_k[i][j]
            if node_k != 0:
                update = node_k + ( B_k_1[i][j]/math.exp( alpha*(t_k-t_k_1) ) )
            else:
                update = 0 + ( B_k_1[i][j]/math.exp( alpha*(t_k-t_k_1) ) )

            Bk[i][j] = update
    return Bk,t_k

def Tie_decay_matrices(alpha, C):
    Net = Adjacency(C)
    B = [[[],-1.0]]
    for i, n_i in enumerate(Net):
        if i==0 :
            continue
        Bi_1 = B[i-1]
        Bi,Ti = ComputeBk(n_i[0],n_i[1] ,alpha, i, Bi_1[0],Bi_1[1])
        B.append([Bi,Ti])
    return B

def Squared_Frobenius_Distance_matrix(B): #B is from tie_decay_matrices(alpha,C)
    D = np.zeros((len(B),len(B)))
    for i, Bi in enumerate(B):
        for j, Bj in enumerate(B):
            if i == j:
                D[i][j] = 0
                D[j][i] = 0
            elif i < j:
                D[i][j] = LA.norm(Bi-Bj)**2
                D[j][i] = D[i][j]
            else:
                continue
    return D

def Degree_matrix(A): #A is a tie-decay adjacency matrix at a time t.
    m = np.shape(A)[0]
    n = np.shape(A)[1]
    if m!=n:
        return print('Error: not a square matrix')
    else:
        degree_matrix = np.zeros((m,n))
    
        for i in range(0,m):
            r = sum(A[i])
            degree_matrix[i][i] = degree_matrix[i][i] + r

        return degree_matrix

def Laplacian_matrix(A): #A is also a tie-decay adjacency matrix at a time t.
    L = Degree_matrix(A) - A
    return L

def Find_eigenpair(D):
    A = []
    eigenvalues = LA.eig(D)[0]
    eigenvectors_col = LA.eig(D)[1]
    
    idx = eigenvalues.argsort()[::-1]   
    eigenvalues = eigenvalues[idx]
    eigenvectors_col = eigenvectors_col[:,idx]
    
    eigenvectors_row = np.transpose(eigenvectors_col)
    for i in range(len(eigenvalues)):
        A.append([eigenvalues[i],eigenvectors_row[i]])
    return A

def Squared_Laplacian_distance(A,B): # A,B are tie-decay matrices. Already squared.
    L = Laplacian_matrix(A)
    K = Laplacian_matrix(B)
    M, N = Find_eigenpair(L), Find_eigenpair(K)
    m = len(M)
    #n = len(N)
    value = 0
    
    for i in range(0,m):
        value = value + ((np.real(M[i][0]) - np.real(N[i][0]))**2)
        
    #value = np.sqrt(value)
    return value

def Squared_Laplacian_distance_matrix(B): # B = Tie_decay_matrices(alpha, C)
    LD = np.zeros((len(B),len(B)))
    
    for i, Bi in enumerate(B):
        for j, Bj in enumerate(B):
            if i == j:
                LD[i][j] = 0
                LD[j][i] = 0
            elif i < j:
                LD[i][j] = Squared_Laplacian_distance(Bi,Bj)
                LD[j][i] = LD[i][j]
            else:
                continue
            
    return LD

def Centered_Distance_matrix(D): # D is a squared distance matrix
    n = np.shape(D)[0]
    m = np.shape(D)[1]
    I = np.identity(n)
    One_mat = np.ones((n,n))
    H = I - (1/n * One_mat)
    D_scaled = -0.5 * np.matmul(H,np.matmul(D,H))
    return D_scaled

def aux_delta(list_of_landmarks): # list_of_landmarks is a list containing landmarks
    P = []
    X = []
    
    for L in list_of_landmarks:
        M = Find_eigenpair(L)
        m = len(M)
        Q = []
        Y = []
        for k in range(0,m):
            Q.append(M[k][0])
            Y.append(M[k][1])
        P.append(Q)
        X.append(Y)
    
    return P,X

def Classical_MDS(D,k): # D is a squared distance matrix, k is the embedding dimension
    L = []
    #D = normalisation(D)
    B = Centered_Distance_matrix(D)
    #B = 1/np.std(B_ast) * B_ast
    Lambda = Find_eigenpair(B) # Recall that the eigenvalues are already sorted.
    total_sum = 0
    neg_sum = 0
    confidence = 0
    sqrt_eigenvalue = []
    for i in range(0,len(Lambda)):
        if Lambda[i][0] > 0:
            sqrt_eigenvalue.append(np.sqrt(Lambda[i][0]))

    for i in range(0,k):
        L.append(sqrt_eigenvalue[i] * Lambda[i][1])

    return L, Lambda

def Col_sum(D): # D is a squared distance matrix
    N = np.shape(D)[1]
    Col_D_sum = np.reshape(np.sum(D, axis = 1), (N,1))
    Delta_mu = 1/N * Col_D_sum
    
    return Delta_mu

#---LMDS starts here---#

a_delta = aux_delta(list_of_landmarks)
landmark_eigenvalues = a_delta[0]
M = len(landmark_eigenvalues)
Laplacian_dist_mat = np.zeros((M,M))
Landmark_eigenpair = Classical_MDS(Laplacian_dist_mat,2)[1] #These eigenpairs already dependends on the choice of distance measure.
delta_mu = Col_sum(Laplacian_dist_mat)

N = np.shape(Laplacian_dist_mat)[1]
inv_sqrt_eigenvalues = []
L_prime_before = []
k = 2
for i in range(0,k):
    if Landmark_eigenpair[i][0] > 0:
        inv_sqrt_eigenvalues.append(1 / np.sqrt(Landmark_eigenpair[i][0]))
    
for i in range(0,k):
    vi = inv_sqrt_eigenvalues[i] * Landmark_eigenpair[i][1]
    L_prime_before = np.concatenate((L_prime_before,vi.tolist()), axis = 0)
        
L_prime = np.reshape(L_prime_before, (k,N)) # Pseudo inverse matrix of L

def squared_dist_vector(t,landmark_eigenvalues,dictionary):
    # "dictionary" is a python dictionary whose keys values are time steps 't' and tie-decay matrices at time 't'.
    eigen_pair_arb = Find_eigenpair(Laplacian_matrix(dictionary[t]))
    eigen_val_arb = np.array([i[0] for i in eigen_pair_arb])
    sq_dist_vec = []
    for i in landmark_eigenvalues:
        sq_dist_vec.append(LA.norm(eigen_val_arb-np.array(i))**2)
    return np.array(sq_dist_vec)

def LMDS(t):
    sq_vec = squared_dist_vector(t,landmark_eigenvalues,dictionary)
    sq_vec = sq_vec.reshape(N,1)
    vec = sq_vec - delta_mu
    embed_coord =  -0.5*np.matmul(L_prime,vec)
    return embed_coord
