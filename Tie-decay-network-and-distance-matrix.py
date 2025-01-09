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
