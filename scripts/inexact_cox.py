"""
File created 23th June 2020
Authors: Laura Iacovissi, Federico Matteo
"""

import numpy as np
import pandas as pd
from numba import njit
from modules.ZFW import *
from scipy.sparse.linalg import norm

@njit
def F(w):
    output = 0
    for i in range(n):
        if y[i] == 1:
            sum_jR = np.sum(np.exp((X @ w)[i:]))
            output += y[i]*(-X[i,:] @ w + np.log(sum_jR))
    return 1/n * output

if __name__ == "__main__":
    # set random seed
    np.random.seed(7)

    # load clinical data to extract time of risk
    clinical = pd.read_table('../Data/SurvivalTimes.txt', index_col=0, sep=';')
    clinical = clinical.set_index(clinical["IDs"], drop=True).iloc[:,:-1]
    clinical["death_event"].value_counts()
    # load rna data
    data = pd.read_table("../Data/mydata.txt", sep = ";")
    data = data.T
    # crate datafram with data of interest
    df = data.merge(clinical, left_index= True, right_index = True)
    df = df.sort_values(by = "new_death")
    del data, clinical
    # define feature matrix X, label vector y and vector of times
    X, y, time = df.iloc[:,:-2], df.iloc[:,-1].values, df.iloc[:,-2].values
    X = np.array(X.apply(lambda x: (x - np.mean(x))/np.std(x), axis = 0))
    del df

    global X, y, time

    # initialize prarameters for the algorithm:
    # stating point
    w0 = np.random.rand(d)
    w0 = w0/np.sum(w0) * np.random.rand(1) *10
    # Lipschitz constant
    L = 2/X.shape[0] * np.linalg.norm(X.T @ X)

    # call ZFW with InexactUpdate
    fpred, f, w, mean, t, loss = stochasticZFW(F, d, w0, method = "IRDSA", r=10, T=40, eps=1e-3)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))
