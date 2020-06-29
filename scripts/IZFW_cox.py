"""
File created 29th June 2020
Authors: Laura Iacovissi, Federico Matteo
"""

import numpy as np
import pandas as pd
from numba import njit
from scipy.sparse.linalg import norm

def e(i, d):
    ei = np.zeros(d)
    ei[i] = 1
    return ei


def IRDSA(F, w, m, c, d):
    """
    Improvised Random Direction stochastic approximation
    for gradient estimation

    INPUT:
    - F: objective function
    - w: current weight
    - m: sample dimension
    - d: features dimension
    - c: costant

    """
    z = np.random.normal(0, 1, (d, m))
    F_w = F(w)
    return np.mean([(F(w + c * z[:,i]) - F_w) / c * z[:,i] for i in range(m)], axis = 0)


def InexactUpdate(g, d, v, r, gamma, mu):
    """
    INPUT
    - g: gradient approximation
    - d: dimension
    - v: starting point
    - r: radius
    - gamma: decreasing coefficient
    - mu: threshold
    """

    haty = v
    t = 1
    while True:
        # ARGMIN PROBLEM
        ht1 = g + gamma*(haty - v)
        i_k = np.argmax(np.abs(ht1))
        ei = e(i_k, d) * r
        yt = np.sign(-ht1[i_k]) * ei
        if np.dot(ht1, yt - haty) >= - mu:
            break
        else:
            haty = (t-1)/(t+1) * haty + 2/(t+1)*yt
            t +=1
    return haty


def IZFW(F, d, w0, L, B = 1, r = 1, T = 100, eps = 1e-6):
    """
    INPUT
    - F: loss function
    - d: dimension
    - w0: starting point
    - L: lipschitz
    - B: 1
    - r: radius of the ball
    - T: max iteration
    - eps: tolerance
    """

    alpha = lambda t: 2/(t+2)
    gamma = lambda t: 4*L/t
    mu = lambda t: L*2*r/(t*T)
    m = lambda t: 100 #t * (t+1) / 2*r * np.max([(d+5)*B*T, d+3])
    c = 1 / (np.sqrt(2*T)) * np.max([1/(d+3), np.sqrt(2*r/(d*(T+1)))]) # smoothing parameter now fixed

    loss = []
    v, w = w0, w0
    partial = 0

    for t in range(1, T+1):
        dt = (1-alpha(t)) * w + alpha(t) * v
        g = IRDSA(F, dt, int(np.ceil(m(t))), c, d)
        v = InexactUpdate(g, d, v, r, gamma(t), mu(t)) #ICG
        w_pred = w
        w = (1 - alpha(t)) * w + alpha(t) * v
        partial += w
        F_w = F(w)
        F_values.append(F_w)
        loss_eval = np.abs(F_values[-2] - F_w)
        loss.append(loss_eval)
        print(f"Loss evaluation at time {t}:\t{loss_eval:.4f}\n")
        if loss_eval < eps: break # check stopping condition
    return F(w_pred), F(w), w, partial/T, t, loss, F_values

@njit
def F(w):
    output = 0
    for i in range(X.shape[0]):
        if y[i] == 1:
            sum_jR = np.sum(np.exp((X @ w)[i:]))
            output += y[i]*(-X[i,:] @ w + np.log(sum_jR))
    return 1/X.shape[0] * output

if __name__ == "__main__":
    # define global variables for data
    global X, y, time, n, d

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

    n, d = X.shape

    # initialize prarameters for the algorithm:
    # stating point
    w0 = np.random.rand(d)
    w0 = w0/np.sum(w0) * np.random.rand(1) *10
    # Lipschitz constant
    L = np.linalg.norm(X.T @ X)

    # call ZFW with InexactUpdate
    fpred, f, w, mean, t, loss, f_values = IZFW(F, d, w0, L, B = 1, r = 10, T = 100, eps = 1e-8)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))

    np.save('../Data/results/function_IZFW_cox.npy',f_values)
    np.save('../Data/results/loss_IZFW_cox.npy',loss)
