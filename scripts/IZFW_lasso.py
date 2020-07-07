"""
File created 29th June 2020
Authors: Laura Iacovissi, Federico Matteo
"""

import numpy as np
import pandas as pd
from sklearn import datasets


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
    return haty, t


def IZFW(F, d, w0, L, B = 1, D = 2, r=1, T = 100, eps = 1e-6):
    """
    INPUT
    - F: loss function
    - d: dimension
    - w0: starting point
    - L: lipschitz
    - B: 1
    - D: radius estimate from above
    - r: l1 ball radius
    - r: radius of the ball
    - T: max iteration
    - eps: tolerance
    """

    alpha = lambda t: 2/(t+2)
    gamma = lambda t: 4*L/t
    mu = lambda t: L*D/(t*T)
    m = lambda t: t * (t+1) / D * (d+3) #np.max([(d+5)*B*T, d+3])
    c = 1 / (np.sqrt(2*T)) * np.max([1/(d+3), np.sqrt(D/(d*(T+1)))]) # smoothing parameter now fixed

    loss = []
    F_values = [F(w0)]
    v, w = w0, w0
    partial = 0
    inner = [0]
    for t in range(1, T+1):
        dt = (1-alpha(t)) * w + alpha(t) * v
        g = IRDSA(F, dt, int(np.ceil(m(t))), c, d)
        v, inner_t = InexactUpdate(g, d, v, r, gamma(t), mu(t)) #ICG
        inner.append(inner_t)
        w_pred = w
        w = (1 - alpha(t)) * w + alpha(t) * v
        partial += w
        F_w = F(w)
        F_values.append(F_w)
        loss_eval = np.abs(F_values[-2] - F_w)
        loss.append(loss_eval)
        print(f"Loss evaluation at time {t}:\t{loss_eval:.7f}\n")
        if loss_eval < eps: break # check stopping condition
    return F(w_pred), F(w), w, partial/T, t, loss, F_values, inner


if __name__ == "__main__":
    # define global variables for data
    global X, y, d

    # set random seed
    np.random.seed(1007)

    # load data
    X, y = datasets.load_svmlight_file("../Data/covtype.libsvm.binary.scale.bz2")

    # space dimension
    d = X.shape[1]
    n = y.shape[0]

    # Lipschitz constant computation
    L = 3
    D = 20 # we will start from m = 6, up to T * (T+1) / D * (d+3) = 28785 (for T=100)
    B = 1

    # define the objective function
    F = lambda w: 0.5/n * np.sum(np.power(y - X @ w, 2))

    # initialize prarameters for the algorithm:
    # stating point
    np.random.seed(1007)
    w0 = np.random.rand(d)
    w0 = w0/np.sum(w0) * np.random.rand(1)

    # call stochastic ZFW with InexactUpdate
    fpred, f, w, mean, t, loss, f_values, inner = IZFW(F, d, w0, L, B, D, T=100, eps=1e-6)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))

    np.save('../Data/results/function_IZFW_lasso_long.npy',f_values)
    np.save('../Data/results/loss_IZFW_lasso_long.npy',loss)
    np.save('../Data/results/inner_IZFW_lasso_long.npy',inner)
