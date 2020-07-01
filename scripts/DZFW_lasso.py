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


def KWSA(F, w, m, c, d):
    """
    Kiefer-Wolfowitz stochastic approximation
    for gradient estimation

    INPUT:
    - F: objective function
    - w: current weight
    - m: sample size (null in this case)
    - d: dimension
    - c: costant

    """

    F_wc = np.array([F(w + c * e(i, d)) for i in range(d)])
    return (F_wc - F(w)) / c

def detZFW(F, L, d, w0, r=1, T=100, eps=1e-5):
    """
    INPUT
    - F: loss function
    - L: Lip constant
    - d: dimension
    - w0: starting point
    - r: radius of the ball
    - T: max iteration
    - eps: tolerance
    """

    gamma = lambda t: 2/(t+2)
    c = lambda t: L*gamma(t)/d
    w = w0
    partial = 0
    loss = []
    F_values = [F(w0)]
    for t in range(1, T+1):
        # comupute the gradient approx
        gt = KWSA(F, w, None, c(t), d)
        # compute the linear problem solution on the L1 Ball of radius r
        i_k = np.argmax(np.abs(gt))
        ei = e(i_k, d) * r
        v = np.sign(-gt[i_k]) * ei
        # compute step
        w_pred = w
        w = (1 - gamma(t)) * w + gamma(t) * v
        partial += w
        F_w = F(w)
        F_values.append(F_w)
        loss_eval = np.abs(F_values[-2] - F_w)
        loss.append(loss_eval)
        print(f"Loss evaluation at time {t}:\t{loss_eval:.7f}\n")
        if loss_eval < eps: break # check stopping condition
    return F(w_pred), F(w), w, partial/T, t, loss, F_values


if __name__ == "__main__":
    # define global variables for data
    global X, y, d

    # set random seed
    np.random.seed(1007)

    # load data
    X, y = datasets.load_svmlight_file("../Data/covtype.libsvm.binary.scale.bz2")

    # space dimension
    d = X.shape[1]

    # Lipschitz constant computation
    L = 3

    # define the objective function
    F = lambda w: 0.5/n * np.sum(np.power(y - X @ w, 2))

    # initialize prarameters for the algorithm:
    # stating point
    np.random.seed(1007)
    w0 = np.random.rand(d)
    w0 = w0/np.sum(w0) * np.random.rand(1)

    # call deterministic ZFW
    fpred, f, w, mean, t, loss, f_values = detZFW(F, L, d, w0, r=1, T=1000, eps=1e-8)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))

    np.save('../Data/results/function_DZFW_lasso.npy',f_values)
    np.save('../Data/results/loss_DZFW_lasso.npy',loss)
