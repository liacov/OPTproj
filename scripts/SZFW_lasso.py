"""
File created 23th June 2020
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


def stochasticZFW(F, d,  w0, method = "IRDSA", r=1, T=100, eps=1e-5):
    """
    INPUT
    - F: loss function
    - d: dimension
    - w0: starting point
    - method: zeroth order oracle
    - r: radius of the ball
    - T: max iteration
    - eps: tolerance
    """

    Parameters_dict = {"KWSA": {"m": None,
                                "c": lambda t: 2 / (np.sqrt(d) * np.power(t+8, 1/3)),
                                "p": lambda t: 4 / np.power(t+8, 2/3),
                                "oracle": KWSA},


                       "RDSA": {"m": 1,
                                "c": lambda t: 2 / (np.power(d, 3/2) * np.power(t+8, 1/3)),
                                "p": lambda t: 4 / (np.power(d, 1/3) * np.power(t+8, 2/3)),
                                "oracle": IRDSA},

                       "IRDSA": {"m": 938,
                                "c": lambda t: 2 * np.sqrt(6) / (np.power(d, 3/2) * np.power(t+8, 1/3)),
                                "p": lambda t: 4 / (np.power(1+d/6, 1/3) * np.power(t+8, 2/3)),
                                "oracle": IRDSA}

                        }

    return sZFW(F, d, w0, Parameters_dict[method], r, T, eps)



def sZFW(F, d, w0, params, r, T, eps):
    """
    INPUT
    - F: loss function
    - d: dimension
    - w0: starting point
    - params: dict of parameters for the selected method
    - r: radius of the ball
    - T: max iteration
    - eps: tolerance
    """

    loss = []
    F_values = [F(w0)]
    gamma = lambda t: 2/(t+8)
    w = w0
    dt = np.zeros(d)
    partial = 0
    for t in range(1, T+1):
        # comupute the gradient approx
        gt = params["oracle"](F, w, params["m"], params["c"](t), d)
        dt = (1 - params["p"](t)) * dt + params["p"](t) * gt
        # compute the linear problem solution on the L1 Ball of radius r
        ei = e(np.argmax(np.abs(dt)), d) * r
        v = np.sign(-dt) * ei
        # compute step
        w_pred = w
        w = (1 - gamma(t)) * w + gamma(t) * v
        partial += w
        F_w = F(w)
        F_values.append(F_w)
        loss_eval = np.abs(F_values[-2] - F_w)
        loss.append(loss_eval)
        print(f"Loss evaluation at time {t}:\t{loss_eval:.4f}\n")
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

    # define the objective function
    F = lambda w: 0.5 * np.sum(np.power(y - X @ w, 2))

    # initialize prarameters for the algorithm:
    # stating point
    np.random.seed(1007)
    w0 = np.random.rand(d)
    w0 = w0/np.sum(w0) * np.random.rand(1)

    # call stochastic ZFW
    fpred, f, w, mean, t, loss, f_values = stochasticZFW(F, d, w0, method = "IRDSA", r=1, T=1000, eps=1e-8)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))

    np.save('../Data/results/function_SZFW_IRDSA_lasso.npy',f_values)
    np.save('../Data/results/loss_SZFW_IRDSA_lasso.npy',loss)
