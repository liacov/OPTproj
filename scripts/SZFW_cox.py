"""
File created 23th June 2020
Authors: Laura Iacovissi, Federico Matteo
"""

import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as pl

from numba import njit
from itertools import product
from scipy.sparse.linalg import norm


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


def IRDSA_out(F, w, m, c, seed):
    """
    Improvised Random Direction stochastic approximation
    for gradient estimation - output function

    INPUT:
    - F: objective function
    - w: current weight
    - m: sample dimension
    - c: costant
    - z: random generated vectors

    """
    np.random.seed(seed)
    z = np.random.normal(0, 1, d)
    F_w = F(w)
    return ((F(w + c * z) - F_w) / c * z)



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
    seeds = np.random.randint(low=1, high=2000, size=m)
    p = mp.Pool(mp.cpu_count())
    out = p.starmap(IRDSA_out, product([F],[w],[m],[c],seeds))
    print(out)
    p.close()
    return np.sum(out)/m


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

                       "IRDSA": {"m": 20,
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
    L = 2/X.shape[0] * np.linalg.norm(X.T @ X)

    # call ZFW with InexactUpdate
    fpred, f, w, mean, t, loss, f_values = stochasticZFW(F, X.shape[1], w0, method = "IRDSA", r=10, T=40, eps=1e-8)
    print('\n\n')
    # print resume
    print(f'OUTPUT:\n\nF(w_pred) = {fpred}\n\nF(w) = {f}\n\nw = {w}\n\naverage w = {mean}\n\nT = {t}')
    # print F(stanting point) VS F(w*)
    print(F(w0), F(w))

    np.save('../Data/results/funtion_SZFW_IRDSA_cox.npy',f_values)
    np.save('../Data/results/loss_SZFW_IRDSA_cox.npy',loss)
