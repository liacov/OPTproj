"""
File created 23th June 2020
Authors: Laura Iacovissi, Federico Matteo
"""

import numpy as np
from numba import njit

def e(i, d):
    ei = np.zeros(d)
    ei[i] = 1
    return ei

@njit(parallel=True)
def IRDSA_out(F, w, m, c, z):
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
    F_w = F(w)
    out = np.zeros(z.shape[0])
    for i in range(m):
        out += (F(w + c * z[:,i]) - F_w) / c * z[:,i]
    return out/m

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
    return IRDSA_out(F, w, m, c, z)

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
        loss_eval = np.abs(F(w_pred) - F(w))
        loss.append(loss_eval)
        print(f"Loss evaluation at time {t}:\t{loss_eval:.4f}\n")
        if loss_eval < eps: break # check stopping condition
    return F(w_pred), F(w), w, partial/T, t, loss
