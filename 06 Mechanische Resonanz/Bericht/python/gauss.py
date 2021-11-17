# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:25:42 2021

@author: fritz
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:16:07 2021

@author: fritz
"""

import numpy as np


def get_derivative(func, x, dx, i):
    """
    Parameters
    ----------
    func : Function f(x), x array
    x : array of measurment values
    dx : small value
    i : index, regarding which variable derivative has to be made

    Returns
    -------
    approximate derivative at x
    """

    x0 = x.copy()
    res1 = func(x0)
    x0[i] += dx
    res2 = func(x0)
    return (res2-res1)/dx

def error(func, x, delta_x, dx = 0):
    """
    

    Parameters
    ----------
    func : Function f(x), x array
    x : array of measurment values
    delta_x : array of uncertainties 
    dx : optional. if different dx's want to be used: array of dx's
        default is 1e-6 for all

    Returns
    -------
    error

    """
    x = x.astype(float)
    delta_x = delta_x.astype(float)
    if dx == 0:
        dx = np.ones_like(x)*1e-10
    res = 0
    for i, delta_xi in enumerate(delta_x):
        res += (get_derivative(func, x, dx[i], i)*delta_xi)**2    
    return np.sqrt(res)






    
