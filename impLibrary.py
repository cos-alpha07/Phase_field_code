# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 12:12:31 2025

@author: Dell
"""

import os, sys, time
from datetime import datetime

from math import sqrt, pi
from sympy import diff, symbols, evalf
from numpy import array, zeros, zeros_like, ix_, dot, asarray, empty, amax, sum, arange, append, linspace, nan
from numpy import newaxis, maximum, cross
from numpy.linalg import eigvals
from scipy.linalg import det, inv, norm, solve
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

import pickle
from pandas import DataFrame, Series, read_csv

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.tri as tri
import scienceplots
plt.style.use(['science', 'ieee', 'high-vis', 'no-latex'])