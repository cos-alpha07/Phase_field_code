# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 12:12:31 2025

@author: Dell
"""

import os, sys, time
from datetime import datetime

from math import sqrt

from numpy import array, zeros, zeros_like, ix_, dot, asarray, empty, sum, max, arange, append, linspace, nan
from scipy.linalg import det, inv, norm, solve

import pickle
from pandas import DataFrame, Series

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import scienceplots
plt.style.use(['science', 'ieee', 'high-vis', 'no-latex'])