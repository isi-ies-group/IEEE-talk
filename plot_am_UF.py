# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 00:44:03 2020

@author: Ruben
"""
import matplotlib.pyplot as plt
import pandas as pd

import cpvlib

am_range = pd.Series(range(1, 10, 1))
plt.plot(am_range, cpvlib.get_simple_util_factor(am_range, 1.7, 0.1, -0.1))
plt.xlabel('AM')
plt.ylabel('UF(AM)')