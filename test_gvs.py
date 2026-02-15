#!/usr/bin/env python
# coding: utf-8

#=============================================================================
#import pandas as pd
import numpy as np

#=============================================================================
from scipy.special import mathieu_cem, mathieu_sem
print("mathieu_cem(1, -10, -180) = ")
mathieu_cem(1, -10, -180)[0]
print("----------------------------------------")
    
#==============================================================================
import glob
import os
import csv
import re

#file_pattern = 'python_generated_gvs/mathieu_ce_gvs_q*'

file_pattern = 'matlab_generated_gvs/mathieu_ce_gvs_q*'
#  % The Matlab GV files are layed out like this:
#  % q
#  % v1, m=0, m=1, m=2, ... m=Ne
#  % v2, m=0, m=1, m=2, ... m=Ne    
#  % v3, m=0, m=1, m=2, ... m=Ne        

matching_files = glob.glob(file_pattern)
row_headers = np.linspace(-np.pi, np.pi, 2500)
row_index = range(2500)

for file in matching_files:
    print("Opening file ")
    print(file)
    data = np.loadtxt(file, delimiter=',', skiprows=1)

    q = float(re.findall(r'-?\d+(?:\.\d+)?', file)[0])
    print(q)
    print("Testing q = %f" % q)

    diff_npass = []
    sp_npass = []
    gv_npass = []
    m_npass = []
    v_npass = []

    for v, v_ind in zip(row_headers, row_index):
        for m in range(1, 35): # how do i get the corresponding value?
            sp = mathieu_cem(m, q, v)[0]
            gv = data[v_ind, m + 1]
            print("sp = %f, gv = %f, diff = %e" % (sp, gv, sp-gv))
            if not np.isclose(sp, gv, atol=2e-6): # why is it has to be m+1 
                print("Fail!")
                #diff = np.abs(sp - gv)
                #diff_npass.append(diff)
                #sp_npass.append(sp)
                #gv_npass.append(gv)
                #m_npass.append(m)
                #v_npass.append(v)
                # if np.isclose(q, 1.0):
                #     print(f'mathieu_cem when q:{q}, v:{v}, m:{m} does not pass the test, scipy value:{sp}, gv: {gv}, diff: {diff} ')

#    results_df = pd.DataFrame({'m': m_npass, 'v': v_npass, 'gv': gv_npass, 'sp': sp_npass, 'diff': diff_npass})
#    results_df.to_csv(f"python_not_passed_ce/not_passed_ce_q{q}.csv", index=False)

    print("----------------------------------------")

#=============================================================================
print(len('ScipyMathieuTesting/python_generated_gvs/python_mathieu_se_gvs_q'))

#============================================================================
file_pattern = 'python_generated_gvs/python_mathieu_se_gvs_q*'
matching_files = glob.glob(file_pattern)
row_headers = np.linspace(-np.pi, np.pi, 2500)
row_index = range(2500)

for file in matching_files:
    print("Opening file ")
    print(file)
    data = np.loadtxt(file, delimiter=',', skiprows=1)
    q = float(file[64:-4])
    print("Testing q = %f" % q)


    diff_npass = []
    sp_npass = []
    gv_npass = []
    m_npass = []
    v_npass = []

    for v, v_ind in zip(row_headers, row_index):
        for m in range(1, 36): # how do i get the corresponding value?
            sp = mathieu_sem(m, q, v)[0]
            gv = data[v_ind, m ]
            if not np.isclose(sp, gv, atol=2e-6): # why is it has to be m+1
                print("Fail!")                
                #diff = np.abs(sp - gv)
                #diff_npass.append(diff)
                #sp_npass.append(sp)
                #gv_npass.append(gv)
                #m_npass.append(m)
                #v_npass.append(v)
                # print(f'mathieu_sem when q:{q}, v:{v}, m:{m} does not pass the test, scipy value:{sp}, gv: {gv}, diff: {diff} ')

#    results_df = pd.DataFrame({'m': m_npass, 'v': v_npass, 'gv': gv_npass, 'sp': sp_npass, 'diff': diff_npass})
#    results_df.to_csv(f"python_not_passed_se/not_passed_se_q{q}.csv", index=False)
    print("----------------------------------------")

    
print('Finished running test')


# print(mathieu_sem(1, -0.001, -180))

