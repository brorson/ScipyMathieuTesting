from scipy.special import mathieu_cem, mathieu_sem
import numpy as np
import glob
import os
import csv
import pandas as pd


file_pattern = 'ScipyMathieuTesting/python_mathieu_ce_gvs_m*'
matching_files = glob.glob(file_pattern)
row_headers = np.linspace(-np.pi, np.pi, 2500)
row_index = range(2500)

m = 6
print(matching_files)
for file in matching_files:
    print(file)
    data = np.loadtxt(file, delimiter=',', skiprows=1)
    
    diff_npass = []
    sp_npass = []
    gv_npass = []
    q_npass = []
    v_npass = []
    q_ind_npass = []

    Q = np.linspace(0.1, 50, 100)
    for v, v_ind in zip(row_headers, row_index):
        for q, q_ind in zip(Q, range(1, 100)): # how do i get the corresponding value?
            sp = mathieu_cem(m, q, v)[0]
            gv = data[v_ind, q_ind]
            if not np.isclose(sp, gv, atol=2e-2): # why is it has to be m+1 
                diff = np.abs(sp - gv)
                diff_npass.append(diff)
                sp_npass.append(sp)
                gv_npass.append(gv)
                q_npass.append(q)
                q_ind_npass.append(q_ind)
                v_npass.append(v)
                # print(f'mathieu_sem when q:{q}, v:{v}, m:{m} does not pass the test, scipy value:{sp}, gv: {gv}, diff: {diff} ')

    results_df = pd.DataFrame({'q': q_npass, 'q_ind_npass': q_ind_npass, 'v': v_npass, 'gv': gv_npass, 'sp': sp_npass, 'diff': diff_npass})
    results_df.to_csv(f"ScipyMathieuTesting/not_passed_ce_m{m}.csv", index=False)


