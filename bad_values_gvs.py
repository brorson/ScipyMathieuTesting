import pandas as pd
import glob
import numpy as np

file_pattern = 'python_not_passed_ce/not_passed_ce_q*'
matching_files = glob.glob(file_pattern)


col = 'diff'

for file in matching_files:

    df = pd.read_csv(file)
   
    q = float(file[36:-4])
    print(f'q = {q}')
    threshold_value = 0.01

    if pd.api.types.is_numeric_dtype(df[col]): # Check if the column is numeric
        df[f'{col}_Flag'] = df[col] > threshold_value

    big_diffs_rows = df.loc[df[col] == True]
    print(big_diffs_rows)


file_pattern = 'python_not_passed_se/not_passed_se_q*'
matching_files = glob.glob(file_pattern)


col = 'diff'

for file in matching_files:

    df = pd.read_csv(file)
   
    q = float(file[36:-4])
    print(f'q = {q}')
    threshold_value = 0.01

    if pd.api.types.is_numeric_dtype(df[col]): # Check if the column is numeric
        df[f'{col}_Flag'] = df[col] > threshold_value

    big_diffs_rows = df.loc[df[col] == True]
    print(big_diffs_rows)