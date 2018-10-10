import numpy as np
import pandas as pd

inPath = '/Users/i6o/Documents/COSC 528/Project 2/'
df = pd.read_excel(inPath+'UTK-peers.xlsx')
first_row_with_all_NaN = df.shape[0]-df[df.isnull().all(axis=1) == True].shape[0]
df = df.loc[0:first_row_with_all_NaN-1]
df = df.drop(['IPEDS#','HBC','Carm R1','Med School Res $'],axis=1)
print df.isnull().any()
