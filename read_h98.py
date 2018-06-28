#-*-Python-*-
# Created by mcclenaghanj at 31 May 2018  14:51
import pandas as pd

dataset = pd.read_csv('~/tenn_plots/ASP-db1.csv')
X = dataset.iloc[:, 3].values
pinj= pd.to_numeric(dataset['PL'].values,errors ='coerce')
amin = pd.to_numeric(dataset['AMIN'].values,errors ='coerce')
kappa = pd.to_numeric(dataset['KAPPA'].values,errors ='coerce')
delta = pd.to_numeric(dataset['DELTA'].values,errors ='coerce')
zeff = pd.to_numeric(dataset['ZEFF'].values,errors ='coerce')
ip = pd.to_numeric(dataset['IP'].values,errors ='coerce')
nel = pd.to_numeric(dataset['NEL'].values,errors ='coerce')
tauth = pd.to_numeric(dataset['TAUTH'].values,errors ='coerce')
rgeo = pd.to_numeric(dataset['RGEO'].values,errors ='coerce')
bt = pd.to_numeric(dataset['BT'].values,errors ='coerce')
tauc92 = pd.to_numeric(dataset['TAUC92'].values,errors ='coerce')
tauc93 = pd.to_numeric(dataset['TAUC93'].values,errors ='coerce')
meff = pd.to_numeric(dataset['MEFF'].values,errors ='coerce')
#standard = pd.to_numeric(dataset['IAE2004I'].values,errors ='coerce')
standard = pd.to_numeric(dataset['DB3V5'].values,errors ='coerce')
phase = dataset['PHASE'].values
date = dataset['DATE'].values
tok = dataset['TOK'].values
aspect = rgeo/amin

X = np.zeros([8,len(pinj)])
X[0,:] = pinj/(1.0e6)
X[1,:] = 1.0/aspect
X[2,:] = kappa
X[3,:] = ip/(1.0e6)
X[4,:] = nel/(1.0e19)
X[5,:] = rgeo
X[6,:] = bt
X[7,:] = meff
y = tauth*tauc92

Xtmp = []
ytmp = []
weight = []
import math
from collections import Counter

OMFIT['inputs_h98'] = OMFITtree('')

q = 0
for i in range(0,len(pinj)):
    avoidNaNandzeros = (
        X[0,i]!=0. and  X[1,i]!=0. and X[2,i]!=0. and
        X[3,i]!=0. and  X[4,i]!=0. and X[5,i]!=0. and
        X[6,i]!=0.   and X[7,i]!=0.   and
        (not math.isnan(X[0,i])) and (not math.isnan(X[1,i])) and (not math.isnan(X[2,i])) and
        (not math.isnan(X[3,i])) and (not math.isnan(X[4,i])) and (not math.isnan(X[5,i])) and
        (not math.isnan(X[6,i])) and (not math.isnan(X[7,i])))
    if (standard[i]>0):
        if ((phase[i][0:2]=='HS' or phase[i][0:2]=='HG')):
            if True:
                OMFIT['eped_inputs_h98'] = OMFITtree('')

            if False:

                OMFIT['inputs_h98'][q] = copy.deepcopy(OMFIT['inputs'])
                OMFIT['inputs_h98'][q]['Paux'] = pinj[i]/(1.0e6)
                OMFIT['inputs_h98'][q]['aspect'] = aspect[i]
                OMFIT['inputs_h98'][q]['Bt'] = abs(bt[i])
                OMFIT['inputs_h98'][q]['Ip'] = abs(ip[i]/(1.0e6))
                OMFIT['inputs_h98'][q]['R'] = rgeo[i]
                OMFIT['inputs_h98'][q]['kappa'] = kappa[i]
                if math.isnan(zeff[i]):
                    OMFIT['inputs_h98'][q]['zeff']  = 2.0
                else:
                    OMFIT['inputs_h98'][q]['zeff'] = zeff[i]
                OMFIT['inputs_h98'][q]['delta'] = delta[i]
                OMFIT['inputs_h98'][q]['ne_ped'] = nel[i]/(1.0e20)/1.2
                OMFIT['inputs_h98'][q]['tokamak'] = tok[i]
                OMFIT['inputs_h98'][q]['tau_exp'] = y[i]
                q+=1
