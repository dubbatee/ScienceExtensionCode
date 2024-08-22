import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
from scipy import stats

### imported metadata from the OGLE IV Catalogue of Delta Scuti Variables into pandaS dataframes ###

#Galactic Bulge Delta Scuti Variable Dataframe - Raw
rawblgds = pd.read_csv ('blgdsdatafinal.csv')
rawblgds.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']

#Galactic Disk Delta Scuti Variable Dataframe - Raw
rawdiskds = pd.read_csv ('diskdsdatafinal.csv')
rawdiskds.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']

### imported metadata from the OGLE IV Catalogue of Classical Cepheid Variables into pandaS dataframes ###

#Galactic Bulge Classical Cepheid Variable Dataframe - Raw 
rawblgceph = pd.read_csv ('blgcephdatafinal.csv')
rawblgceph.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1']

#Galactic Disk Classical Cepheid Variable Dataframe - Raw
rawdiskceph = pd.read_csv ('diskcephdatafinal.csv')
rawdiskceph.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1']

### data cleaning/pre-processing ###

# ROUND 1: apparent magnitude thresholding - based on OGLE IV telescope saturation and sensitivity limits

#Delta Scuti Variables

rawblgds = rawblgds[rawblgds['I'] > 13]
raw1blgds = rawblgds[rawblgds['I'] < 21.5]

rawdiskds = rawdiskds[rawdiskds['I'] > 13]
raw1diskds = rawdiskds[rawdiskds['I'] < 21.5]

#Classical Cepheid Variables

rawblgceph = rawblgceph[rawblgceph['I'] > 13]
raw1blgceph = rawblgceph[rawblgceph['I'] < 21.5]

rawdiskceph = rawdiskceph[rawdiskceph['I'] > 13]
raw1diskceph = rawdiskceph[rawdiskceph['I'] < 21.5]

testingset = pd.read_csv ('TESTINGSET2.csv')
testingset.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2']

# ROUND 2: separating the queried results into Fundamental / First Overtone data frames

print('receiving dataframes')

from bailerjonesqueryds import cleanblgds, cleandiskds
from bailerjonesqueryceph import cleanblgceph, cleandiskceph

### Bulge Delta Scuti Dataframe Mode separation ###

Fundamentalblgds = cleanblgds.drop(cleanblgds[cleanblgds['P2'] != (-99.99)].index)

FirstOvertoneblgds = cleanblgds.drop(cleanblgds[cleanblgds['mode'] == ('singlemode')].index)

### Disk Delta Scuti Dataframe Mode separation ###

Fundamentaldiskds = cleandiskds.drop(cleandiskds[cleandiskds['P2'] != (-99.99)].index)

FirstOvertonediskds = cleandiskds.drop(cleandiskds[cleandiskds['mode'] == ('singlemode')].index)

### Bulge Cepheid Dataframe Mode separation ### 

Fundamentalblgceph = cleanblgceph.drop(cleanblgceph[cleanblgceph['mode'] != ('F')].index)

FirstOvertoneblgceph = cleanblgceph.drop(cleanblgceph[cleanblgceph['mode'] == ('F')].index)

### Disk Cepheid Dataframe Mode separation ###

Fundamentaldiskceph = cleandiskceph.drop(cleandiskceph[cleandiskceph['mode'] != ('F')].index)

FirstOvertonediskceph = cleandiskceph.drop(cleandiskceph[cleandiskceph['mode'] == ('F')].index)

### Export dataframes to the Distance dataframes ###

print('exporting dataframes to cephdistances and dsdistances')





