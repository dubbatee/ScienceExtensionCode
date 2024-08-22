import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import statistics
from scipy import stats

from dsdatasetgeneration import cleanlmcds, cleansmcds 

### Distances to the LMC and SMC ### Used to Generate the P-L Relations

DistSMC = 62440
DistLMC = 49590

### Creating new dataframes for the Fundamental Mode of Pulsation and the First Overtone of Pulsation ###

#SMC Dataframes

Fundamentalsmcds = cleansmcds.drop(cleansmcds[cleansmcds['P2'] != (-99.99)].index)

FirstOvertonesmcds = cleansmcds.drop(cleansmcds[cleansmcds['mode'] == ('singlemode')].index)

#LMC Dataframes

Fundamentallmcds = cleanlmcds.drop(cleanlmcds[cleanlmcds['P2'] != (-99.99)].index)

FirstOvertonelmcds = cleanlmcds.drop(cleanlmcds[cleanlmcds['mode'] == ('singlemode')].index)

print(len(Fundamentalsmcds))
print(len(FirstOvertonesmcds))

print(len(Fundamentallmcds))
print(len(FirstOvertonelmcds))

### Generating P-L Relations, Absolute Magnitude (M) vs log Period (logP) ###

# Functions required:

#1. Function that returns the log10 of the Period of pulsation 

def logperiod(P):
    return np.log10(P)

#2. Function that calculates the Absolute Magnitude given the assumed constant distance to the specific Magellanic cloud

def Msmc(I):
    return(I - 5*np.log10(DistSMC/10)) # For SMC

def Mlmc(I):
    return(I - 5*np.log10(DistLMC/10)) # For LMC

## Generating P-L relation graphs and acquiring the coefficients of the relationships ##

# Fundamental Mode P-L Relations

#SMC P-L Relationship 
y = ((Msmc(Fundamentalsmcds['I'])).tolist())
x = ((logperiod(Fundamentalsmcds['P1'])).tolist())

plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * a + b for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(-1, -2, 'M = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('Delta Scuti SMC Fundamental Mode: Absolute Magnitude vs logP')
plt.ylim(-4, 4)
plt.xlim(-1.6,-0.4)
plt.show()



#LMC P-L Relationship
w = ((Mlmc(Fundamentallmcds['I'])).tolist())
q = ((logperiod(Fundamentallmcds['P1'])).tolist())

plt.scatter(q, w, marker=".")
c, d = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolue I-band Magnitude')
LOBF = [i * c + d for i in q]
plt.plot(q, LOBF, color = "red")
plt.text(-0.8, -2, 'M = ' + format(c.round(3)) + 'logP ' + "+ " + format(d.round(3)))
plt.title('Delta Scuti LMC Fundamental Mode: Absolute Magnitude vs logP')
plt.ylim(-4, 4)
plt.xlim(-1.6,-0.4)
plt.show()



# First Overtone P-L Relations

#SMC P-L Relationship 
y = ((Msmc(FirstOvertonesmcds['I'])).tolist())
x = ((logperiod(FirstOvertonesmcds['P2'])).tolist())

plt.scatter(x, y, marker=".")
m, n = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * m + n for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(-0.9, -2, 'M = ' + format(m.round(3)) + 'logP ' + "+ " + format(n.round(3)))
plt.title('Delta Scuti SMC First Overtone: Absolute Magnitude vs logP')
plt.ylim(-4, 4)
plt.xlim(-1.6,-0.4)
plt.show()



#LMC P-L Relationship
w = ((Mlmc(FirstOvertonelmcds['I'])).tolist())
q = ((logperiod(FirstOvertonelmcds['P2'])).tolist())

plt.scatter(q, w, marker=".")
k, l = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * k + l for i in q]
plt.plot(q, LOBF, color = "red")
plt.text(-0.9, -2, 'M = ' + format(k.round(3)) + 'logP ' + "+ " + format(l.round(3)))
plt.title('Delta Scuti LMC First Overtone: Absolute Magnitude vs logP')
plt.ylim(-4, 4)
plt.xlim(-1.6,-0.4)
plt.show()


