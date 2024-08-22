import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
from scipy import stats

from cephdatasetgeneration import cleanlmcceph, cleansmcceph 

### Distances to the LMC and SMC ### Used to Generate the P-L Relations

DistSMC = 62440
DistLMC = 49590

### Creating new dataframes for the Fundamental Mode of Pulsation and the First Overtone of Pulsation ###

#SMC Dataframes

Fundamentalsmcceph = cleansmcceph.drop(cleansmcceph[cleansmcceph['mode'] != ('F')].index)

FirstOvertonesmcceph = cleansmcceph.drop(cleansmcceph[cleansmcceph['mode'] == ('F')].index)

#LMC Dataframes

Fundamentallmcceph = cleanlmcceph.drop(cleanlmcceph[cleanlmcceph['mode'] != ('F')].index)

FirstOvertonelmcceph = cleanlmcceph.drop(cleanlmcceph[cleanlmcceph['mode'] == ('F')].index)

#print((Fundamentalsmcceph))
#print(len(FirstOvertonesmcceph))

#print(len(Fundamentallmcceph))
#print(len(FirstOvertonelmcceph))

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
y = ((Msmc(Fundamentalsmcceph['I'])).tolist())
x = ((logperiod(Fundamentalsmcceph['P1'])).tolist())
#Delta Scuti SMC Fundamental Mode: Absolute Magnitude vs logP
plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * a + b for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(0.4, -6, 'M = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('Cepheid SMC Fundamental Mode: Absolute Magnitude vs logP')
plt.show()



#LMC P-L Relationship
w = ((Mlmc(Fundamentallmcceph['I'])).tolist())
q = ((logperiod(Fundamentallmcceph['P1'])).tolist())

plt.scatter(q, w, marker=".")
c, d = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolue I-band Magnitude')
LOBF = [i * c + d for i in q]
plt.plot(q, LOBF, color = "red")
plt.text(0.4, -5, 'M = ' + format(c.round(3)) + 'logP ' + "+ " + format(d.round(3)))
plt.title('Cepheid LMC Fundamental Mode: Absolute Magnitude vs logP')
plt.show()



# First Overtone P-L Relations

#SMC P-L Relationship 
y = ((Msmc(FirstOvertonesmcceph['I'])).tolist())
x = ((logperiod(FirstOvertonesmcceph['P1'])).tolist())

plt.scatter(x, y, marker=".")
m, n = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolue I-band Magnitude')
LOBF = [i * m + n for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(-0.2, -4, 'M = ' + format(m.round(3)) + 'logP ' + "+ " + format(n.round(3)))
plt.title('Cepheid SMC First Overtone: Absolute Magnitude vs logP')
plt.show()



#LMC P-L Relationship
w = ((Mlmc(FirstOvertonelmcceph['I'])).tolist())
q = ((logperiod(FirstOvertonelmcceph['P1'])).tolist())

plt.scatter(q, w, marker=".")
k, l = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolue I-band Magnitude')
LOBF = [i * k + l for i in q]
plt.plot(q, LOBF, color = "red")
plt.text(-0.2, -5, 'M = ' + format(k.round(3)) + 'logP ' + "+ " + format(l.round(3)))
plt.title('Cepheid LMC First Overtone: Absolute Magnitude vs logP')
plt.show()




