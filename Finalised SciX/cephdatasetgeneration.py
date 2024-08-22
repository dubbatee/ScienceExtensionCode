import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
from scipy import stats


### imported metadata from the OGLE IV Catalogue of Classical Cepheid Variables into pandaS dataframes ###

#Small Megellanic Cloud Classical Cepheid Variable Dataframe
rawsmcceph = pd.read_csv ('smccephdata.csv')
rawsmcceph.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1']

#Large Megellanic Cloud Classical Cepheid Variable Dataframe
rawlmcceph = pd.read_csv ('lmccephdata.csv')
rawlmcceph.columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1']


### data cleaning/pre-processing ###

#apparent magnitude thresholding - based on OGLE IV telescope saturation and sensitivity limits

print('Raw CEPH SMC size', len(rawsmcceph))
print('Raw CEPH LMC size', len(rawlmcceph))

rawsmcceph = rawsmcceph[rawsmcceph['I'] > 13]
rawsmcceph = rawsmcceph[rawsmcceph['I'] < 21.5]

rawlmcceph = rawlmcceph[rawlmcceph['I'] > 13]
rawlmcceph = rawlmcceph[rawlmcceph['I'] < 21.5]

print('Stars in CEPH SMC after sensitivity thresh-holding', len(rawsmcceph))
print('Stars in CEPH LMC after sensitivity thresh-holding', len(rawlmcceph))

#removal of Milkyway Halo Classical Cepheids - Equal distance approximation cannot be used with Milkyway Classical Cepheids - creating new dataframes

# coefficients of P-L relationship prior to cleansing - required to determine whether to remove star - the loops below take the true apparent magnitude of a star, and compare it with the apparent magnitude the star should have based on the Line of best fit
# if the star's apparent magnitude is lower than the line of best fit apparent magnitude by more than 1.5, the code keeps this star in the raw dataframe (dfISMC / dfILMC)
# if the star's apparent magnitude is greater than this thresholded value (1.5 lower than the LOBF apparent magnitude), the loop adds this star with its ID, apparent magnitude and period to a "cleansed" dataframe (cleanlmcds / cleansmcds)

#1. First step is plotting the line of best fit for the stars: Apparent Magnitude (m) vs Log(10) Period and acquiring the coefficients of the gradient and y intercept. The period used is the fundamental mode period.

# Function that calculates the Log10 of the Period of Pulsation

def logperiod(P):
    return(np.log10(P))

# Graphing the m vs logP relation and acquiring the coefficients of the relation: 

# a is the gradient, b is the y intercept for the SMC relationship

y = (((rawsmcceph['I'])).tolist())
x = ((logperiod(rawsmcceph['P1'])).tolist())

plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Apparent I-band Magnitude')
LOBF = [i * a + b for i in x]
#LOBF2 = [i * a + (b - 1.5) for i in x]
plt.plot(x, LOBF, color = "red")
#plt.plot(x, LOBF2, color = 'blue')
plt.text(-0.25, 14, 'm = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('Raw Cepheid SMC: Apparent Magnitude vs logP')
plt.show()

# c is the gradient, d is the y intercept for the LMC relationship

w = (((rawlmcceph['I'])).tolist())
q = ((logperiod(rawlmcceph['P1'])).tolist())

plt.scatter(q, w, marker=".")
c, d = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Apparent I-band Magnitude')
LOBF = [i * c + d for i in q]
#LOBF2 = [i * c + (d - 1.5) for i in q]
plt.plot(q, LOBF, color = "red")
#plt.plot(q, LOBF2, color = 'blue')
plt.text(-0.25, 13.5, 'm = ' + format(c.round(3)) + 'logP ' + "+ " + format(d.round(3)))
plt.title('Raw Cepheid LMC: Apparent Magnitude vs logP')
plt.show()

#2. Visualisation of the new dataset's m vs logP relationship

cleansmcceph = rawsmcceph.copy()
cleanlmcceph = rawlmcceph.copy()

y = (((cleansmcceph['I'])).tolist())
x = ((logperiod(cleansmcceph['P1'])).tolist())

plt.scatter(x, y, marker='.')
m, n = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Apparent I-band Magnitude')
LOBF = [i * m + n for i in x]
#LOBF2 = [i * a + (b - 1.5) for i in x]
plt.plot(x, LOBF, color = "red")
#plt.plot(x, LOBF2, color = 'blue')
plt.text(-0.25, 14, 'm = ' + format(m.round(3)) + 'logP ' + "+ " + format(n.round(3)))
plt.title('Cleansed Cepheid SMC: Apparent Magnitude vs logP')
plt.show()

w = (((cleanlmcceph['I'])).tolist())
q = ((logperiod(cleanlmcceph['P1'])).tolist())

plt.scatter(q, w, marker='.')
k, l = np.polyfit(q, w, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Apparent I-band Magnitude')
LOBF = [i * k + l for i in q]
#LOBF2 = [i * c + (d - 1.5) for i in q]
plt.plot(q, LOBF, color = "red")
#plt.plot(q, LOBF2, color = 'blue')
plt.text(-0.25, 13.5, 'm = ' + format(k.round(3)) + 'logP ' + "+ " + format(l.round(3)))
plt.title('Cleansed Cepheid LMC: Apparent Magnitude vs logP')
plt.show()

### Functions to generate result data ###

def Msmc(I):
    return(I - 5*np.log10(62440/10)) # For SMC

def Mlmc(I):
    return(I - 5*np.log10(49590/10)) # For LMC

def logperiod(P):
    return np.log10(P)


print('Mean m of SMC Cepheid', statistics.mean(cleansmcceph['I']))
print('Mean LogP of SMC Cepheid', statistics.mean(logperiod(cleansmcceph['P1'])))
print('Mean M of SMC Cepheid', statistics.mean(Msmc(cleansmcceph['P1'])))
print('Mean m of LMC Cepheid', statistics.mean((cleanlmcceph['I'])))
print('Mean LogP of LMC Cepheid', statistics.mean(logperiod(cleanlmcceph['P1'])))
print('Mean M of LMC Cepheid', statistics.mean(Msmc(cleanlmcceph['P1'])))

#1. Function that returns the log10 of the Period of pulsation 

def logperiod(P):
    return np.log10(P)

#2. Function that calculates the Absolute Magnitude given the assumed constant distance to the specific Magellanic cloud

def Msmc(I):
    return(I - 5*np.log10(62440/10)) # For SMC

def Mlmc(I):
    return(I - 5*np.log10(49590/10)) # For LMC

## Generating P-L relation graphs and acquiring the coefficients of the relationships ##

# Whole Dataset P-L relations

#SMC P-L Relationship 
y = ((Msmc(cleansmcceph['I'])).tolist())
x = ((logperiod(cleansmcceph['P1'])).tolist())
#Delta Scuti SMC Fundamental Mode: Absolute Magnitude vs logP
plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * a + b for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(0.4, -6, 'M = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('Cepheid SMC: Absolute Magnitude vs logP')
plt.show()

y = ((Msmc(cleanlmcceph['I'])).tolist())
x = ((logperiod(cleanlmcceph['P1'])).tolist())
#Delta Scuti SMC Fundamental Mode: Absolute Magnitude vs logP
plt.scatter(x, y, marker=".")
a, b = np.polyfit(x, y, 1)
plt.xlabel('Log10 of Period(days)')
plt.ylabel('Absolute I-band Magnitude')
LOBF = [i * a + b for i in x]
plt.plot(x, LOBF, color = "red")
plt.text(0.4, -6, 'M = ' + format(a.round(3)) + 'logP ' + "+ " + format(b.round(3)))
plt.title('Cepheid LMC: Absolute Magnitude vs logP')
plt.show()