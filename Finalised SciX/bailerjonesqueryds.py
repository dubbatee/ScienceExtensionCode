import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.coordinates as coord
#from intragen import testingset
from intragen import raw1blgds, raw1diskds
import statistics
import numpy as np

### Querying the Bailer Jones Catalogue for distances ###

print('Beginning DS Query')

def dsquery(dataframe):

    # Create a new empty column for distances
    dataframe['Distance'] = None
    dataframe['SourceID'] = None

    # Iterate through each row in the DataFrame
    for index, row in dataframe.iterrows():
        ra = row['Ra']
        decl = row['Decl']

        try:
            # Connect to Vizier catalog
            v = Vizier(catalog="I/352", columns=['*', 'RA_ICRS', 'DE_ICRS'])

            # Query for nearby objects within 0.0001 degrees radius
            result = v.query_region(
                coord.SkyCoord(ra=ra, dec=decl, unit=(u.deg, u.deg), frame='icrs'),
                radius=0.0001 * u.deg
            )

            # Check for results
            if len(result) > 0:
                # If results found, extract distance from the first object
                distance = result[0]['rgeo'][0]
                dataframe.at[index, 'Distance'] = distance
                id = result[0]['Source'][0]
                dataframe.at[index, 'SourceID']= id
                
            else:
                # If no results, set distance to 0
                dataframe.at[index, 'Distance'] = 0
                dataframe.at[index, 'SourceID'] = 0

        except Exception as e:
            # Handle potential exceptions during Vizier query
            print(f"Error querying Vizier for row {index}: {e}")
            dataframe.at[index, 'Distance'] = 0  # Set distance to 0 for error cases
            dataframe.at[index, 'SourceID'] = 0

    return dataframe

### BLG Delta Scuti QUERY / 1st ROUND CLEANSING

rawblgdsdistances = dsquery(raw1blgds.copy())  # Avoid modifying original DataFrame
#print(rawblgdsdistances)

cleanblgds = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2', 'Distance', 'SourceID']) #creating a new dataframe for the cleansed data
for index, row in rawblgdsdistances.iterrows(): #looping through raw data dataframe
    if np.float64(row['Distance']) != 0:
        cleanblgds.loc[len(cleanblgds.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['P2'], row['Distance'], row['SourceID']] # adds the star at the end of the dataframe (last index)
    else:
        continue

#print(cleanblgds)
print('stars removed from DS BLG', len(cleanblgds) - len(raw1blgds))
print('CHECK 1')


### DISK Delta Scuti QUERY / 1st ROUND CLEANSING ### 
#print(raw1diskds)

rawdiskdsdistances = dsquery(raw1diskds.copy())  # Avoid modifying original DataFrame
#print(rawdiskdsdistances)

cleandiskds = pd.DataFrame(columns = ['ID', 'mode', 'Ra', 'Decl', 'I', 'V', 'V-I', 'P1', 'P2', 'Distance', 'SourceID']) #creating a new dataframe for the cleansed data
for index, row in rawdiskdsdistances.iterrows(): #looping through raw data dataframe
    if np.float64(row['Distance']) != 0:
        cleandiskds.loc[len(cleandiskds.index)] = [row['ID'], row['mode'], row['Ra'], row['Decl'], row['I'], row['V'], row['V-I'], row['P1'], row['P2'], row['Distance'], row['SourceID']] # adds the star at the end of the dataframe (last index)
    else:
        continue

#print(cleandiskds)
print('stars removed from DS DISK', len(cleandiskds) - len(raw1diskds))
print('CHECK 2')

print('Beginning Ceph Querying shortly')

### Functions to generate result data ###

def logperiod(P):
    return np.log10(P)


print('Mean m of BLG Delta Scuti', statistics.mean(cleanblgds['I']))
print('Mean LogP of BLG Delta Scuti', statistics.mean(logperiod(cleanblgds['P1'])))

print('Mean m of DISK Delta Scuti', statistics.mean((cleandiskds['I'])))
print('Mean LogP of DISK Delta Scuti', statistics.mean(logperiod(cleandiskds['P1'])))


