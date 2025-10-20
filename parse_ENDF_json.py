import pandas as pd
import json
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

# Load your JSON file
with open('json_data/235u_cross_sections.json', 'r') as f:
    data = json.load(f)

# Extract points from each dataset
datasets = {}
polynomials = {}

plt.figure()

x_min =  0.025
x_max = 30.000

y_min = 1.e-2 
y_max = 15. 

for func in data['funcs']:
    
    dataset_name = func['fName']
    points = func.get('pts')

    # Convert to DataFrame
    if not points: 
        continue
    
    # Get the dataframe
    df = pd.DataFrame(points)

    # Get numpy arrays
    pts_MeV = df['x'].to_numpy()
    pts_CS  = df['y'].to_numpy()

    if pts_CS.any() <= 1e-6: 
        continue
    
    polynomials[dataset_name] = Polynomial.fit( np.log(pts_MeV), np.log(pts_CS), deg=4 ) 

    plt.plot(pts_MeV, pts_CS, color='black', label=dataset_name)

    
plt.xlim(x_min, x_max)
plt.xscale('log')

plt.ylim(y_min, y_max)
plt.yscale('log')

plt.show()


# Now you have separate DataFrames for each dataset
# Access them like: datasets['dataset'], datasets['dataset2'], etc.

# Example: print all datasets
for name, df in datasets.items():
    print(f"\n{name}:")
    print(df)