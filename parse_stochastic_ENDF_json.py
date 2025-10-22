import pandas as pd
import json
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

path_to_json = 'json_data/235u_fission.json'
path_out = '235u'

# Load your JSON file
with open(path_to_json, 'r') as f:
    data = json.load(f)

# Extract points from each dataset
datasets = {}
polynomials = {}

plt.figure()

# we're going to have these points be evenly-spaced in log-space according to MeV
MeV_min = 5e-5
MeV_max = 0.00237537

CS_min = 0.5
CS_max = 5.e2

# this means we will have 50 pts per decade (each power of 10)
point_log_spacing = np.log(10.)/50. 

i_color=0
colors = ['black', 'blue', 'red', 'green', 'violet']


# take the json data, and put it into a dictionary, labeled with the dataset name & the numpy arrays of cross sections 
for func in data['funcs']:
    
    dataset_name = func['fName']
    points = func.get('pts')
    
    if not points: 
        continue

    # Get the dataframe
    df = pd.DataFrame(points)

    # Get numpy arrays
    pts_MeV = df['x'].to_numpy()
    pts_CS  = df['y'].to_numpy()

    # fit a 1d-line to this polynomial 
    linear_fit = Polynomial.fit( np.log(pts_MeV), np.log(pts_CS), deg=1)

    m = linear_fit.coef[1]
    b = linear_fit.coef[0]

    pts_CS_fit = np.exp(-1.6) * np.exp( np.log(pts_MeV) * -0.5)

    print(f'points_CS_fit[0] = {pts_CS_fit[0]:.4e}')

    datasets[dataset_name] = np.stack([pts_MeV, pts_CS], axis=-1)

    print(f'parsed dataset for "{dataset_name}", which has {len(pts_MeV)} points.')

    """ ax_scat.scatter(pts_MeV, pts_CS, color='blue', alpha=0.05, marker=".")
    ax_scat.plot(pts_MeV, pts_CS_fit, color='red', alpha=1.)
 """
    plt.hist(np.log10(pts_CS) - np.log10(pts_CS_fit), bins=100)



    
""" plt.xlim(MeV_min, MeV_max)
plt.xscale('log')

plt.ylim(CS_min, CS_max)
plt.yscale('log')
 """# plt.legend()
""" 
ax_scat.set_xlabel('Incident Neutron Energy (MeV)')
ax_scat.set_ylabel('Fission Cross Section (b)')

ax_hist.set_xlabel('Cross section - Model (b)') """

plt.set_xlabel('')

#plot_name = '{}_cross_sections.png'.format(path_out)
#plt.savefig(plot_name)

plt.show()