import pandas as pd
import json
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

path_to_json = 'json_data/235-test.json'
path_out = '235u'

# Load your JSON file
with open(path_to_json, 'r') as f:
    data = json.load(f)

# Extract points from each dataset
datasets = {}
polynomials = {}

plt.figure()

# we're going to have these points be evenly-spaced in log-space according to MeV
MeV_min = 0.5
MeV_max = 50.

CS_min = 1.e-1
CS_max = 1.e3 

# this means we will have 50 pts per decade (each power of 10)
point_log_spacing = np.log(10.)/200. 

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

    datasets[dataset_name] = np.stack([pts_MeV, pts_CS], axis=-1)

    print(f'parsed dataset for "{dataset_name}", which has {len(pts_MeV)} points.')

    plt.plot( pts_MeV, pts_CS, color=colors[i_color], linestyle=":")
    i_color += 1

print('~~~~~~~~~~~ made evenly spaced arrays.')

# now, we construct the points, and output them in arrays
evenly_spaced_arrays = {}
i_color =0 
for data_name in datasets.keys():

    rawdata = datasets[data_name]
    
    MeV_lowest  = rawdata[0,0]
    MeV_highest = rawdata[-1,0]

    n_points = int(np.ceil(np.log(MeV_highest / MeV_lowest) / point_log_spacing))

    pts_output_MeV = np.logspace( np.log10(MeV_lowest), np.log10(MeV_highest), n_points)
    pts_output_CS  = np.zeros(n_points)
    
    print(f'shape of raw data "{data_name}": {np.shape(rawdata)}')
    print(f'energy domain [min, max] => [{MeV_lowest}, {MeV_highest}], n.points = {n_points}')

    for i in range(0, n_points):
        MeV = pts_output_MeV[i]

        # if this datapoint is outside of the range this json data covers, then set it to the last datapoint. 
        if MeV < MeV_lowest: 
            pts_output_CS[i] = rawdata[0,1] 
            continue 

        if MeV > MeV_highest: 
            pts_output_CS[i] = rawdata[-1,1]
            continue 

        for j in range(1, np.shape(rawdata)[0]-1):
            json_MeV = rawdata[j,0]
            if MeV < json_MeV: 
                
                json_MeV_last = rawdata[j-1,0]

                json_CS      = rawdata[j,1]
                json_CS_last = rawdata[j-1,1]

                # do a linear interpolation
                CS = json_CS_last + (MeV - json_MeV_last)* ((json_CS - json_CS_last)/(json_MeV - json_MeV_last))

                pts_output_CS[i] = CS
                break


    evenly_spaced_arrays[data_name] = np.stack([pts_output_MeV, pts_output_CS], axis=-1)

    plt.scatter( pts_output_MeV, pts_output_CS, color=colors[i_color], marker='+', label=data_name)
    i_color += 1 


print('~~~~~~~~~~~ made evenly spaced arrays.')

with open('cpp_array_data/' + path_out + '.txt', 'w') as f:

    f.write(f'// This data is created by "parse_ENDF_json.py", using the dataset "{path_to_json}".\n')
    f.write(f'// The points are evenly, GEOMETRICALLY spaced incident neutron energy, with the first point correpsonding to\n')
    f.write(f'// Mev_min, and the last point corresponding to MeV_max.\n')
    f.write(f'// Each point in the vector is the cross section, in barns.\n\n\n\n')

    for data_name in evenly_spaced_arrays.keys():

        str_info = '//species: "{}" reaction: "{}" from file: "{}"\n'.format(path_out, data_name, path_to_json)

        f.write(str_info)

        rawdata = evenly_spaced_arrays[data_name]

        str_start = 'fEDCS["{}_{}"]'.format(path_out, data_name) + ' = new EnergyDependentCS(' + '{:.8e}, {:.8e}, '.format(rawdata[0,0], rawdata[-1,0]) + '{\n'
        f.write(str_start)

        # we need to do the last line separately, as we don't want to add a comment after the last element. 
        for i in range(0, np.shape(rawdata)[0]-1): 
            str_line = '   {:+.8e},\n'.format(rawdata[i,1])
            f.write(str_line)

        str_line = '   {:+.8e}\n'.format(rawdata[-1,1])
        f.write(str_line)

        f.write('});\n\n\n')

plt.xlim(MeV_min, MeV_max)
plt.xscale('log')

plt.ylim(CS_min, CS_max)
plt.yscale('log')
plt.legend()

plot_name = '{}_cross_sections.png'.format(path_out)
plt.savefig(plot_name)

plt.show()