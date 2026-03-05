from keplarianElements import KeplerianElements
import keHelperFunctions

import math
import datetime

import tabulate

def main():
    
    vectors_file = 'vectors.yaml'
    vector_data = keHelperFunctions.read_in_yaml(vectors_file)

    ke1 = KeplerianElements(vector_data['vectors'][f'vector1']['x_pos'],
                               vector_data['vectors'][f'vector1']['y_pos'],
                               vector_data['vectors'][f'vector1']['z_pos'],
                               vector_data['vectors'][f'vector1']['x_velocity'],
                               vector_data['vectors'][f'vector1']['y_velocity'],
                               vector_data['vectors'][f'vector1']['z_velocity'])

    latc = vector_data['vectors']['vector1']['latc']
    latd = vector_data['vectors']['vector1']['latd']
    lon = vector_data['vectors']['vector1']['lon']
    altd = vector_data['vectors']['vector1']['altd']


    # Date provided: 28 Feb 2026 18:22:45 UTC
    utc_time = datetime.datetime(2026, 2, 28, 18, 22, 45)


    # Vector is inertial, but corresponds to ECEF
    # No need to do RNP transformation for this case
    # Pretend that sun/moon vectors are in this frame


    # Compute 3rd-body acceleration for sun


    # Compute 3rd-body acceleration for moon


    # Compute drag acceleration with for F10=100
    # Remember: model uses F10scaled = F10/100


    # Compute Solar Radiation Pressure


    # Compute total acceleration on SV



if __name__ == '__main__':
    main()
