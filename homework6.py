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
    drag_coefficient = 2.0
    drag_area = 10 # m^2
    radiation_coefficient = 1.41
    solar_pressure_area = 20 # m^2
    vehicle_mass = 1000 # km
    omega_earth = 72.921151467 * math.pow(10, -6) # rad/sec
    nautical_mile = 1.852 # km
    sun_mu = ke1.mu * 332946.09358859973
    moon_mu = ke1.mu / 81.3005764441083 


    # Vector is inertial, but corresponds to ECEF
    # No need to do RNP transformation for this case
    # Pretend that sun/moon vectors are in this frame
    jd = keHelperFunctions.convert_date_to_jd(utc_time.year, utc_time.month, utc_time.day, utc_time.hour, utc_time.minute, utc_time.second)
    mjd = keHelperFunctions.convert_jd_to_mjd(jd)
    besselian_year = keHelperFunctions.convert_mjd_to_besselian_year(mjd)

    sun_vector = keHelperFunctions.determine_sun_vector_lf(jd)
    moon_vector = keHelperFunctions.determine_moon_vector_lf(jd)


    print(f'JD: {jd}')
    print(f'Sun Vector: {sun_vector}')
    print(f'Moon Vector: {moon_vector}')
    print()

    # Compute 3rd-body acceleration for sun
    sun_acceleration = keHelperFunctions.third_body_acceleration(sun_mu, ke1.r_vector, sun_vector)
    print(f'Acceleration of the sun: {sun_acceleration}')


    # Compute 3rd-body acceleration for moon
    moon_acceleration = keHelperFunctions.third_body_acceleration(moon_mu, ke1.r_vector, moon_vector)
    print(f'Acceleration of the moon: {moon_acceleration}')


    # Compute drag acceleration with for F10=100
    # Remember: model uses F10scaled = F10/100
    

    # Compute Solar Radiation Pressure


    # Compute total acceleration on SV



if __name__ == '__main__':
    main()
