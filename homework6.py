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
    solar_pressure = 4.57 * math.pow(10, -6) # Given on slide 22
    omega_earth = 72.921151467 * math.pow(10, -6) # rad/sec
    omega_earth_vector = [0, 0, omega_earth]
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
    sun_acceleration = keHelperFunctions.compute_third_body_acceleration(sun_mu, ke1.r_vector, sun_vector)
    print(f'Acceleration of the sun: {sun_acceleration}')


    # Compute 3rd-body acceleration for moon
    moon_acceleration = keHelperFunctions.compute_third_body_acceleration(moon_mu, ke1.r_vector, moon_vector)
    print(f'Acceleration of the moon: {moon_acceleration}')


    # Compute drag acceleration with for F10=100
    # Remember: model uses F10scaled = F10/100
    velocity_relative_to_atmosphere = keHelperFunctions.compute_velocity_relative_to_atmosphere(ke1.r_vector, ke1.r_dot_vector, omega_earth)
    lat, lon, alt = keHelperFunctions.compute_lat_lon_alt(ke1.r_vector)
    print(f'LAT: {lat}')
    print(f'LON: {lon}')
    print(f'ALT: {alt*.001}')

    alt_km = alt*.001 # Convert altitude from meters to km

    atmospheric_density = keHelperFunctions.compute_atmospheric_density(utc_time, alt_km, ke1.r_vector, sun_vector)

    atmospheric_drag = keHelperFunctions.compute_atmospheric_drag(drag_coefficient, drag_area, vehicle_mass, atmospheric_density, ke1.r_vector, ke1.r_dot_vector, omega_earth_vector)

    print(f'Atmospheric Drag Acceleration: {atmospheric_drag}')


    # Compute Solar Radiation Pressure
    solar_radiation_acceleration = keHelperFunctions.compute_solar_radiation(solar_pressure_area,
                                                                             radiation_coefficient, 
                                                                             vehicle_mass, 
                                                                             1, # Assumed that there isn't an ecclipse
                                                                             solar_pressure,
                                                                             sun_vector)

    print(f'Solar Radiation Acceleration: {solar_radiation_acceleration}')

    # Compute total acceleration on SV
    accelerants = []
    accelerants.append(sun_acceleration)
    accelerants.append(moon_acceleration)
    accelerants.append(solar_radiation_acceleration)
    accelerants.append(atmospheric_drag)

    total_accelerations = keHelperFunctions.compute_total_acceleration_sv(accelerants)
    print(f'Total Acceleration on SV: {total_accelerations}')



if __name__ == '__main__':
    main()
