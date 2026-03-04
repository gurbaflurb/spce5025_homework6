import yaml
import math

import numpy as np


# Read in a yaml that has all the initial vectors for position and velocity
def read_in_yaml(file_name):
    with open(file_name, 'r') as f:
        data = yaml.load(f.read(), Loader=yaml.SafeLoader)
        return data



def convert_arbitrary_perifocal_to_eci(a, e, inclination, raan, aop, nu) -> tuple:
        '''Converts manually provided perifocal values to ECI coordinates. Uses radians and not degrees. Passing in degrees will mess up the calculation'''

        # Hard coded because a nice solution for this exists (Pulled from slide 74)
        x = [ math.cos(raan)*math.cos(aop) - math.sin(raan)*math.sin(aop)*math.cos(inclination), -math.cos(raan)*math.sin(aop) - math.sin(raan)*math.cos(aop)*math.cos(inclination), math.sin(raan)*math.sin(inclination)]
        y = [ math.sin(raan)*math.cos(aop)+math.cos(raan)*math.sin(aop)*math.cos(inclination), -math.sin(raan)*math.sin(aop)+math.cos(raan)*math.cos(aop)*math.cos(inclination), -math.cos(raan)*math.sin(inclination)]
        z = [ math.sin(inclination)*math.sin(aop), math.sin(inclination)*math.cos(aop), math.cos(inclination)]

        return (x, y, z)



def find_arbitrary_position_and_velocity_vector(a: float, eccentricity: float, nu: float) -> tuple:
    '''Returns a tuple in the form position_vector, velocity_vector'''
    mu = 398600441800000.0 # From WGS84
    
    perifocal = a*(1.0-math.pow(eccentricity, 2))
    radius = perifocal/(1.0 + eccentricity*math.cos(nu))

    return ([radius*math.cos(nu), radius*math.sin(nu), 0.0],
            [-math.sqrt(mu/perifocal)*math.sin(nu), math.sqrt(mu/perifocal)*(eccentricity+math.cos(nu)), 0.0])



def keplarian_rk4(r_vector, r_dot_vector, step, mu):
     '''Takes in the position vector (r_vector), velocity vector (r_dot_vector), r (norm of r_vector), step (Size of step between each round), and mu (Probably from WGS 84).
     Mu is typically provided as meters cubed over seconds squared.'''

     # Commented out sections can help to bug things 

     # Given on slide 22
     r0norm = np.linalg.norm(r_vector)
     rd_a_pt = (-mu/math.pow(r0norm, 3))*r_vector

     r_a = r_vector + (step/2)*r_dot_vector
     rd_a = r_dot_vector + (step/2)*rd_a_pt
     k1 = [r_dot_vector, rd_a_pt]

#      print(f'k1  : {k1[0]}')
#      print(f'k1 a: {k1[1]}')
#      print(f'ra: {r_a}')
#      print(f'va: {rd_a}')
#      print(f'Norm of r: {r0norm}')
#      print()

     # Given on slide 23
     ranorm = np.linalg.norm(r_a)
     rd_b_pt = (-mu/math.pow(ranorm, 3))*r_a

     r_b = r_vector+(step/2)*rd_a
     rd_b = r_dot_vector+(step/2)*rd_b_pt
     k2 = [rd_a, rd_b_pt]

#      print(f'k2  : {k2[0]}')
#      print(f'k2 a: {k2[1]}')
#      print(f'ra: {r_b}')
#      print(f'va: {rd_b}')
#      print(f'Norm of r: {ranorm}')
#      print()

     # Given on slide 24
     rbnorm = np.linalg.norm(r_b)
     rd_c_pt = (-mu/math.pow(rbnorm, 3))*r_b

     r_c = r_vector+(step)*rd_b # These are wrong on slide 24, we multiply by the step instead of step / 2
     rd_c = r_dot_vector+(step)*rd_c_pt # These are wrong on slide 24, we multiply by the step instead of step / 2
     k3 = [rd_b, rd_c_pt]

#      print(f'k3  : {k3[0]}')
#      print(f'k3 a: {k3[1]}')
#      print(f'rc: {r_c}')
#      print(f'vc: {rd_c}')
#      print(f'Norm of r: {rbnorm}')
#      print()


     # Given on slide 25
     k4 = [rd_c, (-mu/math.pow(np.linalg.norm(r_c), 3))*r_c]
#      print(f'k4  : {k4[0]}')
#      print(f'k4 a: {k4[1]}')
#      print()

     step_position_solution = r_vector + step*( ((k1[0])/6) + ((k2[0])/3) + ((k3[0])/3) + ((k4[0])/6) )

     step_velocity_solution = r_dot_vector + step*( ((k1[1])/6) + ((k2[1])/3) + ((k3[1])/3) + ((k4[1])/6) )

     return (step_position_solution, step_velocity_solution)



def keplarian_rk4_oblate_earth(r_vector, r_dot_vector, step):
     '''Takes into account the Earth is Oblate. Takes in the position vector (r_vector), velocity vector (r_dot_vector), r (norm of r_vector), and step (Size of step between each round)
     Mu in this case is not provided by the user but instead defined by the WGS 84 standard for this specific implementation.'''

     # Defined by WGS 84
     earth_radius = 6378137 # This value is in meters
     earth_j2 = 0.00108262998905194
     mu = 398600441800000 # This value is in meters cubed per second squared

     # K1
     # Given on slide 52
     r = math.sqrt(math.pow(r_vector[0], 2) + math.pow(r_vector[1], 2) + math.pow(r_vector[2], 2))

     a_pt_s = r_vector[2]/r
     a_pt_1 = -(mu/math.pow(r, 3))
     a_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(a_pt_s, 2)))*r_vector
     
     a = a_pt_1 * a_pt_2

     r_a = r_vector + (step/2)*r_dot_vector
     rd_a = r_dot_vector + (step/2)*a
     k1 = [r_dot_vector, a]

     # K2
     r = math.sqrt(math.pow(r_a[0], 2) + math.pow(r_a[1], 2) + math.pow(r_a[2], 2))

     b_pt_s = r_a[2]/r
     b_pt_1 = -(mu/math.pow(r, 3))
     b_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(b_pt_s, 2)))*r_a
     
     a = b_pt_1 * b_pt_2

     r_b = r_vector+(step/2)*rd_a
     rd_b = r_dot_vector+(step/2)*a
     k2 = [rd_a, a]

     # K3
     r = math.sqrt(math.pow(r_b[0], 2) + math.pow(r_b[1], 2) + math.pow(r_b[2], 2))

     c_pt_s = r_b[2]/r
     c_pt_1 = -(mu/math.pow(r, 3))
     c_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(c_pt_s, 2)))*r_b
     
     a = c_pt_1 * c_pt_2

     r_c = r_vector+(step)*rd_b # These are wrong on slide 24, we multiply by the step instead of step / 2
     rd_c = r_dot_vector+(step)*a # These are wrong on slide 24, we multiply by the step instead of step / 2
     k3 = [rd_b, a]

     # K4
     # Given on slide 25
     r = math.sqrt(math.pow(r_c[0], 2) + math.pow(r_c[1], 2) + math.pow(r_c[2], 2))

     d_pt_s = r_c[2]/r
     d_pt_1 = -(mu/math.pow(r, 3))
     d_pt_2 = (1+((3*earth_j2)/2)*math.pow(earth_radius/r, 2)*(1-5*math.pow(d_pt_s, 2)))*r_c
     
     a = d_pt_1 * d_pt_2

     k4 = [rd_c, a]

     step_position_solution = r_vector + step*( ((k1[0])/6) + ((k2[0])/3) + ((k3[0])/3) + ((k4[0])/6) )

     step_velocity_solution = r_dot_vector + step*( ((k1[1])/6) + ((k2[1])/3) + ((k3[1])/3) + ((k4[1])/6) )

     return (step_position_solution, step_velocity_solution)



def compute_f_g_f_dot_g_dot(nu, delta_nu, p, e, pos_vec, mu) -> tuple:
     '''Takes in the current True Anomaly (Nu), and the delta True Anomaly (delta_nu) in degrees.
     Also takes in p, and the eccentricity (e).
     Also takes in the position vector (Generally the r_vector but can be any given position vector).
     Also takes in a value of Mu (Generally from WGS 84).'''
     radian_nu = math.radians(nu)
     radian_delta_nu = math.radians(delta_nu)

     r = p/(1+e*math.cos(radian_delta_nu + radian_nu))

     # Compute f
     f = 1 - (r/p)*(1 - math.cos(radian_delta_nu))
     
     # Compute g
     g = ((r*np.linalg.norm(pos_vec))/(math.sqrt(mu * p))) * math.sin(radian_delta_nu)

     # Compute g_dot
     g_dot = 1 - (np.linalg.norm(pos_vec)/p) * (1 - math.cos(radian_delta_nu))

     # Compute f_dot based on f, g, and g_dot
     f_dot = (f * g_dot - 1)/g

     return (f, g, f_dot, g_dot)



def compute_arbitrary_new_position(f, g, current_position, velocity_components):
     '''Determine the new position vector from a provided f and g function'''
     return np.array(current_position) * f + np.array(velocity_components) * g



def compute_arbitrary_new_velocity(f_dot, g_dot, current_position, velocity_components):
     '''Determine the new velocity vector from a provided f_dot and g_dot function'''
     return np.array(current_position) * f_dot + np.array(velocity_components) * g_dot



def convert_date_to_jd(year: int, month: int, day: int, hour: int, min: int, sec: int) -> float:
     '''Converts a provided date to Julian Date (JD). Returns a float with the JD in days'''
     timeut = hour + ( min / 60 ) + ( sec / 3600 )

     jd = (367 * year) - np.floor( 7 * ( year + np.floor( ( month + 9 ) / 12 ) ) / 4 ) -  np.floor( 3 * ( np.floor( ( year + ( month - 9 ) / 7 ) / 100 ) + 1 ) / 4 ) + np.floor( ( 275 * month ) / 9 ) + day + 1721028.5

     return jd + ( timeut / 24 )

def convert_jd_to_mjd(jd: float) -> float:
     '''Converts a provided julian date to Modified Julian Date (MJD). Returns a float with the MJD in days'''
     return jd - 2400000.5

def convert_mjd_to_besselian_year(mjd: float):
     '''Converts a Modified Julian Date (MJD) to Besselian Year'''
     return 2000.0 + (mjd - 51544.03)/365.242199

def convert_besselian_to_ut2_ut1(besselian_years: float):
     return 0.022 * math.sin(2*math.pi*besselian_years) - 0.012 * math.cos(2*math.pi*besselian_years) - 0.006 * math.sin(4*math.pi*besselian_years) + 0.007 * math.cos(4*math.pi*besselian_years)

def convert_ut1_utc(mjd, ut2_ut1):
     '''Takes in the Modified Julian Date and the UT2-UT1 date and returns the UT1-UTC date'''

     # From https://datacenter.iers.org/data/latestVersion/bulletinA.txt
     # Dated: 26 February 2026                                    Vol. XXXIX No. 009
     return 0.0640 + 0.00003 * (mjd - 61091) - (ut2_ut1)

def get_tai_utc() -> float:
     '''Taken from: https://datacenter.iers.org/data/latestVersion/bulletinA.txt  which at the time of writing (26 February 2026) defines TAI-UTC as 37.000000 seconds'''
     return 37.000000

def get_tt_utc() -> float:
     '''Returns the terrestrial time '''
     return get_tai_utc() + 32.184

def get_j2000_jd_epoch() -> float:
     '''Returns the Julian Date epoch for the j2000 specified epoch: 1 Jan 2000 12:00:00'''
     return convert_date_to_jd(2000, 1, 1, 12, 0, 0)

def determine_seconds_since_j2000_epoch(jd) -> float:
     '''Takes in a Julian Date and converts it to seconds since the J2000 epoch'''
     return (jd - get_j2000_jd_epoch()) * 86400


def determine_tai_s(seconds):
     '''Takes in a time in seconds since the j2000 epoch and converts to the TAI seconds'''
     return seconds + get_tai_utc()

def determine_tai(seconds):
     '''Takes in the TAI seconds and returns the TAI'''
     return get_j2000_jd_epoch() + seconds/86400

def determine_tbd(tai):
     '''Takes in the TAI and computes the Terrestrial Barycentric Dynamic (TBD) time'''
     return tai + 32.184/86400

def determine_tt(jd):
     '''Not to be confused with get_tt_utc(). This will calculate the Terrestrial Time for a given Julian Date'''
     return jd + get_tt_utc()/86400

def determine_sun_vector_lf(jd):
     '''Uses a low fidelity model from the Astronomical Almanac for 2017 for the apparent right ascension and declination of the sun'''
     au = 149597870700 # Definition of a single AU in meters
     n = jd - 2451545

     # L
     mean_longitude_of_sun = math.radians(280.460+0.9856474*n)
     mean_longitude_of_sun = mean_longitude_of_sun % (2*math.pi)
     if mean_longitude_of_sun < 0:
          mean_longitude_of_sun = mean_longitude_of_sun + 2*math.pi

     # g
     mean_anomaly_of_sun = math.radians(357.528+0.9856003*n)
     mean_anomaly_of_sun = mean_anomaly_of_sun % (2*math.pi)
     if mean_anomaly_of_sun < 0:
          mean_anomaly_of_sun = mean_anomaly_of_sun + 2*math.pi

     # Lambda
     ecliptic_longitude = mean_longitude_of_sun + math.radians(1.915)*math.sin(mean_anomaly_of_sun) + math.radians(0.020)*math.sin(2*mean_anomaly_of_sun)
     # Beta
     ecliptic_latitude = 0

     #epsilon
     obliquity_of_ecliptic = math.radians(23.439-0.0000004*n)

     # Alpha
     f = 180/math.pi
     t = math.pow(math.tan(obliquity_of_ecliptic/2), 2)
     right_ascension = ecliptic_longitude - f*t*math.sin(2*ecliptic_longitude)+(f/2)*math.pow(t, 2)*math.sin(4*ecliptic_longitude)
     # right_ascension = math.atan(math.cos(obliquity_of_ecliptic)*math.tan(ecliptic_longitude))

     distance_from_sun = 1.00014 - 0.01671*math.cos(mean_anomaly_of_sun) - 0.00014*math.cos(2*mean_anomaly_of_sun)

     # Sun Coordinates
     x = (distance_from_sun * math.cos(ecliptic_longitude))*au
     y = (distance_from_sun * math.cos(obliquity_of_ecliptic)*math.sin(ecliptic_longitude))*au
     z = (distance_from_sun * math.sin(obliquity_of_ecliptic)*math.sin(ecliptic_longitude))*au

     return [x, y, z]

def determine_moon_vector_lf(jd):
     '''Uses a low fidelity model from the Astronomical Almanac for 2017 for the apparent right ascension and declination of the sun'''
     er = 6378137 # meters
     t = (jd-2451545)/36525


     lambdalambda = math.radians(218.32+481267.881*t) + math.radians(6.29) * math.sin(math.radians(135.0 + 477198.87*t)) - math.radians(1.27) * math.sin(math.radians(259.3 - 413335.36*t))+ math.radians(0.66) * math.sin(math.radians(235.7 + 890534.22*t)) + math.radians(0.21) * math.sin(math.radians(269.9 + 954397.74*t)) - math.radians(0.19) * math.sin(math.radians(357.5 + 35999.05*t)) -  math.radians(0.11) * math.sin(math.radians(186.5 + 966404.03*t))

     beta = math.radians(5.13) * math.sin(math.radians(93.3 + 483202.02*t)) + math.radians(0.28)*math.sin(math.radians(228.2+960400.89*t))-math.radians(0.28) * math.sin(math.radians(318.3 + 6003.15*t))-   math.radians(0.17)*math.sin(math.radians(217.6 - 407332.21*t))

     new_pi = math.radians(0.9508)+math.radians(0.0518)*math.cos(math.radians(135.0+477198.87*t))+math.radians(0.0095)*math.cos(math.radians(259.3-413335.36*t))+math.radians(0.0078)*math.cos(math.radians(235.7+890534.22*t))+math.radians(0.0028)*math.cos(math.radians(269.9+954397.74*t))


     sd = 0.2724*new_pi
     r = 1/math.sin(new_pi)


     l = math.cos(beta)*math.cos(lambdalambda)
     
     m = 0.9175*math.cos(beta)*math.sin(lambdalambda)-0.3978*math.sin(beta)
     
     n = 0.3978*math.cos(beta)*math.sin(lambdalambda)+0.9175*math.sin(beta)

     
     alpha = math.atan(m/l)

     delta = math.asin(n)

     x = r*l*er
     y = r*m*er
     z = r*n*er

     return [x, y, z]
