from keplarianElements import KeplerianElements
import keHelperFunctions

import math
import datetime

import tabulate

def main():
    
    utc_time = datetime.datetime(2026, 2, 28, 18, 22, 45)
    print(f'Input time: {utc_time.date()} {utc_time.time()}')

    # Compute Julian Date
    jd = keHelperFunctions.convert_date_to_jd(utc_time.year, utc_time.month, utc_time.day, utc_time.hour, utc_time.minute, utc_time.second)
    print(f'Julian Date: {jd}')
    mjd = keHelperFunctions.convert_jd_to_mjd(jd)
    print(f'Modified Julian Date: {mjd}')
    besselian_year = keHelperFunctions.convert_mjd_to_besselian_year(mjd)
    print(f'Besselian Year: {besselian_year}')
    ut2_ut1 = keHelperFunctions.convert_besselian_to_ut2_ut1(besselian_year)
    print(f'UT2-UT1: {ut2_ut1}')


    # Compute UT1-UTC
    ut1_utc = keHelperFunctions.convert_ut1_utc(mjd, ut2_ut1)
    print(f'UT1-UTC: {ut1_utc}')


    # Compute TAI Julian Date
    tai_utc = keHelperFunctions.get_tai_utc()
    print(f'TAI-UTC: {tai_utc}')
    tt_utc = keHelperFunctions.get_tt_utc()
    print(f'TT-UTC: {tt_utc}')


    # Compute TT/TDB Julian Date
    seconds_since_j2000 = keHelperFunctions.determine_seconds_since_j2000_epoch(jd)
    print(f'Seconds since J2000: {seconds_since_j2000}')
    tai_s = keHelperFunctions.determine_tai_s(seconds_since_j2000)
    print(f'TAI Seconds: {tai_s}')
    tai = keHelperFunctions.determine_tai(tai_s)
    print(f'TAI: {tai}')
    tbd = keHelperFunctions.determine_tbd(tai)
    print(f'TBD: {tbd}')
    tt = keHelperFunctions.determine_tt(jd)
    print(f'TT: {tt}')

    # Using low-fidelity formula, compute TOD Sun vector
    sun_coords = keHelperFunctions.determine_sun_vector_lf(tt)
    print(f'Sun Vector: {sun_coords}')


    # Using low-fidelity formula, compute TOD Moon vector
    moon_coords = keHelperFunctions.determine_moon_vector_lf(tt)
    print(f'Moon Vector: {moon_coords}')


if __name__ == '__main__':
    main()
