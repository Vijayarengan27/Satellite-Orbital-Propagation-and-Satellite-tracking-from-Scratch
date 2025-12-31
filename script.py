import math
import numpy as np

import datetime


def epoch_time_difference(daystring, current_time, julian=False):
    """returns the time for the given daystring from TLE in %Y-%m-%d %H:%M:%S"""

    y = int(str(daystring)[:2])
    day = float(str(daystring)[2:])

    if y <= 57:  # to differentiate between satellites launched in the 1900's and 2000's
        y += 2000
    else:
        y += 1900

    epoch_datetime = datetime.datetime(y, 1, 1) + datetime.timedelta(days=day - 1)

    if julian:
        y = epoch_datetime.year
        m = epoch_datetime.month
        d = epoch_datetime.day

        z = y + int((m + 9) / 12)
        j_o = (367 * y) - (7 * z / 4) + int(275 * m / 9) + d + 1721013.5
        ut = epoch_datetime.hour + epoch_datetime.minute/60 + epoch_datetime.second/3600
        return j_o, ut

    naive_current_time = current_time.replace(tzinfo=None)  # done to make the time offset naive corresponding to the
    # epoch date time

    time_diff = float((naive_current_time - epoch_datetime).total_seconds())
    return time_diff


def mean_anomaly(m_o, n, t_o, current_time):
    """ The function calculates mean anomaly of a satellite given its initial mean anomaly, mean motion and the
    reference time in string. Mean motion is the angular velocity. The units for all in seconds and radians"""

    time_diff = epoch_time_difference(t_o, current_time)  # in seconds
    n_rad_s = (n * 2 * math.pi)/(24 * 60 * 60)
    m = m_o + (n_rad_s * time_diff)

    m = m % (2 * math.pi)  # normalization

    return m


def true_anomaly_position(eccentric_anomaly, e, n):
    """given eccentric anomaly and eccentricity, this returns the true anomaly of the satellite and its position in
     its orbital frame"""

    t_a = 2 * math.atan2(((1 + e) ** 0.5) * math.tan(eccentric_anomaly/2), (1 - e) ** 0.5)
    t_a = t_a % (2 * math.pi)

    a = (398600.4418 * ((86400 / (2 * math.pi * n)) ** 2)) ** (1/3)  # time period = (2 pi / root of mu) semi major axis to the power 3/2
    r_s = a * (1 - e * math.cos(eccentric_anomaly))

    r = np.array([r_s * math.cos(t_a), r_s * math.sin(t_a), 0])

    h = (398600.4418 * a * (1 - (e ** 2))) ** 0.5  # calculate angular momentum to find out velocity vector
    z = 398600.4418 / h

    v = np.array([-z * math.sin(t_a), z * (e + math.cos(t_a)), 0])

    return r, t_a, v


def adjusting_earth_rotation(t_o, current_time):
    """ to adjust earth rotation as we predict in future from the measurement, we rotate the earth for that time
    to find the correct greenwich meridian position"""

    time_diff = epoch_time_difference(t_o, current_time)
    omega = 7.292115 * (10 ** -5)  # earth's angular velocity in rad/s

    # find Julian days j_o
    j_o, ut = epoch_time_difference(t_o, current_time, julian=True)

    # find T_o corresponding to it
    T_o = (j_o - 2451545) / 36525

    # find greenwich sidereal time theta_g_o, IN DEGREES
    total_theta_g_o = 100.4606184 + (36000.77004 * T_o) + (0.000387933 * (T_o ** 2)) - (2.583 * (10 ** -8) * (T_o ** 3))

    theta_g_o = math.radians((total_theta_g_o % 360) + (360.98564724 * ut/24))  # check whether to modulo here

    total_theta_g = theta_g_o + (omega * time_diff)
    theta_g = total_theta_g % (2 * math.pi)

    return theta_g


def turning_coordinates(theta_g, r):
    """To change coordinate system of the given position vector to a given theta for getting a coordinate system for
    a meridian or a longitude"""
    r_3_theta_g = np.array([[math.cos(theta_g), math.sin(theta_g), 0],  # rotation along z-axis to align the greenwich
                       [- math.sin(theta_g), math.cos(theta_g), 0],  # meridian for getting the overhead satellite position
                       [0, 0, 1]])
    return np.dot(r_3_theta_g, r)  # returns the position with respect to the rotated frame


def calc_observer_position(phi_dash, theta, H, t_o, current_time):
    """Given latitude, longitude and altitude, this function returns the accurate observer position vector using the
    ellipsoidal earth model. All in radians here and kilometers for altitude. t_o is added to get local sidereal time"""

    # theta_l = adjusting_earth_rotation(t_o, current_time) + theta # not needed as observer is fixed in ECEF system
    theta_gl = theta

    f = 1/298.257223563  # factor of oblateness or flattening
    Re = 6378.137  # radius at equator

    phi = math.atan(((1 - f) ** 2) * math.tan(phi_dash))  # geodetic latitude in radians
    e = (2 * f) - (f ** 2)
    z = 1 - (e * (math.sin(phi) ** 2))

    r = np.array([((Re/(z ** 0.5)) + H) * math.cos(phi) * math.cos(theta_gl), ((Re/(z ** 0.5)) + H) * math.cos(phi) * math.sin(theta_gl),
                  ((Re * ((1 - f) ** 2)/(z ** 0.5)) + H) * math.sin(phi)])  # position vector

    return r


def calc_elevation(phi_dash, theta, r, r_o, v):
    """ first get the observer position,r_o, using latitude, longitude and altitude and subtract it from the satellite's
    position r, then calculate its elevation,if its 90 degrees, then it is directly overhead to the observer position"""

    rho = r - r_o
    # convert rho to sez or nez topo-centric horizontal system and get the angle of elevation

    omega = 7.292115 * (10 ** -5)  # earth's angular velocity in rad/s
    v_obs = omega * r_o
    v_rho = v - v_obs

    f = 1/298.257223563  # factor of oblateness or flattening
    phi = math.atan(((1 - f) ** 2) * math.tan(phi_dash))  # use geodetic latitude for accurate

    q_sez = np.array([[-math.sin(phi)*math.cos(theta), -math.sin(phi)*math.sin(theta), math.cos(phi)],
                      [-math.sin(theta), math.cos(theta), 0],
                      [math.cos(phi)*math.cos(theta), math.cos(phi)*math.sin(theta), math.sin(phi)]])

    rho_sez = np.dot(q_sez, rho)
    v_rho_sez = np.dot(q_sez, v_rho)

    angle_of_elevation = math.degrees(math.asin(rho_sez[2]/np.linalg.norm(rho_sez)))
    angle_of_azimuth = math.degrees(math.atan2(rho_sez[1], rho_sez[0]))  # look for quadrant ambiguity
    if angle_of_azimuth < 0:
        angle_of_azimuth += 360
    real_a_of_a = math.degrees(math.atan2(rho_sez[1], - rho_sez[0]))
    print('the real angle of azimuth from 0 degree north is', real_a_of_a)
    print('the velocity of satellite from observer position in km/s in sez frame is', v_rho_sez)
    print('the distance of satellite from the observer in kms is', np.linalg.norm(rho_sez))
    return angle_of_elevation, angle_of_azimuth


def calc_satellite_lat_long(r):
    """Calculate the latitude and longitude of the satellite overhead on the earth from the position of satellite in
    ecef"""

    longitude = math.atan2(r[1], r[0])  # direct tan inverse gives the longitude
    a = 6378.137  # radius of earth
    e_squared = 0.00669438  # square of earth's eccentricity
    e_dash_squared = 0.00673949
    b = 6356.752  # polar radius
    d = r[0] ** 2 + r[1] ** 2
    theta = math.atan2(r[2] * a, d * b)
    latitude = math.atan2(r[2] + (e_dash_squared * b * (math.sin(theta)) ** 3), d - (e_squared * a * (math.cos(theta)) ** 3))  # wrong here
    return math.degrees(latitude % (2 * math.pi)), math.degrees(longitude % (2 * math.pi))


# if __name__ == '__main__':

    # m = (mean_anomaly(math.radians(30),11,25076.5))
    # t = (eo.eccentric_and_true_from_mean(m,0.055))
    # print(true_anomaly_position(t[0], 0.055, 11))
    # print(t[1])