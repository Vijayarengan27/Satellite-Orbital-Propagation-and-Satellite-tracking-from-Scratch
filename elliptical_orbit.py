import matplotlib.pyplot as plt
import numpy as np
import math


def mean_anomaly_from_radii(r_p, r_a, time):
    """This function calculates mean anomaly from getting the time elapsed from perigee, the perigee and apogee
    distances of a geocentric orbit"""

    e = (r_a - r_p) / (r_a + r_p)  # calculate eccentricity
    mu = 398600  # geocentric orbit with negligible satellite mass

    h = (r_p * mu * (1 + e)) ** 0.5  # angular momentum which is constant throughout

    T = (2 * math.pi / (mu ** 2)) * ((h /(1 - e**2)**0.5)**3)  # time period of the orbit

    mean_anomaly = 2 * math.pi * time / T

    return T, mean_anomaly, e


def eccentric_and_true_from_mean(mean_a, e):
    """ This function uses newton's method of finding the root of the equation to find the eccentric anomaly, if Mean
    anomaly and eccentricity are given. The initialization of ith value of eccentric anomaly is also provided. """

    # initialization
    if mean_a <= math.pi:
        eccentric_a_i = mean_a + e/2
    else:
        eccentric_a_i = mean_a - e/2

    f_i = eccentric_a_i - e * math.sin(eccentric_a_i) - mean_a
    f_dash_i = 1 - e * math.cos(eccentric_a_i)

    while abs(f_i/f_dash_i) > 10 ** -8:
        eccentric_a_i_1 = eccentric_a_i - (f_i/f_dash_i)
        eccentric_a_i = eccentric_a_i_1
        f_i = eccentric_a_i - e * math.sin(eccentric_a_i) - mean_a
        f_dash_i = 1 - e * math.cos(eccentric_a_i)

    tan_true_a = (((1+e)/(1-e))**0.5) * math.tan(eccentric_a_i/2)
    true_a = 2 * np.atan(tan_true_a)  # check here
    return eccentric_a_i, true_a


# def calc_angular_momentum(r,v):
#     h = np.cross(r,v)  # cross product of matrices of vector compenents of position and velocity(tangential)
#     return h
#
#
# def calc_elliptical_parameters(r,v,e):
#     mu =  # take from the book, assuming the satellite's mass is negligible
#
#     h = calc_angular_momentum(r, v)
#
#     r_p = (h **2 / mu) * (1/(1 + e))  # perigee radius
#     r_a = (h **2 / mu) * (1/(1 - e))  # apogee radius
#
#     a = (r_p + r_a) / 2  # semi-major axis
#     b = a * ((1 - e **2) ** 0.5)  # semi-minor axis
#
#     f = a*e  # distance of focus from the origin
#
#     return a, b
#
# def plot_elliptical_orbit(a, b):
#
#     # plot this function after x ** 2 / a ** 2 + y ** 2/ b ** 2 = 1


# if __name__ == '__main__':
#     print(mean_anomaly_from_radii(9600, 21000,10800))
#     print(eccentric_and_true_from_mean(3.6029, 0.3725))