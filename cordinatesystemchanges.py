import math
import numpy as np


def orbital_elements(r, v):
    """ Given the position and velocity vectors, orbital elements can be found using this function"""

    R = np.linalg.norm(r)  # magnitude of position vector, distance
    V = np.linalg.norm(v)  # magnitude of velocity vector, speed

    v_r = np.dot(r, v) / R  # radial velocity, sign gives the direction of motion towards or away from perigee
    print(v_r)
    print(len(v), len(r))
    h = np.cross(r, v)  # specific angular momentum
    H = np.linalg.norm(h)  # magnitude of specific angular momentum

    i = math.degrees(math.acos(h[2]/H))  # angle of inclination of the orbit to the equatorial frame

    k_cap = np.array([0, 0, 1])  # unit vector in k direction
    n = np.cross(k_cap, h)  # this gives the node line
    N = np.linalg.norm(n)

    if n[1] >= 0:
        ohm = math.degrees(math.acos(n[0]/N))  # right ascension angle of the orbit if y-component is positive
    else:
        ohm = 360 - math.degrees(math.acos(n[0]/N))  # if its negative, since right ascension is between 0 and 180, subtract it from 360

    e = (np.cross(v, h) - (398600 * (r / R)))/398600  # eccentricity vector
    E = np.linalg.norm(e)  # magnitude of eccentricity

    if e[2] >= 0:
        omega = math.degrees(math.acos(np.dot(n, e)/(N*E)))  # argument of perigee
    else:
        omega = 360 - math.degrees(math.acos(np.dot(n, e)/ (N * E)))

    if v_r >= 0:
        theta = math.degrees(math.acos(np.dot(e, r)/(E * R)))  # true anomaly
    else:
        theta = 360 - math.degrees(math.acos(np.dot(e, r)/(E * R)))  # true anomaly

    return h, n, e, i, ohm, omega, theta


def geocentric_perifocal_frames(state, right_ascension, inclination, argument_of_perigee, to_geocentric=True):
    """ Given a state vector (column vector), and the correct angles (degrees), this function converts the state vector to
    geocentric or peri-focal frames whenever mentioned in the binary variable 'to_geocentric' """

    # ohm = math.radians(right_ascension)
    # i = math.radians(inclination)
    # omega = math.radians(argument_of_perigee)

    ohm = right_ascension
    i = inclination
    omega = argument_of_perigee

    r_3_ra = np.array([[math.cos(ohm), math.sin(ohm), 0],  # rotation along z-axis to align node line and x-axis
           [- math.sin(ohm), math.cos(ohm), 0],
           [0, 0, 1]])

    r_1_i = np.array([[1, 0, 0],  # rotation along x-axis to align the peri-focal angular momentum vector  and z-axis
             [0, math.cos(i), math.sin(i)],
             [0, - math.sin(i), math.cos(i)]])

    r_3_ap = np.array([[math.cos(omega), math.sin(omega), 0],  # rotation along z-axis to align the eccentricity
              [- math.sin(omega), math.cos(omega), 0],  # vector (perigee) and y-axis
              [0, 0, 1]])

    q = np.dot(np.dot(r_3_ap, r_1_i), r_3_ra)  # notice the order of multiplication

    if to_geocentric:
        return np.dot(q.T, state)  # transpose of the matrix should be multiplied for the transformation of coordinates
    return np.dot(q, state)


# if __name__ == "__main__":
#     a = np.array([6285, 3628, 0])
#     print(geocentric_perifocal_frames(a, 40, 30, 60))






