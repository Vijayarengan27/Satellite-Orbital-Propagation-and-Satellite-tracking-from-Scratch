import datetime
import math
import numpy as np

import requests
import elliptical_orbit as eo
import script as sc
import cordinatesystemchanges as coc


def get_tle_data():
    """ To get the tle data from celestrack. Use it once a day as it updates once or twice a day"""

    url = "https://celestrak.org/NORAD/elements/stations.txt"
    response = requests.get(url)
    tle_data = response.text
    print(tle_data)
    with open('tle.txt', 'w') as f:
        for line in tle_data.splitlines():
            f.write(line)
            f.write('\n')
    print('Data is added to the file')
    f.close()


# get_tle_data()  # run this to get the tle data of all satellites

stations = []
station_tle_dict = {}
with open('tle.txt', 'r') as tle_data:
    for i, line in enumerate(tle_data):
        if (i + 1) % 3 == 1:
            stations.append(line.strip())

        elif (i + 1) % 3 == 2:
            k = stations[-1] + '_line_1'
            station_tle_dict[k] = line
        else:
            k = stations[-1] + '_line_2'
            station_tle_dict[k] = line
# needed values

station = stations[20]  # select any station here
print('The selected satellite is', station)
epoch_time_string = station_tle_dict[f"{station}_line_1"].split()[3]
eccentricity = float(station_tle_dict[f"{station}_line_2"].split()[4]) * 10 ** -6
initial_mean_anomaly = np.radians(float(station_tle_dict[f"{station}_line_2"].split()[6]))
mean_motion = float(str(station_tle_dict[f"{station}_line_2"].split()[7])[:-5])
argument_of_perigee = math.radians(float(station_tle_dict[f"{station}_line_2"].split()[5]))
inclination = math.radians(float(station_tle_dict[f"{station}_line_2"].split()[2]))
right_ascension = math.radians(float(station_tle_dict[f"{station}_line_2"].split()[3]))
current_time = datetime.datetime.now(tz=datetime.timezone.utc)




# Flow of process:
# get the mean anomaly,
m_a = sc.mean_anomaly(initial_mean_anomaly, mean_motion,epoch_time_string,current_time)

# then eccentric and true anomaly, (take true anomaly from the next function)
e_a, t_a_1 = eo.eccentric_and_true_from_mean(m_a, eccentricity)

# then satellite's position in the
# peri-focal frame,
r_p, t_a, v = sc.true_anomaly_position(e_a, eccentricity, mean_motion)

# then change it to geocentric frame, then calculate overall theta change of greenwich to the
# propagated time,then rotate the position of satellite further to that theta.
r_g = coc.geocentric_perifocal_frames(r_p, right_ascension, inclination, argument_of_perigee)
v_g = coc.geocentric_perifocal_frames(v, right_ascension, inclination, argument_of_perigee)

theta_g = sc.adjusting_earth_rotation(epoch_time_string, current_time)

r_geo = sc.turning_coordinates(theta_g, r_g)
v_geo = sc.turning_coordinates(theta_g, v_g)

# Get the accurate observer position with geodetic latitude and calculate the difference in satellite and observer
# position.
r_obs_geo = sc.calc_observer_position(math.radians(10.82), math.radians(78.7),0, epoch_time_string, current_time)


# Then convert this range vector to nez or sez frame (topo-centric horizontal).
# Then at last find the angle of elevation and when it is 90 degrees, it is overhead.
a_e, a_a = sc.calc_elevation(math.radians(10.82), math.radians(78.7), r_geo, r_obs_geo, v_geo)
print('the angle of elevation in sez frame from the observer is', a_e)
print('the angle of azimuth from 0 degree south from the observer is', a_a)

# print('the latitude and longitude of ground track of satellite is', sc.calc_satellite_lat_long(r_geo))








