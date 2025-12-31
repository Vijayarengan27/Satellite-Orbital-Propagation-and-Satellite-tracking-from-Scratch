# Orbital Mechanics

A small project to propagate satellite position from the TLE(two line element) and then track it from observer position using 
first principles without using black-box libraries.

## Overview

We take TLE of about 32 satellites from open source NORAD website and then propagate its motion and then find its relative position
from the observer. Their data is present in tle.txt file.

This project builds an end-to-end satellite tracking pipeline from scratch:

- Parsing TLE data
- Time handling (Julian date)
- Orbital anomaly calculations
- Earth rotation and coordinate transformations
- Observer-relative tracking in the SEZ (South-East-Zenith) frame

Finished work:
- Orbit calculation
- Simple propagation
- Satellite tracking from observer position in SEZ (South East Zenith) frame

## Accuracy

The cross-checking is done by using data from in-the-sky.org which shows the current satellite position map. The distance to satellite 
measurement is approximately within +/- 150 kms range and the azimuth angle and the angle of elevation is approximately within +/- 2 degrees.
The reason for this is because of the air drag and then the gravitational perturbations which needs to be accounted. It is in the
future work pipeline.

## How to run

- In main.py, uncomment the line 26 just the first time running to get the latest TLE data.
- Then select the required satellites available from the TLE.txt. For example, if you want to see ISS (NAUKA), it is the
3rd satellite in the data, so select 2 in line 43.
- Change your latitude and longitude in line 80 respectively in calculate_observer_position function call.
- Then run main.py
- To cross verify, check the angle of elevation, azimuth angle in sez (South East Zenith) and the distance with in-the-sky.org,where you
can search for a satellite and changing your observer latitude and longitude.

Here is the example of tracking verification of ISS NAUKA,
## From in-the-sky.org

![in-the-sky.org screenshot](images/in-the-sky.org%20example.png)

## From the code

![Code output](images/Output%20example.png)

Look how the distance of satellite from observer, angle of elevation in sez frame and azimuth angles are all similar.

## Future Work

Improve the accuracy further, polish the code, and create a visualization of satellites.

## References

1. Howard D. Curtis, *Orbital Mechanics for Engineering Students*, Elsevier  
   — used as the primary reference for orbital propagation and coordinate transformations.

2. in-the-sky.org  
   — used for cross-verification of satellite position, azimuth, and elevation data.

Thanks for reading, suggestions and feedback are always welcomed.





