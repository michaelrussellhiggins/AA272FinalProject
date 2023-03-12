# AA272 Final Project - Applying GPS/IMU Sensor Fusion to Athlete Tracking

The primary objective of this project was to use GPS/IMU sensor fusion to track American football receivers and cornerbacks through a number of differnt routes. Sensor data was collected using the GNSSLogger app, found on the Google Play Store.

The Data folder contains much of the raw (.txt) and processed GPS (.csv) and IMU (.xlsx) data.

The ClipData folder contains information regarding the routes run, including the approximated ground truth data and the times corresponding to the beginning and end of each route run.

The Ephem folder contains the ephemeris data corresponding to the time of the experiment, and is sourced from https://cddis.nasa.gov.

The HelperScripts folder contains miscellaneous functions used by several scripts.

The opensource folder contains open source scripts developed by Google for use with their GNSSLogger app. The original version of this folder can be found at https://github.com/google/gps-measurement-tools. SOme of these scripts were modified for the purposes of this project.

Many of the remaining scripts relate to the processing of the raw data into a usable format for the sensor fusion algorithms. Much of the data processing is handled by the DataExtraction.m script, which calls many of the scripts included in the top level of the repository.
