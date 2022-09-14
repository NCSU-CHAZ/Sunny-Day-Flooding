# Script: pullTime_plot.py
# Author Name: Jenero Knowles
# Date of Authorship: 07.29.2022
# Last Modified: 08.26.2022
# Purpose: Create plots of ADCIRC results against observed data.
#
# Software Requirements: pullTime.py, observation station name and ID, observed data file, measured data file
#
# Output: JPEG image of a plot of ADCIRC results against observed data.
# -*- coding: utf-8 -*-

""" Import packages to script. """

import os
import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt

""" Read mesh and element from ADCIRC file. 
        Parameter (var, units) - mesh
        water level (zeta, m) - fort.63.nc
        wind speed (windx/windy, m/s) - fort.74.nc
        water velocity (u-vel/v-vel, m/s) - fort.64.nc 
        mean wave direction (swan_DIR, degrees) - swan_DIR.63.nc
        significant wave height (swan_HS, m) - swan_HS.63.nc
        peak wave period (swan_TPS, s) - swan_TPS.63.nc  """

# Remove " " from full paths
m_data = input("Enter the full path of modeled/measured data csv file (pullT_zeta_3.csv): ")


''' Assign plot variables with units from filename.'''

if "zeta" in m_data:
    var = "Water Level"
    units = "m MSL"
    ob_stat = input("Enter the ID and name of the observation station (8638610 Sewells Point): ")
    # Remove " " from full paths
    ob_data = input("Enter the full path of observed data csv file (CO-OPS_8638610_met.csv): ")
    """if water levels use this for observed data"""
    x_o = []
    y_o = []
    df = pd.read_csv(ob_data, parse_dates={"Comb_Date": ["Date", "Time (GMT)"]})        # Combine Date&Time
    df = df.drop(df.columns[[1, 2]], axis=1)        # Remove all other columns not water levels
    x_o = df["Comb_Date"].tolist()                  # Add date items to list for plot
    y_o = df["Verified (m)"].tolist()               # Add verified water levels to list for plots

    """ modeled water levels from time step. """
    x_m = []
    y_m = []

    outputFileName = os.path.splitext(m_data)[0] + "_rev.csv"
    with open(m_data, 'r') as inFile, open(outputFileName, 'w') as outfile:     # Copy csv to make edits in the form
        r = csv.reader(inFile)                                                  # of a dataframe
        w = csv.writer(outfile)
        w.writerow(['Timestamp', 'DateTime', 'Station No.', 'Modeled'])         # Add titles to columns to align cells
        for row in r:
            if any(field.strip() for field in row):
                w.writerow(row)

    df_2 = pd.read_csv(outputFileName)
    df_2['Timestamp'].replace(' ', np.nan, inplace=True)
    df_2['DateTime'] = pd.to_datetime(df_2['DateTime'])     # Change datatype from string to Panda Timestamp
    df_2 = df_2.drop("Station No.", axis=1)                 # Remove all other columns not water levels
    x_m = df_2["DateTime"].tolist()
    y_m = df_2["Modeled"].tolist()

elif "wind" in m_data:
    var = "Wind Speed"
    units = "m/s"
    ob_stat = input("Enter the ID and name of the observation station (8632200 Kiptopeke Beach): ")
    # Remove " " from full paths
    ob_data = input("Enter the full path of observed data csv file (CO-OPS_8632200_met.csv): ")
    """if wind speed, use this for observed data"""
    x_o = []
    y_o = []
    df = pd.read_csv(ob_data, parse_dates={"Comb_Date": ["Date", "Time (GMT)"]})  # Combine date & time into one column
    df = df.drop(df.columns[[2, 3, 4, 5, 6, 7]], axis=1)    # Remove all other columns not related to wind speed
    for ix, x in enumerate(df["Wind Speed (m/s)"]):         # Check Nans in the form of empty cells and fill with -9
        if x == "-":
            df.loc[ix, "Wind Speed (m/s)"] = float(-9)      # Change -9 t0 -99 for larger range if needed
    df["Wind Speed (m/s)"] = df["Wind Speed (m/s)"].astype(float)       # Change datatype from object to float
    df.drop(df[df["Wind Speed (m/s)"] <= 0].index, inplace=True)
    x_o = df["Comb_Date"].tolist()
    y_o = df["Wind Speed (m/s)"].tolist()

    """ modeled wind speed from time step. """
    x_m = []
    y_m = []
    outputFileName = os.path.splitext(m_data)[0] + "_rev.csv"
    m_data_y = input("Enter the full path of modeled/measured data csv file or y-file (pullT_windy_3.csv): ")
    outputFileName_y = os.path.splitext(m_data_y)[0] + "_rev.csv"

    with open(m_data, 'r') as inFile, open(outputFileName, 'w') as outfile:     # Copy csv to make edits in the form
        r = csv.reader(inFile)                                                  # of a dataframe
        w = csv.writer(outfile)
        w.writerow(['Timestamp', 'DateTime', 'Station No.', 'Modeled_x'])       # Add titles to columns to align cells
        for row in r:
            if any(field.strip() for field in row):
                w.writerow(row)

    with open(m_data_y, 'r') as inFile, open(outputFileName_y, 'w') as outfile:  # Copy csv to make edits in the form
        r = csv.reader(inFile)                                                   # of a dataframe
        w = csv.writer(outfile)
        w.writerow(['Timestamp', 'DateTime', 'Station No.', 'Modeled_y'])   # Add titles to columns to align cells
        for row in r:
            if any(field.strip() for field in row):
                w.writerow(row)

    df_x = pd.read_csv(outputFileName)
    df_x['Timestamp'].replace(' ', np.nan, inplace=True)
    df_x['DateTime'] = pd.to_datetime(df_x['DateTime'])             # Change datatype from string to Panda Timestamp
    df_x = df_x.drop(["Station No.", "Timestamp"], axis=1)          # Remove all other columns not water levels
    df_y = pd.read_csv(outputFileName_y)
    df_y['Timestamp'].replace(' ', np.nan, inplace=True)
    df_y['DateTime'] = pd.to_datetime(df_y['DateTime'])             # Change datatype from string to Panda Timestamp
    df_y = df_y.drop(["Station No.", "Timestamp"], axis=1)          # Remove all other columns not water levels

    df_2 = df_x
    df_2 = df_2.join(df_y["Modeled_y"])                             # Add the y-component to the dataframe
    df_2['Magnitude'] = np.linalg.norm(df_2[['Modeled_x', 'Modeled_y']].values, axis=1)  # Compute magnitude of x&y
    #   df_2['Magnitude'] = (df_2['Modeled_x']**2 + df_2['Modeled_y']**2)**(1/2)
    x_m = df_2["DateTime"].tolist()
    y_m = df_2["Magnitude"].tolist()

elif "vel" in m_data:
    var = "Water Velocity"
    units = "m/s"
elif "DIR" in m_data:
    var = "Mean Wave Direction"
    units = "Degrees"
elif "HS" in m_data:
    var = "Significant Wave Height"
    units = "m"
    ob_stat = input("Enter the ID and name of the observation station (41001 East Hatteras): ")
    # Remove " " from full paths
    ob_data = input("Enter the full path of observed data csv file (41001h2011.txt): ")
    """if wave height, use this for observed data"""
    x_o = []
    y_o = []
    df = pd.read_csv(ob_data, delim_whitespace=2, skiprows=1)
    df = df.astype(str)
    df["Date"] = df[["#yr", "mo", "dy"]].agg('-'.join, axis=1)
    df["Time"] = df[["hr", "mn"]].agg(':'.join, axis=1)
    df["DateTime"] = df["Date"] + ' ' + df["Time"]
    df["DateTime"] = pd.to_datetime(df["DateTime"])
    df.drop(df[df["DateTime"] >= "2011-09-01"].index, inplace=True)
    df.drop(df[df["DateTime"] < "2011-08-21"].index, inplace=True)
    df["m"] = df["m"].astype(float)
    df.drop(df[df["m"] >= 50].index, inplace=True)
    x_o = df["DateTime"].tolist()
    y_o = df["m"].tolist()

    """ modeled wave heights from time step. """
    x_m = []
    y_m = []

    outputFileName = os.path.splitext(m_data)[0] + "_rev.csv"

    with open(m_data, 'r') as inFile, open(outputFileName, 'w') as outfile:
        r = csv.reader(inFile)
        w = csv.writer(outfile)
        w.writerow(['Timestamp', 'DateTime', 'Station No.', 'Modeled'])
        for row in r:
            if any(field.strip() for field in row):
                w.writerow(row)

    df_2 = pd.read_csv(outputFileName)
    df_2['Timestamp'].replace(' ', np.nan, inplace=True)
    df_2['DateTime'] = pd.to_datetime(df_2['DateTime'])
    df_2 = df_2.drop("Station No.", axis=1)
    x_m = df_2["DateTime"].tolist()
    y_m = df_2["Modeled"].tolist()

else:
    var = "Peak Wave Period"
    units = "s"


""" Plot comparison lines on chart. """

plt.plot(x_m, y_m, label="Modeled", linestyle="-")
plt.scatter(x_o, y_o, label="Observed", s=2, c="grey")
plt.gcf().autofmt_xdate()
plt.xlabel('Date')
plt.ylabel('{} ({})'.format(var, units))
if var == "Wind Speed":
    plt.ylim(bottom=0)
plt.title('{} at {}'.format(var, ob_stat))
plt.legend()
plt.grid()
plt.savefig("OutputFiles" + "/" + "{}_{}.png".format(var, ob_stat), dpi=300, bbox_inches='tight')
plt.show()
