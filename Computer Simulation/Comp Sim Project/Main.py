# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 17:48:46 2022

@author: Lukah
"""

from LaunchSim import LaunchSim as LS
from BodySim import BodySim as BS
from Body import Body
import math
import os

def main():
    
    Bodies = []
    timestep = 18000 #Speed of Simulation\Time Interval between Updates (sec)
    time = 40000 #Length of Simulation (multiply by timstep to get time of simulation)
    time2 = 1000 #shorter length used for many repeated simulations 
    earthEscapeVel = 11200 #an example value
    sunEscapeVel = 615000 #an example value

    #change to path for this file's location (not including file)
    directory = "C:\\Users\Lukah\Documents\Comp Sim Project\Comp Sim Project"

    #make this directory the working directory (just to minimise issues, not really necessary)
    os.chdir(directory) 
    #assigning variables for folder locations of this project
    folder1 = os.path.join(directory, "Planet Data")
    folder2 = os.path.join(directory, "Satellite Data")
    folder3 = os.path.join(directory, "Output Data")

    #loop through files in folder
    for file in os.listdir(folder1):
        
        filename = os.path.join(folder1, file)
        Bodies.append(Body(filename, timestep)) #append a Body object to a list using the file

    #Create new Body object file for later use
    newfile = Body.writeFile(folder2, "Satellite1", "2.180 * 10 ** 3", "(1.496 * 10 ** 11 + 6.378 * 10 ** 6, 0)", "(0, 0)", "w", "5 * (10 ** 6)")
    newFilePath = os.path.join(folder2, newfile)
    satellite = Body(newFilePath, timestep)

    #Bodies.append(b)
    a = BS(Bodies, timestep, maxTime=time, repeat=False, directory=folder3)
    a.display()
    a.energyGraph()
    
    #l = LS(Bodies, satellite, timestep, [earthEscapeVel, sunEscapeVel/10], [0, 0.5 * math.pi], maxTime=time2)
    #l.launchTo("Mars", "Earth")

main()