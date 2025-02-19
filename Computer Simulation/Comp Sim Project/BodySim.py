# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 17:49:29 2022

@author: Lukah
"""

from pickle import NONE
from Body import Body
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from numpy.linalg import norm


class BodySim(object):
    """
    Simulation of Body objects.

    Parameters
    ----------
    object : object
        Class inheritance
    """

    def __init__(self, bodies, timestep, maxTime=100, repeat=True, directory=None):
        """
        Initialisation of instance of BodySim object for Simulation of Body objects.

        Parameters
        ----------
        bodies : list
            List of Body objects
        timestep : float
            Timestep of simulation (time between updates, in seconds)
        maxTime : int, optional
            Maximum number of updates run by simulation, by default 100
        repeat : bool, optional
            Whether or not the simulation has a maximum number of itterations/updates, or keeps going until proccess is stopped, by default True
        directory : str, optional
            Folder name for storing data, by default None
        """
        #reject anything with less than two items, or is not iterable
        if len(bodies) < 2:

            try:

                iterable = iter(bodies)
            
            except:

                raise TypeError("Type " + str(type(bodies) + " not iterable"))

        #loop through input list of bodies
        for i in range(len(bodies)):

            #check if variable is of type Body
            if isinstance(bodies[i], Body):

                continue

            else:

                #raised if one of the arguements was not a Body object
                raise TypeError("Type \"" + str(type(bodies[i])) + "\", not of type \"Body\"")

        self.bodies = bodies
        self.timestep = timestep
        self.maxTime = maxTime
        self.repeat = repeat
        self.directory = directory

        #initialising some stuff
        self.time = 0
        self.energyData = []
        self.orbitalPeriods = [0.0] * len(self.bodies)
        
        print("New Simulation\n\n")

        #initialise acceleration
        self.acceleration(1)
        
        
    def acceleration(self, x=2):
        """
        Updates the acceleration for each Body using an integration scheme, and stores the values.

        Parameters
        ----------
        x : int, optional
            0=timstep before, 1=current timstep, 2=next timestep. The default is 2.

        Returns
        -------
        None.

        """
        
        G = Body.G
        self.displacements = np.zeros((len(self.bodies), len(self.bodies), 2))

        #creating a matrix to contain position vectors of one body on another
        for i in range(len(self.displacements)):
            
            for j in range(len(self.displacements[i])):
                
                self.displacements[i][j] = self.bodies[j].pos - self.bodies[i].pos

        #loop through all bodies (twice)
        for i in range(len(self.bodies)):
            
            for j in range(len(self.bodies)):
                
                #if calculating values with the same bodies or with no difference in positions
                if i == j or norm(self.displacements[i][j]) == 0.0:
                    
                    self.bodies[j].acc[x] += [0.0, 0.0]
                
                else:
                    
                    #update acceleration
                    self.bodies[i].acc[x] += (self.bodies[j].mass / (norm(self.displacements[j][i]) ** 3)) * self.displacements[j][i]
            
            self.bodies[i].acc[x] = self.bodies[i].acc[x] * (-G)

        
    def stepforward(self):
        """
        Calculates the values of bodies for the next timestep in the simulation, using Beeman integration in lockstep.
        """
        #loop through all bodies and update their positions
        for i in range(len(self.bodies)):

            self.bodies[i].position()

        #update all accelerations
        self.acceleration(2)

        #loop through all bodies and update their velocities
        for i in range(len(self.bodies)):

            self.bodies[i].velocity()



    def cut(self):
        """
        Cuts the firt index of the position, velocity and acceleration arrays, and append a new blank array to have new values added to. 

        Returns
        -------
        None.

        """
        #erase the old values and append new array to maintain shape
        for i in range(len(self.bodies)):
            
            self.bodies[i].acc = np.concatenate((self.bodies[i].acc, [[0.0, 0.0]]), 0)
            self.bodies[i].acc = np.delete(self.bodies[i].acc, 0, 0)
          
    
    def update(self, x=None, o=None):
        """
        Updates the position of a index Body in Bodies.

        Parameters
        ----------
        x : int, optional
            Used for tageting a specific bodies attribute in list, default is None. If no value is given the function will not return anything.
        o : None, optional
            Used for determining whether or not to calculate orbital period data.

        Returns
        -------
        tuple, optional
            Position of Body in Bodies.

        """
        #update positions, velocities and accelerations of all bodies
        self.stepforward()

        #do this to update stuff
        self.cut()

        #useful incase this data is unwanted
        if o != None:

            self.orbitalPeriod()

        self.time += 1

        #if a new simulation has started this will clear the file where the energy data is stored
        if self.time <= 1:
            
            if self.directory != None:

                self.writeData(self.directory, "Total Energy", self.totalEnergy(), reset=True)

        else:
            
            #if a different file location is desired for the data to be stored
            if self.directory != None:

                self.writeData(self.directory, "Total Energy", self.totalEnergy())

            else:

                self.energyData.append(self.totalEnergy())

        if x != None:

            return (self.bodies[x].pos[0], self.bodies[x].pos[1])


    def animate(self, i):
        """
        Changes the center position of each Patch in patches by updating them with update().

        Parameters
        ----------
        i : int
            Used for itterating.

        Returns
        -------
        patches : list
            List of all patches generated for Simulation.

        """
        
        #loop through list of patches
        for x in range(len(self.patches)):
            
            #set center of patch to the updated position
            self.patches[x].center = self.update(x, o=1)
                   
        return self.patches
    

    def display(self):
        """
        Displays and runs the animation.

        Returns
        -------
        None.

        """
        self.patches = []

        #plotting begins!
        plt.style.use('dark_background')
        fig = plt.figure()
        ax = plt.axes()
        
        #loop for the number of bodies
        for i in range(len(self.bodies)):
            
            #append a cirlce with a position to the patches list
            self.patches.append(plt.Circle((self.bodies[i].pos[0], self.bodies[i].pos[0]), self.bodies[i].radius * 100, color=self.bodies[i].colour, label=self.bodies[i].name, animated=True))
        
        #loop through list of patches
        for i in range(len(self.patches)):
            
            #add patch to axis
            ax.add_patch(self.patches[i])
            
        #some formatting stuff
        ax.axis('scaled')

        #deteriminging appropriate scale based on orbits of Body
        radii = []
        for i in range(len(self.bodies)):

            radii.append(round(norm(self.bodies[i].initpos),4))

        scale = max(radii)

        ax.set_xlim(-1.3 * scale, 1.3 * scale)
        ax.set_ylim(-1.3 * scale, 1.3 * scale)
        
        #animate
        self.anim = FuncAnimation(fig, self.animate, frames=self.maxTime, repeat=self.repeat, interval = 1, blit=True)
        
        plt.title("Celestial Body Simulation")
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

        
    def orbitalPeriod(self):
        """
        Calculates the Orbital Period of a Body using Kepler's Third Law.

        Returns
        -------
        orbitalPeriods : list
            List of Orbital Periods of Body in Bodies.

        """
        secondsToEarthYears = 3.17098e-8

        #loop through bodies and check to see if they have crossed the x axis (this assumes they start on the x axis and did not instantly cross it)
        for i in range(len(self.bodies)):

            if self.bodies[i].prevPos[1] < 0 and self.bodies[i].pos[1] > 0 and self.orbitalPeriods[i] == 0.0:
                
                #calculate, print and store orbital periods of the bodies 
                oPeriod = self.time * self.timestep * secondsToEarthYears
                print(str(self.bodies[i].name) + " Orbital Period:   " + str(oPeriod) + " (Earth Years)")
                self.orbitalPeriods[i] = oPeriod

            else:

                continue
                    
        return self.orbitalPeriods    


    def totalKineticEnergy(self):
        """
        Calculates the Kinetic energy of all Bodies in Simulation and Sums them.

        Returns
        -------
        totalKineticEnergy : float
            Total Kinetic Energy of all Bodies in Simulation.

        """
        kineticEnergy = []

        #loop through bodies and calculate their kinetic energy
        for i in range(len(self.bodies)):
            
            kineticEnergy.append(0.5 * (self.bodies[i].mass * (norm(self.bodies[i].vel) ** 2)))

        #calculating the total energy of the system by summing the kinetic energy of all bodies
        totalKineticEnergy = np.sum(kineticEnergy)
        
        return totalKineticEnergy
    

    def totalPotenialEnergy(self):
        """
        Calculates the Potential energy of all Bodies in Simulation.

        Returns
        -------
        totalPotentialEnergy : float
            Total Potential Energy of all Bodies in Simulation.

        """
        G = Body.G
        totalPotentialEnergy = 0.0

        #loop through all the bodies (twice)
        for i in range(len(self.bodies)):

            for j in range(len(self.bodies)):
                
                #
                if i == j or norm(self.displacements[i][j]) == 0.0:

                    continue

                else:

                    totalPotentialEnergy += (G * self.bodies[i].mass * self.bodies[j].mass) / (norm(self.displacements[i][j]))

        #multiply sum by a factor of -0.5
        totalPotentialEnergy = totalPotentialEnergy * (-0.5)
        
        return totalPotentialEnergy
    
    
    def totalEnergy(self):
        """
        Returns the Sum of the total Kinetic and Potential Energies of all Bodies in Simulation.

        Returns
        -------
        E : float
            Total Combined Energy of all Bodies in Simulation.

        """
        E = 0.0

        k = self.totalKineticEnergy()
        u = self.totalPotenialEnergy()
        
        #sum of total kinetic and potential energy
        E = k + u
        
        return E


    def writeData(self, directory, name, value, reset=False):
        """
        Writes a value to a line in a file, can also reset the file to a blank file.

        Parameters
        ----------
        directory : String
            Folder Location.
        name : String
            File Name
        value : None
            Value to store in file
        reset : bool, optional
            Decide whether or not the information in the file is to be deleted, by default False

        Returns
        -------
        file : file
            File with stored Values
        """
        #access specified file and append value given to method
        with open(os.path.join(directory, name + ".txt"), 'a') as file:

            file.write(str(value) + "\n")
        
        #can be used to clear data in file
        if reset:

            with open(os.path.join(directory, name + ".txt"), 'a') as file:

                file.truncate(0)

        return file


    def energyGraph(self, folder=None, name="Total Energy.txt"):
        """
        Produces a graph over time of data read from a file, if None give the method will use the classes defualt locations.

        Parameters
        ----------
        folder : str, optional
            Folder Location for Data, by default None
        name : str, optional
            Name of File, by default "Total Energy.txt"
        """
        #if either the method or class instance has been given a file location and name
        if folder != None or self.directory != None:
            
            xValues = 0
            yData = []

            #if the method has been given a file location and name
            if folder != None:

                with open(os.path.join(folder, name), 'r') as file:
                    
                    #for line in file append values to list
                    for line in file:
                        
                        #strip any new line characters that might interupt evaluation of line
                        separate = line.rstrip("\n")
                        yData.append(float(eval(separate)))
                        xValues += 1

            #if the instance of class has been given a file location and name, but the method hasnt
            else:

                with open(os.path.join(self.directory, name), 'r') as file:
                    
                    #for line in file append values to list
                    for line in file:

                        #strip any new line characters that might interupt evaluation of line
                        separate = line.rstrip("\n")
                        yData.append(float(eval(separate)))
                        xValues += 1

            #create x values based on how many data point were read and multiply them by timestep of simulation for time
            xData = np.linspace(0, xValues, xValues)
            xData = xData * self.timestep

            #plot data
            plt.plot(xData, yData, "-g")
            plt.xlabel("Time")
            plt.ylabel("Total Energy")
            plt.title("Energy Graph")
            plt.show()

        else:
            
            #built in list incase no file was created to store data
            timeData = np.linspace(0, self.time, self.time)
            timeData = timeData * self.timestep
            energyData = self.energyData

            #plot data
            plt.plot(timeData, energyData, "-g")
            plt.xlabel("Time")
            plt.ylabel("Total Energy")
            plt.title("Energy Graph")
            plt.show()


    def addBody(self, *newBodies):
        """
        Adds one or more Body objects to the called upon BodySim object.

        Returns
        -------
        list
            List of the called upon objects Body objects after appedning to BodySim.bodies.

        Raises
        ------
        TypeError
            Error raised when a non Body type object is added to BodySim.bodies.
        """
        #loop through the arguement(s) given to method (*args)
        for body in newBodies:

            #check if variable is of type Body
            if isinstance(body, Body):

                #append Body to list of bodies in class instance
                self.bodies.append(body)
                self.orbitalPeriods.append(0.0)
                self.resetValues()
                self.acceleration(1)

            else:
                
                #raised if one of the arguements was not a Body object
                raise TypeError("Type \"" + str(type(body)) + "\", not of type \"Body\"")

        return self.bodies


    def removeBody(self, *delBodies):
        """
        Removes one or more Body objects from the BodySim object this function is called on.

        Returns
        -------
        list
            List of the called upon objects Body objects after appedning to BodySim.bodies.
        """
        bodiesRemoved = 0
        lenDelBodies = 0

        #loop through the arguement(s) given to method (*args)
        for body in delBodies:

            lenDelBodies += 1

            #check if variable is of type Body
            if isinstance(body, Body):

                #loop through list of bodies
                for i in range(len(self.bodies)):
                    
                    #check to see if any of the objects in list match the one to be removed
                    if self.bodies[i].name == body.name:

                        #remove item
                        del self.bodies[i]
                        del self.orbitalPeriods[i]
                        bodiesRemoved += 1

                    else:

                        continue

            else:

                #raised if one of the arguements was not a Body object
                raise TypeError("Type \"" + str(type(body)) + "\", not of type \"Body\"")
            
        #just some nice feedback on how many items were removed from the instance
        if bodiesRemoved != lenDelBodies:

            print(str(bodiesRemoved) + " Body Objects removed, " + str(lenDelBodies - bodiesRemoved) + " item(s) not found")

        return self.bodies


    def proximity(self, target1Name, target2Name):
        """
        Calculates the absolute distance between two Body object in the BodySim object.

        Parameters
        ----------
        target1Name : str
            Name of first target Body.
        target2Name : str
            Name of second target Body

        Returns
        -------
        distance : float
            Distance between the two target Body objects
        """
        #loop through list of bodies
        for i in range(len(self.bodies)):

            #check to see if the body is the same using the "name" variable
            if self.bodies[i].name == target1Name:

                target1 = i

            else:
                continue
                raise NameError("Target 1 not found")

        for i in range(len(self.bodies)):

            if self.bodies[i].name == target2Name:

                target2 = i

            else:
                continue
                raise NameError("Target 2 not found")

        distance = abs(norm(self.displacements[target1][target2]))

        return distance

    def resetValues(self):
        """
        Reset the positions of all Body objects in bodies using their initial values from file

        Returns
        -------
        list
            List of Body ojects
        """
        for i in range(len(self.bodies)):

            self.bodies[i].resetValues()

        return self.bodies