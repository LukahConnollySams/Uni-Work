# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 17:49:17 2022

@author: Lukah
"""
import numpy as np
import os

class Body(object):
    """
    Object that represents a celestial Body.

    Parameters
    ----------
    object : object
        Class inheritence
    """
    G = 6.67 * 10**(-11)
    
    def __init__(self, fileDirectory, timestep):
        """
        Initialises Body object, taking a file location and reading the data for the instanced object variables.

        Parameters
        ----------
        fileDirectory : str
            location of file with information for object.
        timestep : float, int
            Timestep used for integration.
        """
        self.file = fileDirectory
        name, mass, pos, vel, colour, radius = self.readFile(self.file)
        self.name = name
        self.mass = mass
        self.pos = np.array(pos)
        self.prevPos = np.array([0.0, 0.0])
        self.initpos = np.array(pos)
        self.vel = np.array(vel)
        self.initvel = np.array(vel)
        self.acc = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
        self.colour = colour
        self.radius = radius
        self.timestep = timestep
    
    
    def readFile(self, fileIn):
        """
        Reads files and extracts values needed for object Body. 

        Parameters
        ----------
        fileIn : File
            A File with information of am object.

        Returns
        -------
        name : String
            Name of the Body.
        mass : float
            Float for the mass of the Body in kg.
        initpos : tuple
            Tuple for the coordinates of the Body for the simulation.
        initvel : tuple
            Tuple for the components of the Bodies velocity for the simulation.
        colour : str
            Colour of object, as a string for simulation. Anything that can be used for matplotlib colour can be accepted.
        radius : float
            Radius of object for simulation.

        """
        G = Body.G
        file = open(fileIn, "r")
        
        #loop through each line of file
        for line in file:
            
            #extract name variable
            if line.startswith("Name"):
                
                separate = line.split(" = ")
                #separate[1].rstrip(".txt")
                name = separate[1].strip("\n")
            
            #extract mass value
            elif line.startswith("mass"):
                
                separate = line.split(" = ")
                mass = float(eval(separate[1]))
                
            #extract initial position value as tuple
            elif line.startswith("initial position"):
                
                separate = line.split(" = ")

                #separate string tuple
                values = separate[1].split(",")

                #strip brackets, evaluate values and store in tuple
                one = values[0].lstrip("(")
                two = values[1].strip(")\n")
                initpos = (float(eval(one)), float(two))
                
            #extract initial velocity value as tuple
            elif line.startswith("initial velocity"):
                
                separate = line.split(" = ")

                #separate string tuple
                values = separate[1].split(",")

                #strip brackets, evaluate values and store in tuple
                one = values[0].lstrip("(")
                two = values[1].strip(")\n")
                initvel = (float(eval(one)), float(eval(two)))
                
            #extract colour value
            elif line.startswith("colour"):
                
                separate = line.split(" = ")
                two = separate[1].strip("\n")
                colour = str(two)
                
            #extract radius value
            elif line.startswith("radius"):
                
                separate = line.split(" = ")
                two = separate[1].strip("\n")
                radius = float(eval(two))
                
            else:
                
                #ignore other lines (e.g. blank lines or headers)
                continue
            
        return name, mass, initpos, initvel, colour, radius
    

    def position(self):
        """
        Calculates next position in time using Beeman integration.

        Returns
        -------
        np.array
            Numpy array containing x and y values for object's position. 
        """
        #position formula   
        self.prevPos = self.pos 
        self.pos = self.pos + self.vel * self.timestep + (1/6) * (4 * self.acc[1] - self.acc[0]) * self.timestep ** 2
        
        return self.pos

    def velocity(self):
        """
        Calculates next velocity in time using Beeman integration

        Returns
        -------
        np.array
            Numpy array containing x and y values for object's velocity
        """
        self.vel = self.vel + (1/6) * (2 * self.acc[2] + 5 * self.acc[1] - self.acc[0]) * self.timestep

        return self.vel
        
    
    def kineticEnergy(self):
        """
        Calculates the kinetic energy of a Body

        Returns
        -------
        float
            Kinetic energy of the Body in joules.

        """
        self.kineticEnergy = 0.5 * self.mass * (self.vel ** 2)
        
        return self.kineticEnergy


    def writeFile(directory, name, mass, pos, vel, colour, radius):
        """
        Write file in format for Body object.

        Parameters
        ----------
        directory : str
            String of the new files location
        name : str
            Name of new Body object file
        mass : float, str
            Mass of new Body object file
        pos : tuple, str
            Position of new Body object file
        vel : tuple, str
            Velocity of new Body object file
        colour : str
            Colour of new Body object file
        radius : float, str
            Radius of new Body object file

        Returns
        -------
        str_
            Name of file created.
        """
        #open file and write values in format for storing data ready to be read by this class
        with open(os.path.join(directory,name + ".txt"), 'w') as file:

            line1 = "Name = " + str(name) + "\n"
            line2 = "" + "\n"
            line3 = "mass = " + str(mass) + "\n"
            line4 = "initial position = " + str(pos) + "\n"
            line5 = "initial velocity = " + str(vel) + "\n"
            line6 = "colour = " + str(colour) + "\n"
            line7 = "radius = " + str(radius)

            file.writelines([line1, line2, line3, line4, line5, line6, line7])

        return name + ".txt"

    
    def __str__(self):
        """
        Format for representing/converting this object to a string.

        Returns
        -------
        stringout: str
            String representation of Body Object.
        """
        #essentially printing the first 8 lines of the file
        line0 = "\n\"\"\"\nBody Object: " + "\n"
        line1 = "   Name: " + str(self.name) + "\n"
        line2 = "   " + "\n"
        line3 = "   Mass: " + str(self.mass) + "\n"
        line4 = "   Position: " + str(self.pos) + "\n"
        line5 = "   Initial Position: " + str(self.initpos) + "\n"
        line6 = "   Velocity: " + str(self.vel) + "\n"
        line7 = "   Initial Velocity: " + str(self.initvel) + "\n"
        line8 = "   Acceleration: " + str(self.acc[1]) + "\n"
        line9 = "   Colour: " + str(self.colour) + "\n"
        line10 = "   Radius: " + str(self.radius) + "\n\"\"\"\n"

        stringout = line0 + line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 +line9 + line10

        return stringout

    def resetValues(self):
        """
        Resest the position, velocity and acceleration values of the body object beack to the inital values.

        Returns
        -------
        Body
            This Body object instance.
        """
        self.pos = self.initpos
        self.vel = self.initvel
        self.acc = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])

        return self