# -*- coding: utf-8 -*-
"""
Created on Thurs Mar 31 14:26:23 2022

@author: Lukah
"""
from multiprocessing.sharedctypes import Value
from Body import Body
from BodySim import BodySim
import numpy as np
import math
from numpy.linalg import norm
#from numba import jit

class LaunchSim(BodySim):
    """
    Launch simulation of Body object, similar to BodySim object.

    Parameters
    ----------
    BodySim : object
        BodySim like object

    Raises
    ------
    ValueError
        Raised when the velRange parameter is given too large a value.
    ValueError
        Raised when the thetaRange parameter is given values outside of acceptable ranges
    """
    G = Body.G

    def __init__(self, bodies, launchObject, timestep, velRange, thetaRange, maxTime=100, repeat=False, directory=None):
        """_summary_

        Parameters
        ----------
        bodies : list
            List of Body objects
        launchObject : Body
            Body that is being launched
        timestep : float, int
            Timestep of simulation (time between updates, in seconds)
        velRange : list, tuple
            Range of velocities tested by simulation
        thetaRange : list, tuple
            Range of angles tested by simulation
        maxTime : int, optional
            Maximum number of updates run by simulation, by default 100
        repeat : bool, optional
            Whether or not the simulation has a maximum number of itterations/updates, or keeps going until proccess is stopped, by default False
        directory : _type_, optional
            Folder name for storing data, by default None

        Raises
        ------
        ValueError
            Raised when the input of the launch velocity range is greater than the speed of light
        ValueError
            Raised when the input of the launch angle range is greater 2 pi or less than 0
        """
        #initialise the same as BodySim object
        super(LaunchSim, self).__init__(bodies, timestep, maxTime, repeat, directory)
        #add the launch object to list of bodies
        super(LaunchSim, self).addBody(launchObject)

        self.launchObject = launchObject

        #velocity cant be greater than the speed of light
        if abs(velRange[1]) > 3 * 10 ** 8:

            #error is raised when velocity is too high
            raise ValueError("Value too large, must be smaller than 3e8")
        
        else:
            
            self.velRange = velRange

        #angle range must be within 0 and 2 pi, though techniaclly doesnt matter
        if 0 > thetaRange[0] > 2 * math.pi or 0 > thetaRange[1] > 2 * math.pi:

            #error raise when not within accepted range
            raise ValueError("Values must be within 0 and 2 pi.")

        else:

            self.thetaRange = thetaRange


    def update(self):
        """
        Uses the super() function to return the update() method of the BodySim class on this object.

        Returns
        -------
        method
            BodySim.update()
        """
        return super(LaunchSim, self).update()

    def animate(self, i):
        """
        Uses the super() function to return the animate() method of teh BodySim class on this object.

        Parameters
        ----------
        i : int
            Used for itterating.

        Returns
        -------
        patches : list
            List of all patches generated for Simulation.

        """
        return super().animate(i)


    def display(self):
        """
        Uses the super() function to return the display() method of the BodySim class on this object.

        Returns
        -------
        method
            BodySim.display()
        """
        return super(LaunchSim, self).display()


    def proximity(self, target1Name, target2Name):
        """
        Uses the super() function to return the proximity() method of the BodySim class on this object.

        Parameters
        ----------
        target1Name : str
            Name of a Body object in self.bodies
        target2Name : str
            Name of a Body object in self.bodies

        Returns
        -------
        float
            Distancec between two Body objects in simulation.
        """
        return super(LaunchSim, self).proximity(target1Name, target2Name)

    def resetValues(self):
        """
        Uses the super() function to return the restValues() method of the BodySim class on this object.

        Returns
        -------
        list
            List of body objects.
        """
        return super().resetValues()


    def launchTo(self, targetName, returnTarget, dv=10, dtheta=10):
        """
        Simulates a body simulation where the proximity of the launch object to the target and return target is monitored and the closest value found. 

        Parameters
        ----------
        targetName : str
            Name of target Body
        returnTarget : str
            Name of return target Body
        dv : int, optional
            Number of values between range of angles, by default 10
        dtheta : int, optional
            number of values between range of velocities, by default 10

        Returns
        -------
        float
            Speed at which the object gets the closest to the target.
        """
        #variables used to store the idex of certain items in list, used later
        k = 0
        l = 0
        m = 0

        #loop through bodies to find specific bodies using "name" attribute
        for i in range(len(self.bodies)):
            
            #if body matches target name
            if targetName == self.bodies[i].name:

                #store index
                k = i
            
            #if body matches launch object's name
            elif self.launchObject.name == self.bodies[i].name:

                #store index
                l = i

            #if body matches the return target'#s name
            elif returnTarget == self.bodies[i].name:

                #store index
                m = i

            #if none match
            else:

                continue
        
        #initaialise the range of speeds being searched through
        minValue = self.velRange[0]
        maxValue = self.velRange[1]
        VRange = np.linspace(minValue, maxValue, dv)

        #initialise list for range of angles to search through
        minAngle = self.thetaRange[0]
        maxAngle = self.thetaRange[1]
        theta = np.linspace(minAngle, maxAngle , dtheta)

        #list thats will be used later
        allproxy1 = []
        allproxy2 = []

        #search for v and take two lowest prox then repeat until within certain sig figures?????????
        
        #this will take some time (mainly depending on how many values are being tested)
        print("This will take some time:")

        #loop through angles
        for i in range(len(theta)):

            #loop through magnitude of velocities
            for j in range(len(VRange)):

                #reset positions, speeds and accelerations of bodies in simulations
                self.resetValues()

                #set velocity values for the launch object
                self.bodies[l].vel[0] = VRange[j] * math.cos(theta[i])
                self.bodies[l].vel[1] = VRange[j] * math.sin(theta[i]) 

                #print this velocity
                print("Speed " + str(self.bodies[l].vel) + " Norm " + str(norm(self.bodies[l].vel)))

                #used to store the proximities of object to the launch object
                proxy1 = []
                x = 0

                #loop runs until simulation is done
                while x < self.maxTime:
                    
                    #update simulation for timestep
                    self.update()

                    #store distances between objects
                    d1 = self.proximity(self.bodies[k].name, self.bodies[l].name)
                    proxy1.append(d1)

                    #after at least two simulations
                    if x >= 1:

                        #if the objects were closer than last timstep
                        if proxy1[1] < proxy1[0]:

                            #delete the last timestep value
                            del proxy1[0]

                        #if the objects were closer last timestep
                        elif proxy1[0] < proxy1[1]:

                            #delete current timestep's value
                            del proxy1[1]

                        #if the distances were equal
                        elif proxy1[0] == proxy1[1]: 

                            #delete last timestep's value
                            del proxy1[0]

                        #if the current timesteps value is zero
                        elif proxy1[1] == 0.0:

                            #delete last timesteps value
                            del proxy1[0]

                            #break loop as the distance should no longer change 
                            break
                    
                    x += 1

                #print to tell user that one value has been completely tested
                print("Velocity simulated " + str((i) * dtheta + (j+1)))
                
                #append the closest distance to list
                allproxy1.append(proxy1[0])

                #when more than one velocity has been tested
                if len(allproxy1) > 1:
                    
                    #if closest value was closer than previous simulated velocity 
                    if allproxy1[1] < allproxy1[0]:

                        #delete previous value
                        del allproxy1[0]
                        #store velocity
                        finalVel = self.bodies[l].vel

                    #if closest value was previous simulated velocity
                    elif allproxy1[0] < allproxy1[1]:

                        #delete current value
                        del allproxy1[1]
                        
                    #if both values are the same
                    elif allproxy1[0] == allproxy1[1]: 
                        
                        #delete previous value
                        del allproxy1[0]
                        
                #print closest distance value
                print("New closest distance " + str(allproxy1))

        #set velocity of launch object to best velocity for closest proximity to target
        self.bodies[l].vel = finalVel
        t = 0
        proxy2 = []

        #resest values for next simulation
        self.resetValues()

        #run simulation again and find closest distance to return target??????????? need to specify time???????_______________________
        while t < self.maxTime:
                    
            #update simulation for timestep
            self.update()

            #store distances between objects
            d2 = self.proximity(self.bodies[m].name, self.bodies[l].name)
            proxy2.append(d2)

            #after at least two simulations
            if t >= 1:

                #if the objects were closer than last timstep
                if proxy2[1] < proxy2[0]:

                    #delete the last timestep value
                    del proxy2[0]

                #if the objects were closer last timestep
                elif proxy2[0] < proxy2[1]:

                    #delete current timestep's value
                    del proxy2[1]

                #if the distances were equal
                elif proxy2[0] == proxy2[1]: 

                    #delete last timestep's value
                    del proxy2[0]

                #if the current timesteps value is zero
                elif proxy2[1] == 0.0:

                    #delete last timesteps value
                    del proxy2[0]

                    #break loop as the distance should no longer change 
                    break
            
            t += 1

        #print output
        print("Speed: ")
        print(finalVel)
        print("Closest Distance to target: ")
        print(allproxy1)
        print("Closest Distance to return target: ")
        print(proxy2)
        
        #return the velocity that got the closest to the target
        return finalVel