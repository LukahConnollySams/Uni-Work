#imports
from array import array
from typing import List
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import random as R
import math as m


class Ising_Model(object):

    def __init__(self, J, T, N, random_lattice=True, sweeps=20, method='kawasaki'):

        self.J = J
        self.T = T
        self.N = N
        self.k_B = 1 
        self.sweeps = sweeps
        self.updates = 0
        self.sus = []
        self.cap = []
        self.m = []
        self.e = []
        self.Ts = []
        self.method = method

        if type(random_lattice) == type(True):
            
            self.lattice = np.random.choice([-1, 1], size = (self.N, self.N))

        else:
            
            self.lattice = random_lattice
        
    def __str__(self) -> str:
        """
        Prints IsingModel.
        """
        return f"Lattice size: {self.N}\n Method: {self.method}\n Random Lattice: {self.random_lattice}\n Lattice: {self.lattice}"

    def too_close(self, xy) -> bool:
        """
        Determines if 2 points, given by input array, in the lattice are neighbours.

        Parameters
        ----------
        xy : list/array
            x and y postions of 2 points in lattice.

        Returns
        -------
        bool
            Returns True if they are neighbours.
        """

        #find the 2 points and determine the difference in the indicies
        p1, p2 = np.array([xy[0], xy[1]]), np.array([xy[2], xy[3]])
        p_diff = p2 - p1

        #if they are next to each other in lattice
        if abs(p_diff[0]) <= 1 and abs(p_diff[1]) <= 1 and abs(p_diff[0]) + abs(p_diff[1]) < 2:

            return True

        #edge case for end or start of lattice
        elif abs(p_diff[0]) == len(self.lattice[0]) or abs(p_diff[1]) == len(self.lattice):

            return True

        #if they arent next to each other in lattice
        else:

            return False

    def same_spin(self, xy) -> bool:
        """
        Determines if 2 points, given by input array, in the lattice are of the same spin.

        Parameters
        ----------
        xy : list/array
            x and y postions of 2 points in lattice.

        Returns
        -------
        bool
            Returns True if they are of the same spin.
        """

        #if the spin of the two points is the same
        if self.lattice[xy[0]][xy[1]] == self.lattice[xy[2]][xy[3]]:

            return True
        
        #if the two points arent the same
        else:

            return False

    def random_p(self) -> array:
        """
        Generates random indices for either 1 or 2 points on a 2D array.

        Returns
        -------
        array
            Array containing indices for 2D array.
 
        """

        #chose 4 random numbers to represent two random points x and y in array
        if self.method == 'kawasaki':

            xy = np.random.choice(self.N, 4)

            #redo random indices if they are neighbours or are of the same spin
            while self.too_close(xy) or self.same_spin(xy):

                xy = np.random.choice(self.N, 4)
                
        #chose 2 random numbers to represent 1 random point's x and y in array
        elif self.method == 'glaubor':

            xy = np.random.choice(self.N, 2) 


        return xy

    def neighbour_sum_products(self, i, j) -> int:
        """
        Determines the sum of the products of the home cell, and the neighbouring cells in a 2D array.

        Parameters
        ----------
        i : int
            First indices of the target cell in array (x coordinate).
        j : int
            Second indices of the target cell in array (y coordinate).

        Returns
        -------
        int/float
            Sum of the products of the home cell and neighbouring cells.
        """

        #spin value of target and list for value appending
        home = self.lattice[i][j]
        neighbours = []

        #append neighbour's spin to list
        neighbours.append(np.roll(self.lattice, 1, axis=1)[i][j])
        neighbours.append(np.roll(self.lattice, -1, axis=1)[i][j])
        neighbours.append(np.roll(self.lattice, 1, axis = 0 )[i][j])
        neighbours.append(np.roll(self.lattice, -1, axis = 0)[i][j])

        #sum the product of neighbour and target
        sum_values = np.sum(home * np.array(neighbours))


        return sum_values

    def energy(self, spin_sum) -> float:
        """
        Determines the Energy of a configuration given the sum of its spin products.

        Parameters
        ----------
        spin_sum : float/int
            sum of the products of the configurations spin.

        Returns
        -------
        float/int
            Energy of configuration.
        """

        #energy of a given spin configuration
        E = -1 * self.J * spin_sum

        return E
    
    def probability_s(self, E) -> float:
        """
        Determine the probability of a certain configuration of spins in a lattice.

        Parameters
        ----------
        E : float/int
            Energy of a configuration.

        Returns
        -------
        float
            Probability of a given configuration.
        """

        #probability of a given energy
        probability = m.exp(-1 * E / (self.k_B * self.T))

        return probability

    def change_spin(self, xy) -> List[int]:
        """
        Changes spin of particle(s) in lattice.

        Parameters
        ----------
        xy : list, optional
           list of lattice position(s), by default None

        Returns
        -------
        list/array
            list of lattice position(s)

        Raises
        ------
        ValueError
            xy is of invalid size for function
        """

        #if 2 separate points ('kawasaki')
        if len(xy) == 4:

            #find indices
            i = xy[0]
            j = xy[1]
            i2 = xy[2]
            j2 = xy[3]

            #swap spins of the 2 points
            self.lattice[i][j] = -1 * self.lattice[i][j]
            self.lattice[i2][j2] = -1 * self.lattice[i2][j2]


        #if only one point ('glaubor')
        elif len(xy) == 2:
            
            #find indices
            i = xy[0]
            j = xy[1]

            #invert spin
            self.lattice[i][j] = self.lattice[i][j] * (-1)

        #when xy isnt the correct size
        else:

            raise ValueError('xy size did not match, list must be of size 2 or 4')


        return xy

    def total_magnetism(self) -> float:
        """
        Caclulates average of M^2 and M

        Returns
        -------
        float
            Values for M^2 and M
        """
        #lists for appending 
        m2 = []
        m = []

        #looping constant
        x = 0

        #loop for a set amount of samples
        while x <= self.sweeps/self.N:

            #select a random point
            xy = self.random_p()
            m2.append(self.lattice[xy[0]][xy[1]]**2)
            m.append(self.lattice[xy[0]][xy[1]])

            #increment for while loop
            x += 1

        #sum and change for average squared as well as the average of the square
        m2 = np.sum(m2) / len(m2)
        m_ = np.sum(m) / len(m)

        return m2, m_

    def susceptibility(self, m2, m) -> float:
        """
        Calculates susceptibility as X = (<m^2> - <m>^2) / (N * k_B * T)

        Parameters
        ----------
        m2 : float/int
            Average of Magnetisation Squared
        m : float/int
            Average Magnetisation

        Returns
        -------
        float/int
            Susceptibility of data.
        """
        #calculates susceptibility from input values
        x = (m2 - m**2) / (self.N * self.k_B * self.T)

        return x

    def total_energy(self) -> float:
        """
        Caclulates average of E^2 and E

        Returns
        -------
        float
            Values for E^2 and E
        """
        #used to store data from loop
        e2 = []
        e_ = []
        x = 0

        #loop for set amount of datapoints for sample
        while x <= self.sweeps/self.N:

            #select random point
            xy = self.random_p()
            #append values to list
            e_sum = self.neighbour_sum_products(xy[0], xy[1])
            e = self.energy(e_sum)
            e2.append(e**2)
            e_.append(e)

            x += 1

        #average the values
        e2 = np.sum(e2) / len(e2)
        e_ = np.sum(e_) / len(e_)

        return e2, e_

    def scaled_heat_capacity(self, e2, e):
        """
        Caluclates the Scaled Heat Capacity as C = (<e^2> - <e>^2) / (N * k_B * T**2)

        Parameters
        ----------
        e2 : float/int
            Average of Energy Squared
        e : float/int
            Average Energy

        Returns
        -------
        float/int
            Scaled Heat Capacity of Data
        """
        #calculates Scaled Heat Capacity from given values
        c = (e2 - e**2) / (self.N * self.k_B * self.T**2)

        return c

    def update(self, i, sweeps_data=False):
        """
        Updates the lattice depending on objects method.

        Parameters
        ----------
        i : int
            Itterable used for animation.
        sweeps_data : bool, optional
            Will calculate extra data and store in object, by default False

        Returns
        -------
        array
            Updated Lattice of object.
        """
        #used to itterate while loop
        x = 0

        #loop for number of sweeps
        while x <= self.sweeps:

            #kawasaki
            if self.method == 'kawasaki':
                
                #get random points
                xy = self.random_p()

                #calculate spin configuration sum
                start_sum1 = self.neighbour_sum_products(xy[0], xy[1])
                start_sum2 = self.neighbour_sum_products(xy[2], xy[3])

                #calculate energy from spin configuration sum
                start_e1 = self.energy(start_sum1)
                start_e2 = self.energy(start_sum2)

                #find total energy of the two points
                start_total_e = start_e1 + start_e2

                #swap spins
                self.change_spin(xy)

                #calculate spin configuration sum
                end_sum1 = self.neighbour_sum_products(xy[0], xy[1])
                end_sum2 = self.neighbour_sum_products(xy[2], xy[3])

                #calculate energy from spin configuration sum
                end_e1 = self.energy(end_sum1)
                end_e2 = self.energy(end_sum2)

                #find total energy of the two points
                end_total_e = end_e1 + end_e2

                #calculate difference in energy and the probability of the swap happening
                e_total_diff = end_total_e - start_total_e
                probability = self.probability_s(e_total_diff)

                #if probability is less than one, the swap happens based on its probability
                if 1 > probability:

                    #generate random value between 0 and 1 to compare with probability
                    ran_value = R.random()

                    if ran_value > probability:
                        
                        #revert spin swap
                        self.change_spin(xy)
                    

            #glaubor
            elif self.method == 'glaubor':

                #get random point
                xy = self.random_p()

                #calculate spin configuration sum
                start_sum = self.neighbour_sum_products(xy[0], xy[1])

                #calculate energy from spin configuration sum
                start_e = self.energy(start_sum)

                #flip spins
                self.change_spin(xy)

                #calculate spin configuration sum
                end_sum = self.neighbour_sum_products(xy[0], xy[1])

                #calculate energy from spin configuration sum
                end_e = self.energy(end_sum)

                #calculate difference in energy and the probability of the swap happening
                e_diff = end_e - start_e
                probability = self.probability_s(e_diff)

                #if probability is less than one, the swap happens based on its probability
                if 1 > probability:
                    
                    #generate random value between 0 and 1 to compare with probability
                    ran_value = R.random()

                    if ran_value > probability:
                        
                        #revert spin flip
                        self.change_spin(xy)

            self.updates += 1

            #self.temperature_plots_values(x)

            #increment for while loop
            x += 1

            if sweeps_data:

                #calculates other data
                self.temperature_plots_values(x)

            else:

                continue

        return self.lattice       

    def temperature_plots_values(self, x):
        """
        Calculates extra data, which is stored in the object

        Parameters
        ----------
        x : int
            Number of current sweeps in step.
        """
        
        #after and equilibrium has been allowed to form
        if self.updates >= 100:
            
            #every 10 sweeps in an update
            if x >= 10:

                #calculate averages of values
                m2, m = self.total_magnetism()
                e2, e = self.total_energy()

                #calculate susceptibility and scaled heat capacity
                sus = self.susceptibility(m2, m)
                cap = self.scaled_heat_capacity(e2, e)
                
                #append to object variables
                self.m.append(m)
                self.e.append(e)
                self.sus.append(sus)
                self.cap.append(cap)
                self.Ts.append(self.T)

            else:

                pass
        
        else:

            pass

    def animate(self, i):
        """
        Returns plot in list length 1 for animation use.

        Parameters
        ----------
        i : int
            Iterrable for animation.

        Returns
        -------
        list
            List of length 1, containing imshow plot from lattice.
        """

        #store imshow() plot of lattice 
        im = plt.imshow(self.update(i), interpolation='none', origin='lower')

        #return plot as list of length 1 for FuncAnimation
        return [im]

    def display(self, interval, frames):
        """
        Animates and displays evolution of lattice under the object's simulation method

        Parameters
        ----------
        interval : int
            Times in ms between plot/animation updates.
        frames : int
            Number of repeats
        """
        #create plt figure
        fig = plt.figure()
        
        #animate lattice updates
        self.anim = FuncAnimation(fig, self.animate, frames=frames, repeat=False, interval=interval, blit=True)
        
        #show plt figure and animation
        plt.show()


def main():

    #read parameters from terminal input
    print('Lattice Size: ')
    N = int(input())

    print('Method: ')
    method = str(input()).lower()

    #forces user to input correct method/spelling
    while True:

        if method == 'glaubor' or method == 'kawasaki':

            break

        else:

            print('Method input not compatible, try again: ')
            method = str(input()).lower()

    print('Sweep Size: ')
    sweeps = int(input())

    print('Simulation Length: ')
    length = int(input())

    print('Random Lattice (y/n): ')
    lattice_type = str(input())

    #avoid incorrect inputs, and input custom lattices (mainly here for fun)
    try:

        if lattice_type == 'y':

            random_lattice = True

        elif lattice_type == 'n':

            print('Input 1D list here  as values separated by a space (will convert to np array with shape given above):')
            random_lattice = str(input())
            random_lattice = random_lattice.split()
            for ind, val in enumerate(random_lattice):
                
                val = int(val)

            random_lattice = np.array(random_lattice).reshape((N, N))

    #default to True if input is incompatible
    except:

        print('Random Lattice has defaulted to "y"')
        random_lattice = True

    #Create IsingModel Object from input parameters
    x = Ising_Model(1, 1, N, random_lattice=random_lattice, sweeps=sweeps, method=method)

    #create one for animation demonstration
    print('Would you like an example animation? (y/n): ')

    example = str(input())

    if example == 'y':


        y = Ising_Model(1, 1, N, sweeps=sweeps, method=method)
        y.display(1, length)

    elif example == 'n':

        pass 

    #used to store data for later reading/writing
    sus_data = []
    cap_data = []
    t_data = []
    m_data = []
    e_data = []
    temps = np.arange(1, 3, step=0.1)

    #loops through given range of temperatures
    for i in range(len(temps)):

        #Create IsingModel Object from input parameters
        x = Ising_Model(1, temps[i], N, random_lattice=random_lattice, sweeps=sweeps, method=method)
        j = 0
        
        #update for given number of times
        while j < length:

            #Run simulation
            x.update(0, sweeps_data=True)
            j += 1

        ##append data to respective lists
        sus_data.append(np.average(x.sus))
        cap_data.append(np.average(x.cap))
        t_data.append(temps[i])
        m_data.append(np.average(x.m))
        e_data.append(np.average(x.e))

        #print comment too make sure youre not going insane from the wait (can be commented out)
        print(f"temperature: {temps[i]} done.")

    #plots of data
    plt.plot(t_data, sus_data, label='Susceptibility')
    plt.title('Susceptibility vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('Susceptibility')
    plt.show()
    plt.plot(t_data, cap_data, label='Scaled Heat Capacity')
    plt.title('Scaled Heat Capacity vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('Scaled Heat Capacity')
    plt.show()
    plt.plot(t_data, m_data, label='Total Magnetisation')
    plt.title('Total Magnetisation vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('Total Magnetisation')
    plt.show()
    plt.plot(t_data, e_data, label='Total Energy')
    plt.title('Total Energy vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('Total Energy')
    plt.show()
    

    #write data to file
    file = open('ising_model_data.txt', 'w')
    #file header
    file.write('#Temp   Susceptibility  SH_Capacity    Magnetisation    Energy\n')

    #each line is written in f string format
    for i in range(len(t_data)):

        file.write(f'{round(t_data[i], 2)} {sus_data[i]} {cap_data[i]} {m_data[i]} {e_data[i]}\n')


#Block imports from running main func
if __name__ == '__main__':
    main()