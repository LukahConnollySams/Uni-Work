from array import array
from typing import Tuple
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from timeit import default_timer as timer
import random as R

class GameOfLife(object):

    def __init__(self, size: int, sweeps: int, n: int=2, prob: Tuple[float, float]=(0.2, 0.3), lattice=None) -> None:

        self.size = size
        self.sweeps = sweeps
        self.n = n
        self.p1 = prob[0]
        self.p2 = prob[1]

        #feature to input custom initial lattice
        if lattice is not None:

            self.lattice = lattice

        else:

            #creates size x size grid with random 0's or 1's
            self.lattice = np.random.choice([0, 1], size=(self.size, self.size))

        #store number of alive cells at initialisation
        self.alive = []
        self.alive.append(self.lattice.flatten().tolist().count(1.0))

    def __str__(self) -> str:
        """
        String of object

        Returns
        -------
        str
            String of object
        """
        return f'{self.alive[-1]} Cells Alive \nLattice: \n{self.lattice}'

    def sum_neighbours(self) -> array:
        """
        Creates an array of same shape as lattice, but with value equal to number of alive neighbours

        Returns
        -------
        array
            number of neighbours of each cells in lattice
        """
        neighbours = []

        #determine values, and append, 8 neighbouring cells
        neighbours.append(np.roll(self.lattice, 1, axis=1))
        neighbours.append(np.roll(self.lattice, -1, axis=1))
        neighbours.append(np.roll(self.lattice, 1, axis=0 ))
        neighbours.append(np.roll(self.lattice, -1, axis=0))
        neighbours.append(np.roll(self.lattice, (1, 1), axis=(1, 0)))
        neighbours.append(np.roll(self.lattice, (1, -1), axis=(1, 0)))
        neighbours.append(np.roll(self.lattice, (-1, 1), axis=(1, 0)))
        neighbours.append(np.roll(self.lattice, (-1, -1), axis=(1, 0)))

        #sum the rolled lattices on top of each other to create lattice shaped array with number of alive nieghbours in each cell
        return np.sum(neighbours, axis=0)

    def random_p(self):
        """
        Generates random indices for 1 point on a 2D array.

        Returns
        -------
        array
            Array containing indices for 2D array.
        """

        #random 2 values from the idices of the list
        return  np.random.choice(self.size, 2)

    def update(self, x) -> array:
        """
        Updates the entire lattice of the game of life

        Parameters
        ----------
        x : int
            itterable used for animation

        Returns
        -------
        array
            Updated lattice
        """

        # used for later calculations
        updated_alive = 0
        #calculate how many neighbouring cells are alive
        sum = self.sum_neighbours()

        #loop through number of sweeps before pushing update
        m = 0
        while m <= self.sweeps:

            #select random point
            xy = self.random_p()
            i, j = xy[0], xy[1]

            #if cell is dead
            if self.lattice[i][j] == 0:

                #make cell alive if n neighbours are also
                if sum[i][j] == self.n:
                    
                    #generate random number and compare to probability
                    if self.p2 > R.random():
                        self.lattice[i][j] = 1 

                        #add to current alive cell count
                        updated_alive +=1
                    
                    #kill cell
                    else:
                        
                        self.lattice[i][j] = 0

            #if cell is alive, kill
            else:
                
                #generate random number and compare to probability
                if self.p1 > R.random():

                    self.lattice[i][j] = 0

                #kill cell
                else:

                    updated_alive +=1

            m += 1


        #append number of alive cells after update
        self.alive.append(updated_alive)


        return self.lattice

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
        im = plt.imshow(self.update(i), interpolation='none', origin='upper')

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
    
    @staticmethod
    def boundary_array(size: int, width: int) -> array:
        """
        Creates an array of ones with a boundary of zeros

        Parameters
        ----------
        size : int
            size of array (nxn).
        width : int
            width of boundary.

        Returns
        -------
        array
            Array with boundary
        """

        #create base array
        base_array: array = np.zeros((size, size))
        
        #create boundary mask
        boundaries: array = np.ones((size, size), dtype=bool)
        boundaries[base_array.ndim * (slice(width, -width),)] = False

        #set boundary on base array to 0
        base_array[~boundaries] = 1.0

        return base_array

    @staticmethod
    def heatmap(size, updates, sweeps=30, res=0.05, data_file_name='heatmap.dat', heatmap_name='heatmap.png') -> array:
        """
        Creates several Objects with varying update probabilities, updates them over a set amount of time, and produces a heatmap of the average number of infected cells in lattice, with a data file.

        Parameters
        ----------
        size : int
            Size of n x n lattice.
        updates : int
            Number of updates for each simulation/object.
        sweeps : int, optional
            Number of cells to update in an update, by default 30
        res : float, optional
            Size of probability steps, by default 0.05
        data_file_name : str, optional
            Name of output file for data, by default 'heatmap.dat'
        heatmap_name : str, optional
            Name of output file for heatmap plot, by default 'heatmap.png'

        Returns
        -------
        array
            Values used in the heatmap display.
        """
        
        #Initialise lists of probabilities for loops
        p1 = np.arange(0, 1 + res, step=res)
        p2 = np.arange(0, 1 + res, step=res)

        #for output data calculations later
        heatmap = np.zeros((len(p1), len(p2)))

        #open file for data and add header
        file = open(data_file_name, 'w')
        file.write('#p1    p2    Alive%' + '\n')

        #loop through probabilities
        for i in range(len(p1)):

            for j in range(len(p2)):
                
                #start timer to track function whilst running
                start = timer()

                #create GameOfLife object with different probabilities for each loop
                x = GameOfLife(size, sweeps, prob=(p1[i], p2[j]))

                #update for the number of times given to function
                itt: int = 0
                while itt <= updates:

                    x.update(1)
                    itt += 1
                
                #after simulation is finished, add average to array and write to file
                heatmap[i][j] = np.average(x.alive[-updates:]) / size**2
                file.write(str(p1[i]) + '    ' + str(p2[j]) + '    ' + str(heatmap[i][j]) + '\n')

                #finish timer
                end = timer()

                print(f'Simulation no. {i*int(1/res + 1) + j} finished in {round(end-start, 2)}s with {heatmap[i][j]}% of total cells Alive')
        
        file.close()

        #plot heatmap
        plt.imshow(heatmap, origin='lower', extent=[0, 1, 0, 1], cmap='inferno')
        plt.colorbar()
        plt.title('Average Alive Fraction Heatmap')
        plt.xlabel('p2')
        plt.ylabel('p1')
        plt.savefig(heatmap_name)
        plt.show()

        return heatmap
    
    
    @staticmethod 
    def var_heatmap(size, updates, sweeps=30, res=0.05, data_file_name='var_data.dat', heatmap_name='var_heatmap.png'):
        """
        Creates several Objects with varying update probabilities, updates them over a set amount of time, and produces a plot of the variance of infected cells in lattice, with a data file.

        Parameters
        ----------
        size : int
            Size of n x n lattice.
        updates : int
            Number of updates for each simulation/object.
        sweeps : int, optional
            Number of cells to update in an update, by default 30
        res : float, optional
            Size of probability steps, by default 0.05
        data_file_name : str, optional
            Name of output file for data, by default 'var_data.dat'
        heatmap_name : str, optional
            Name of output file for heatmap plot, by default 'var_heatmap.png'

        Returns
        -------
        list
            Values used in the variance plot.
        """
        
        #Initialise list of probabilities for loop
        p1 = np.arange(0, 1 + res, step=res)
        p2 = np.arange(0, 1 + res, step=res)

        heatmap = np.zeros((len(p1), len(p2)))

        #open file for data and add header
        file = open(data_file_name, 'w')
        file.write('#p1    p2    variance\n')

        #loop through probabilities
        for i in range(len(p1)):

            for j in range(len(p2)):
                

                #start timer to track function whilst running
                start = timer()

                #create GameOfLife object with different probabilities for each loop
                x = GameOfLife(size, sweeps, prob=[p1[i], p2[j]])

                #update for the number of times given to function
                itt=0
                while itt <= updates:

                    x.update(1)
                    itt += 1

                #calculate variance
                heatmap[i][j] = np.var(np.array(x.alive) / size**2)
                
                #append values to respective lists and data files
                file.write(f'{p1[i]}    {p2[j]}    {heatmap[i][j]}' + '\n')

                #finish timer
                end = timer()
                print(f'Simulation no. {i*int(1/res + 1) + j} ended successfully in {round(end - start, 2)}s, alive variance: {heatmap[i][j]}')
        

        file.close()

        #plot variance for each probability
        plt.imshow(heatmap, origin='lower', extent=[0, 1, 0, 1], cmap='inferno', aspect='auto')
        plt.plot(np.argmax(heatmap, axis=1)*res, p1, label='Maximal Variance')
        plt.colorbar()
        plt.title('Average Alive Fraction Heatmap')
        plt.xlabel('p2')
        plt.ylabel('p1')
        plt.legend()
        plt.savefig(heatmap_name)
        plt.show()

        #find equation from polyfit
        bestfit = np.polyfit(p1, np.argmax(heatmap, axis=1)*0.05, 2)

        #format and print equation
        equation = str(round(bestfit[0], 2)) + "x**2 + " + str(round(bestfit[1], 2)) + "x + " + str(round(bestfit[2], 2)) 
        print(equation)

        return heatmap

def main():

    _ = GameOfLife.heatmap(100, 5000, 50)
    _ = GameOfLife.var_heatmap(100, 5000, 50)

#to ignore error on importing 
if __name__ == '__main__':

    main()