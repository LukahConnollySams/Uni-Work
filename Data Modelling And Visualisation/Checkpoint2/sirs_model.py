
from game_of_life import GameOfLife
import numpy as np
import random as R
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
from timeit import default_timer as timer

class SIRS(GameOfLife):


    def __init__(self, n, probabilities, sweeps, immune=None) -> None:

        self.n = n
        self.s_rate = probabilities[0]
        self.i_rate = probabilities[1]
        self.r_rate = probabilities[2]
        self.sweeps = sweeps

        if immune is not None:
            
            self.immune = immune

        else:

            self.immune = 0
            
        self.lattice = np.random.choice([0, 1, 2], size = (self.n, self.n))

        self.sus = []
        self.sus.append(self.lattice.flatten().tolist().count(0.0))
        self.inf = []
        self.inf.append(self.lattice.flatten().tolist().count(1.0))
        self.rec = []
        self.rec.append(self.lattice.flatten().tolist().count(2.0))

    def __str__(self) -> str:
        return super().__str__()
    
    def random_p(self):
        """
        Generates random indices for 1 point on a 2D array.

        Returns
        -------
        array
            Array containing indices for 2D array.
        """

        #random 2 values from the idices of the list
        return  np.random.choice(self.n, 2)

    def sum_neighbours(self):
        """
        Generates a list for each cell in the object's lattice of the state of their 4 neighbours combined into one array of shape (n, n).

        Returns
        -------
        List
            Lists of nieghbour's states of each cell in lattice.
        """

        #empty list for storing later values
        neighbours = np.zeros((self.n, self.n)).tolist()

        #get state of each neighbour
        n1 = np.roll(self.lattice, 1, axis=1)
        n2 = np.roll(self.lattice, -1, axis=1)
        n3 = np.roll(self.lattice, 1, axis=0 )
        n4 = np.roll(self.lattice, -1, axis=0)

        #loop for each cell
        for i in range(self.n):

            for j in range(self.n):

                #change to a list of neighbours for each cell
                neighbours[i][j] = [n1[i][j], n2[i][j], n3[i][j], n4[i][j]]

        return neighbours

    def update(self, itter):
        """
        Updates a set amount of cells in the SIRS lattice. 

        Parameters
        ----------
        itter : int
            Used for animation itteration, not needed to run main component of function.

        Returns
        -------
        Array
            Updated lattice of the object
        """

        #get stat of neighbours for each cell
        neighbours = self.sum_neighbours()
        
        #update set amount of cells set by sweeps
        x = 0
        while x <= self.sweeps:
            
            #get random point
            xy = self.random_p()

            #check value in cell
            home = self.lattice[xy[0]][xy[1]]

            #if cell is susceptible
            if home == 0:
                
                #check if there are infected neighbours
                if neighbours[xy[0]][xy[1]].count(1.0) > 0:
                    
                    if self.s_rate > R.random():

                        #change to infected
                        self.lattice[xy[0]][xy[1]] = 1

            #if cell infected
            elif home == 1:

                if self.i_rate > R.random():
                    
                    #change to recovered                    
                    self.lattice[xy[0]][xy[1]] = 2

                

            #if cell is recovered
            elif home == 2:

                if self.r_rate > R.random():

                    #change to immune
                    if self.immune != 0 and self.lattice.flatten().tolist().count(3.0)/self.n**2 < self.immune:
                        
                        self.lattice[xy[0]][xy[1]] = 3

                    #change to susceptible
                    else:

                        self.lattice[xy[0]][xy[1]] = 0

            else:

                continue

            

            x += 1

        #add number of each state to respective lists
        self.sus.append(self.lattice.flatten().tolist().count(0.0))
        self.inf.append(self.lattice.flatten().tolist().count(1.0))
        self.rec.append(self.lattice.flatten().tolist().count(2.0))
            

        return self.lattice

    def animate(self, i):
        return super().animate(i)

    def display(self, interval, frames):
        return super().display(interval, frames)

    @staticmethod
    def avg_inf_heatmap(n, updates, patience=100, sweeps=30, res=0.05, data_file_name='SIRS_avg_heatmap.dat', heatmap_name='SIRS_avg_heatmap.png'):
        """
        Creates several Objects with varying update probabilities, updates them over a set amount of time, and produces a heatmap of the average number of infected cells in lattice, with a data file.

        Parameters
        ----------
        n : int
            Size of n x n lattice.
        updates : int
            Number of updates for each simulation/object.
        patience : int, optional
            How many updates to preform before equilibrium, by default 100
        sweeps : int, optional
            Number of cells to update in an update, by default 30
        res : float, optional
            Size of probability steps, by default 0.05
        data_file_name : str, optional
            Name of output file for data, by default 'SIRS_avg_heatmap.dat'
        heatmap_name : str, optional
            Name of output file for heatmap plot, by default 'SIRS_avg_heatmap.png'

        Returns
        -------
        array
            Values used in the heatmap display.
        """
        
        #Initialise lists of probabilities for loops
        p1 = np.arange(0, 1 + res, step=res)
        p3 = np.arange(0, 1 + res, step=res)

        #for output data calculations later
        infected_heatmap = np.zeros((len(p1), len(p3)))

        #open file for data and add header
        file = open(data_file_name, 'w')
        file.write('#p1    p3    infected%')

        #loop through probabilities
        for i in range(len(p1)):

            for j in range(len(p3)):
                
                #start timer to track function whilst running
                start = timer()
                #create SIRS object with different probabilities for each loop
                x = SIRS(n, [p1[i], 0.5, p3[j]], sweeps)

                #update for a period of patience time
                for _ in range(patience):

                    x.update(1)

                #update for the number of times given to function
                itt=0
                while itt <= updates:

                    x.update(1)
                    itt += 1
                
                #after simulation is finished, add average to array and write to file
                infected_heatmap[i][j] = np.average(x.inf[-updates:]) / n**2
                file.write(str(p1[i]) + '    ' + str(p3[j]) + '    ' + str(infected_heatmap[i][j]) + '\n')
                #finish timer
                end = timer()

                print(f'Simulation no. {i*int(1/res + 1) + j} finished in {round(end-start, 2)}s with {infected_heatmap[i][j]}% of total cells infected')
        
        file.close()

        #plot heatmap
        plt.imshow(infected_heatmap, origin='lower', vmin=0, vmax=1.0, extent=[0, 1, 0, 1], cmap='inferno')
        plt.colorbar()
        plt.title('Infection and Susceptibility Probability Heatmap')
        plt.xlabel('p3')
        plt.ylabel('p1')
        plt.savefig(heatmap_name)
        plt.show()

        return infected_heatmap

    @staticmethod
    def var_inf_heatmap(n, updates, repeats=5, patience=100, sweeps=30, res=0.05, data_file_name='SIRS_var_data.dat', heatmap_file_name='SIRS_var_heatmap.png'):
        """
        Creates several Objects with varying update probabilities, updates them over a set amount of time, and produces a plot of the variance of infected cells in lattice, with a data file.

        Parameters
        ----------
        n : int
            Size of n x n lattice.
        updates : int
            Number of updates for each simulation/object.
        patience : int, optional
            How many updates to preform before equilibrium, by default 100
        sweeps : int, optional
            Number of cells to update in an update, by default 30
        res : float, optional
            Size of probability steps, by default 0.05
        data_file_name : str, optional
            Name of output file for data, by default 'SIRS_var_heatmap.dat'
        heatmap_name : str, optional
            Name of output file for heatmap plot, by default 'SIRS_var_heatmap.png'

        Returns
        -------
        list
            Values used in the variance plot.
        """
        
        #Initialise list of probabilities for loop
        p1 = np.arange(0.2, 0.5 + res, res)

        #open file for data and add header
        file = open(data_file_name, 'w')
        file.write('#p1    variance_avg    variance_stdv\n')

        #for output data calculations later
        var_avg = []
        var_stdv = []

        #loop through probabilities
        for i in range(len(p1)):
            
            #temporary list for data storage
            temp_var = []

            #start timer to track function whilst running
            start = timer()
            file.write(str(p1[i]))

            #loop through number of repeats
            for j in range(repeats):

                #create SIRS object with different probabilities for each loop
                x = SIRS(n, [p1[i], 0.5, 0.5], sweeps)

                #update for a period of patience time
                for __ in range(patience):

                        x.update(1)

                #update for the number of times given to function
                itt=0
                while itt <= updates:

                    x.update(1)
                    itt += 1

                #after simulation is finished, add variance to list and write to file
                temp_var.append((np.average(np.square(x.inf[-updates:])) - np.average(x.inf[-updates:])**2) / n**2)
            
            #append values to respective lists and data files
            var_avg.append(np.average(temp_var))
            var_stdv.append(np.std(temp_var))
            file.write('    ' + str(var_avg[i]) + '    ' + str(var_stdv[i]) + '\n')

            #finish timer
            end = timer()
            print(f'{repeats} repeat(s) of simulation {i+1} ended successfully in {round(end - start, 2)}s, infected variance: {var_avg[i]}')
        

        file.close()

        #plot variance for each probability
        plt.errorbar(p1, var_avg, yerr=var_stdv)
        plt.title('Infected variance as a function of p1')
        plt.xlabel('p1')
        plt.ylabel('Infected variance')
        plt.savefig(heatmap_file_name)
        plt.show()

        return p1, var_avg

    @staticmethod
    def immune_frac(n, updates, patience=100, sweeps=30, res=0.05, data_file_name='SIRS_frac_inf_data.dat', heatmap_file_name='SIRS_frac_inf_heatmap.png'):
        """
        Creates a plot and data file of the fraction of immunity and the average number of infected cells.

        Parameters
        ----------
        n : int
            Size of n x n lattice.
        updates : int
            Number of updates for each simulation/object.
        patience : int, optional
            How many updates to preform before equilibrium, by default 100
        sweeps : int, optional
            Number of cells to update in an update, by default 30
        res : float, optional
            Size of probability steps, by default 0.05
        data_file_name : str, optional
            Name of output file for data, by default 'SIRS_avg_heatmap.dat'
        heatmap_name : str, optional
            Name of output file for heatmap plot, by default 'SIRS_avg_heatmap.png'

        Returns
        -------
        List
            Lists of data used in plot
        """

        #list for data storage
        imm_frac = np.arange(0, 1 + res, res)
        inf_frac = []

        #create file with header
        file = open(data_file_name, 'w')
        file.write('#imm_frac    infrac\n')

        #loop through immunity fraction values
        for i in range(len(imm_frac)):
            
            #start time and initiailse object
            start = timer()
            x = SIRS(n, [0.5, 0.5, 0.5], sweeps, immune=imm_frac[i])

            #loop for an initial amountof time
            for _ in range(patience):

                x.update(1)

            #loop fot eh number of updates given
            itt = 0
            while itt <= updates:

                x.update(1)
                itt += 1
            
            #append data to lists and files
            inf_frac.append(np.average(x.inf[-updates:]) / n**2)
            file.write(f'{inf_frac[i]}    {np.average(x.inf[-updates:]) / n**2}\n')
            
            #end timer
            end = timer()
            #print statement to track progress over long wait times
            print(f'Simulation number {i+1}, ended successfully in {end-start}s')

        file.close()

        #plot data
        plt.plot(imm_frac, inf_frac)
        plt.title('Infected Fraction vs Fraction of Immunity')
        plt.xlabel('Fraction of Immunity')
        plt.ylabel('Infected Fraction')
        plt.savefig(heatmap_file_name)
        plt.show()

        return imm_frac, inf_frac

def main():

    #read commandline arguements
    inputs = sys.argv

    #redirects user to call help if unsure how to run file
    if len(inputs) <= 1:

        print("Use 'help' for info on how to run file")

    #help calls for instructions to be printed to user for understanding of accepted formats
    elif inputs[1] == 'help':

        print("To run file call one of below with their respective format:") 
        print('For an example animation: file animation lattice_size(int) p1(float) p2(float) p3(float) sweeps(int) time_between_frames(int) number_of_updates(int)')
        print('For Average Heatmap: file avg_heatmap lattice_size(int) no_of_updates(int)')
        print('For Variance Plot: file var_plot lattice_size(int) no_of_updates(int)')
        print('For fraction of immunity plot: file immune_frac lattice_size(int) no_of_updates(int)')

    #creates example animation using .display() with system arguements
    elif inputs[1] == 'animation':

        x = SIRS(inputs[2], [inputs[3], inputs[4], inputs[5]], inputs[6])
        x.display(inputs[7], inputs[8])

    #creates a plot of the average number of infected cells over multiple SIRS simulations using system arguements
    elif inputs[1] == 'avg_heatmap':

        _ = SIRS.avg_inf_heatmap(inputs[2], inputs[3])

    #creates a plot of the average variance of the infected cells over mutliple SIRS simulations using system arguements
    elif inputs[1] == 'var_plot':

        _, __ = SIRS.var_inf_heatmap(inputs[2], inputs[3], res=inputs[4])

    #creates a plot of the immune fraction and the number of infected cells over mutliple SIRS simulations using system arguements
    elif inputs[1] == 'immune_frac':

        _, __ = SIRS.immune_frac(inputs[2], inputs[3])

if __name__ == '__main__':

    main()