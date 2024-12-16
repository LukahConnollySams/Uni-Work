
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from timeit import default_timer as timer
import sys
from scipy.signal import find_peaks


class GameOfLife(object):

    def __init__(self, n, lattice_type=None) -> None:

        self.n = n
        self.lattice_type = lattice_type
        
        if lattice_type == 'glider':
            
            #create array with simple glider in top left
            self.lattice = np.zeros((self.n, self.n))
            pos = [(0, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

            #takes positions above and applies them to lattice as alive cells
            for x in pos:

                self.lattice[x] = 1

        elif lattice_type == 'oscillator':
            
            #create array with sample oscillator in top left
            self.lattice = np.zeros((self.n, self.n))
            pos = [(1, 0), (1, 1), (1, 2)]

            #takes positions above and applies them to lattice as alive cells
            for x in pos:

                self.lattice[x] = 1

        else:
            
            #creates nxn with random 0's or 1's, with lattice_type as probability for cells = 1
            self.lattice = np.random.choice([0, 1], size = (self.n, self.n))

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

    def sum_neighbours(self,):
        """
        Creates an array of same shape as lattice, but with value equal to number of alive neighbours

        Returns
        -------
        list
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

    def update(self, x):
        """
        Updates the entire lattice of the game of life

        Parameters
        ----------
        x : int
            itterable used for animation

        Returns
        -------
        list
            Updated lattice
        """

        # used for later calculations
        updated_alive = 0
        #calculate how many neighbouring cells are alive
        sum = self.sum_neighbours()

        for i in range(self.n):

            for j in range(self.n):

                

                #if cell is dead
                if self.lattice[i][j] == 0:

                    #make cell alive if 3 neighbours are also
                    if sum[i][j] == 3:

                        self.lattice[i][j] = 1 

                        #add to current alive cell count
                        updated_alive +=1

                #if cell is alive
                else:

                    #kill if surrounded by more than 3 alive neighbours
                    if sum[i][j] > 3:

                        self.lattice[i][j] = 0

                    #kill if surrounded by less than 2 alive neighbours
                    elif sum[i][j] < 2:

                        self.lattice[i][j] = 0

                    #stay alive if surrounded by 2 or 3 neighbours
                    else:

                        self.lattice[i][j] = 1 

                        #add to current alive cell count
                        updated_alive +=1


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
    
    def glider_com(self, steps, no_peak=2):
        """
        Calculates the Centre of Mass Velocity of a glider in a lattice.
        Parameters
        ----------
        steps : int
            number of updates for simulation
        no_peak : int, optional
            number of peaks and troughs to use for velocity calculations, by default 2

        Returns
        -------
        float
            Average centre of mass velocity of glider
        list
            Distances of the centre of mass of the glider
        list
            Update numbers/tiems of the glider positions

        Raises
        ------
        ValueError
            If object is not of type 'glider'
        """

        #raise error if not a glider type lattice (for fun)
        if self.lattice_type != 'glider':

            raise ValueError(f"Object has wrong lattice type for method, recieved {self.lattice_type}, must be 'glider'")

        #lists used later
        r_cm = []
        time = []

        #loop for the number of simulation steps
        x = 0
        while x <= steps:

            gliders_pos = []

            for i in range(self.n):

                for j in range(self.n):

                    #find positions of all alive cells and append distance to list
                    if self.lattice[i][j] == 1.0:

                        gliders_pos.append(np.linalg.norm([i, j]))
            
            
            #calculate centre of mass 
            r_cm.append(np.sum(gliders_pos) / len(gliders_pos))
            time.append(x)

            self.update(1)
        
            x += 1

        #find peaks and troughs (for boundary conditions) and fomrat appropriately
        peaks = list(find_peaks(r_cm)[0])
        troughs = list(find_peaks(-1*np.array(r_cm))[0])
        troughs.insert(0, 0)
        peaks = peaks[:no_peak]
        troughs = troughs[:no_peak]

        vels = []

        #calculate velocities within boundaries
        for i in range(len(troughs)):

            dist = r_cm[peaks[i]] - r_cm[troughs[i]]
            vels.append(dist / (peaks[i] - troughs[i]))

        #average the velocities
        com_vel = np.average(vels)

        return com_vel, r_cm, time
    
    @staticmethod
    def equil(sim_no, n, patience): 
        """
        Calculates the time it takes for the simulation to reach equilibrium.

        Parameters
        ----------
        sim_no : int
            number of updates to perform
        n : int
            size of lattice
        patience : int
            number of updates wait before determining simulation is in equilibrium

        Returns
        -------
        list
            list of equilibrium times
        """
        equilibriums = []
        
        #loop through number of simulations
        for i in range(sim_no):

            #start timer for proceeding code
            start = timer()
            
            #create an instance of the class
            x = GameOfLife(n)

            #update lattice initially for a number of times equal to patience 
            for j in range(patience):

                x.update(1)

            #checks amount of unique values in the final section of the self.alive variable with a size of patience 
            while not len(set(x.alive[-patience:])) <= 2:
               
                x.update(1)

            #append the actual length of simulation to reach equilibrium
            equilibriums.append((len(x.alive) - patience + 1))

            #end timer and print to check for errors
            end = timer()
            print(f'simulation number: {i+1}, ended successfully in {end - start}s')
            

        return equilibriums


def main():

    #read commandline arguements
    inputs = sys.argv

    #redirects user to call help if unsure how to run file
    if len(inputs) <= 1:

        print("Use 'help' for info on how to run file")

    #help calls for instructions to be printed to user for understanding of accepted formats
    elif inputs[1] == 'help':

        print("To run file call one of below with their respective format:")  
        print("For example animation: file animation lattice_size(int) lattice_type(str, or float) time_between_frames(int) number_of_updates(int)")
        print("For equilibration plot: file equilibrium no_of_simulations(int) lattice_size(int) patience(int)")
        print("For COM Velocity plot: file com lattice_size(int) glider no_of_updates(int)")

    #creates example animation using .display() with system arguements
    elif inputs[1] == 'animation':

        x = GameOfLife(int(inputs[2]), inputs[3])
        x.display(int(inputs[4]), int(inputs[5]))

    #creates a plot of the equilibration times over multiple GameOfLife simulations using system arguements
    elif inputs[1] == 'equilibrium':

        equils = GameOfLife.equil(int(inputs[2]), int(inputs[3]), int(inputs[4]))

        plt.hist(equils, bins=20, density=True)
        plt.title('Equilibrium Point Densities')
        plt.xlabel('Equilibrium')
        plt.ylabel('Density')
        plt.show()

    #creates a centre of mass plot for an example glider using system arguements
    elif inputs[1] == 'com':
        
        x = GameOfLife(int(inputs[2]), inputs[3])
        vel, r_cm, time = x.glider_com(int(inputs[4]))
        print(f'Glider speed: {round(vel, 2)}')
        plt.plot(time, r_cm)
        plt.title('Glider COM over time')
        plt.xlabel('Time')
        plt.ylabel('r_cm')
        plt.show()


#to ignore error on importing 
if __name__ == '__main__':

    main()