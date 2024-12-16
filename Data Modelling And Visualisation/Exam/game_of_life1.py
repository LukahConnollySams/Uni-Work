from array import array
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


class GameOfLife(object):

    def __init__(self, size: int, n: int, prob: float=0.01, lattice=None) -> None:

        self.size = size
        self.n = n

        #feature to input custom initial lattice
        if lattice is not None:

            self.lattice = lattice

        else:

            #creates size x size grid with random 0's or 1's, with prob as probability for cells = 1
            self.lattice = np.random.choice([0, 1], size=(self.size, self.size), p=(1 - prob, prob))

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

        #for every cell in lattice
        for i in range(self.size):

            for j in range(self.size):

                #if cell is dead
                if self.lattice[i][j] == 0:

                    #make cell alive if n neighbours are also
                    if sum[i][j] == self.n:

                        self.lattice[i][j] = 1 

                        #add to current alive cell count
                        updated_alive +=1

                #if cell is alive, kill
                else:

                    self.lattice[i][j] = 0


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

    def display(self, interval, frames) -> None:
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


def main():
    
    #n=2 simulations
    x = GameOfLife(100, 2)
    x.display(1, 5000)

    #plot data
    plt.plot(np.arange(0, len(x.alive)), x.alive)
    plt.title('No. of Alive Cells for n=2')
    plt.xlabel('Number of Updates')
    plt.ylabel('No. of alive Cells')
    plt.savefig('alive_cells_n=2.png')
    plt.show()

    #n=3 simulations
    x = GameOfLife(100, 3)
    x.display(1, 5000)

    #plot data
    plt.plot(np.arange(0, len(x.alive)), x.alive)
    plt.title('No. of Alive Cells for n=3')
    plt.xlabel('Number of Updates')
    plt.ylabel('No. of alive Cells')
    plt.savefig('alive_cells_n=3.png')
    plt.show()

    #box simulation
    box_lattice = GameOfLife.boundary_array(50, 15)

    x = GameOfLife(len(box_lattice), 2, lattice=box_lattice)
    x.display(1, 5000)

    #plot data
    plt.plot(np.arange(0, len(x.alive)), x.alive)
    plt.title('No. of Alive Cells for n=2, Box')
    plt.xlabel('Number of Updates')
    plt.ylabel('No. of alive Cells')
    plt.savefig('alive_cells_n=2_box.png')
    plt.show()


#to ignore error on importing 
if __name__ == '__main__':

    main()