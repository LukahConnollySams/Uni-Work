from array import array
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


class Lattice(object):

    def __init__(self, n: int) -> None:

        self.n = n

        #creates nxn with random 0's or 1's
        self.lattice = np.random.choice([0, 1], size = (self.n, self.n))

    def __str__(self) -> str:
        """
        String of object

        Returns
        -------
        str
            String of object
        """
        return f'{self.lattice}'

    def neighbours(self) -> list:
        """
        Creates an array of same shape as lattice, but with value equal to a list of neighbouring values

        Returns
        -------
        list
            neighbours of each cell in lattice
        """
        neighbours: list = []

        #determine values, and append, 8 neighbouring cells
        neighbours.append(np.roll(self.lattice, 1, axis=1))
        neighbours.append(np.roll(self.lattice, -1, axis=1))
        neighbours.append(np.roll(self.lattice, 1, axis=0 ))
        neighbours.append(np.roll(self.lattice, -1, axis=0))


        return neighbours

    def update(self, x) -> array:
        """
        Updates the entire lattice

        Parameters
        ----------
        x : int
            itterable used for animation

        Returns
        -------
        list
            Updated lattice
        """
        

        return self.lattice

    def animate(self, i: int) -> list:
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

    def display(self, interval: int, frames: int) -> None:
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
        self.anim = FuncAnimation(fig, self.animate, frames=frames, repeat=False, interval=interval, blit=False)

        #show plt figure and animation
        plt.show()
        plt.close()