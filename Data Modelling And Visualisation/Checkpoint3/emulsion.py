from array import array
from lattice import Lattice
import numpy as np
from matplotlib import pyplot as plt


class Emulsion(Lattice):

    a: float = 0.1
    k: float = 0.1
    M: float = 0.1

    def __init__(self, n: int, time_step: float, len_step: float, sweeps: int = 1, thi_init: tuple = (0.5, 0.001)) -> None:

        self.n = n
        self.t_0: int = 0
        self.time_step = time_step
        self.len_step = len_step
        self.sweeps = sweeps
        self.lattice: array = (np.random.rand(self.n, self.n) - 0.5) * thi_init[1] + thi_init[0]
        self.free_e: list = []
        self.free_e.append(self.free_energy_density())
        self.updates: list = []
        self.updates.append(0)

    def __str__(self) -> str:
        return super().__str__()
    
    def neighbours(self, arr: list=None) -> array:

        if arr is not None:

            neighbours: list = []

            #determine values, and append, 4 neighbouring cells
            neighbours.append(np.roll(arr, 1, axis=1))
            neighbours.append(np.roll(arr, -1, axis=1))
            neighbours.append(np.roll(arr, 1, axis=0 ))
            neighbours.append(np.roll(arr, -1, axis=0))

            return np.sum(neighbours, axis=0)

        else:

            return super().neighbours()

    def update(self, x) -> array:

        for i in range(self.sweeps):

            self.thi = self.comp_order()
        
        self.free_e.append(self.free_energy_density())
        self.updates.append(self.sweeps * (len(self.updates) + 1))

        return self.thi

    def animate(self, i: int) -> list:
        return super().animate(i)

    def display(self, interval: int, frames: int) -> None:
        return super().display(interval, frames)
    
    def comp_order(self) -> array:

        new_mu = self.chem_potential()

        thi_new: array = self.lattice + (Emulsion.M * self.time_step / (self.len_step ** 2)) * (self.neighbours(arr=new_mu)- 4 * new_mu)

        return thi_new

    def chem_potential(self) -> array:

        mu_new: array = -Emulsion.a * self.lattice + Emulsion.a * self.lattice ** 3 - (Emulsion.k / (self.len_step ** 2)) * (self.neighbours(arr=self.lattice) - 4 * self.lattice)

        return mu_new

    def free_energy_density(self) -> float:

        free_e: array = (Emulsion.a / 2) * (0.5 * self.lattice ** 4 + -1 * self.lattice ** 2) + (Emulsion.k / 2) * (self.neighbours(arr=self.lattice) - 4 * self.lattice) ** 2

        return np.sum(free_e.flatten())


def main():


    x = Emulsion(100, 1.5, 1, sweeps=250, thi_init=(0.5, 0.001))
    x.display(1, 200)
    plt.close()

    plt.plot(x.updates, x.free_e)
    plt.title('Free Energy over Time')
    plt.xlabel('Number of updates')
    plt.ylabel('Free Energy')
    plt.show()

if __name__ == '__main__':

    main()