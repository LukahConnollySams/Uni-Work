from array import array
from lattice import Lattice
import numpy as np
from matplotlib import pyplot as plt


class Poisson(object):

    def __init__(self, n: int, mode: str, time_step: float, len_step: float, tolerance: float=0.01) -> None:

        self.n: int = n
        self.time_step = time_step
        self.len_step = len_step
        self.e_pot: array = self.boundary_array((n, n, n))
        self.rho: array = np.zeros((n, n, n))
        self.rho[n//2, n//2, n//2] = 1
        self.e_field: array = self.e_field_update()
        self.updates: list = []
        self.updates.append(0)
        self.mode: str = mode
        self.tolerance: float = tolerance
        self.omega: float = 1.88

    def __str__(self) -> str:
        return f'Electrostatic Potential:\n{self.e_pot}\n Electric Field:\n{self.e_field}'

    def neighbours(self, arr: list = None) -> array:
        
        neighbours =[]

        arr = arr[arr.ndim * (slice(1, -1),)]

        neighbours.append(np.roll(arr, 1, axis=1))
        neighbours.append(np.roll(arr, -1, axis=1))
        neighbours.append(np.roll(arr, 1, axis=0 ))
        neighbours.append(np.roll(arr, -1, axis=0))
        neighbours.append(np.roll(arr, 1, axis=2))
        neighbours.append(np.roll(arr, -1, axis=2))

        neighbours = np.sum(neighbours, axis=0)

        return np.pad(neighbours, 1, mode='constant')

    def electrostatic_potential_update(self) -> array:

        return (self.neighbours(arr=self.e_pot) + self.rho) / 6

    def e_field_update(self) -> array:

        self.e_field = -1 * np.gradient(self.e_pot)

        return self.e_field

    def update(self, steps) -> None:

        x = 0

        while True:

            old_lattice = self.e_pot

            if self.mode == 'jacobi':
                self.e_pot = self.electrostatic_potential_update()
                self.updates.append(len(self.updates))

            elif self.mode == 'gs':

                for i in range(1, self.n-1):
                    for j in range(1, self.n-1):
                        for k in range(1, self.n-1):

                            self.e_pot[i, j, k] = ((self.e_pot[i + 1, j, k] + self.e_pot[i - 1, j, k] + 
                                                    self.e_pot[i, j + 1, k] + self.e_pot[i, j - 1, k] + 
                                                    self.e_pot[i, j, k + 1] + self.e_pot[i, j, k - 1] + self.rho[i, j, k]) / 6) * self.omega + (1 - self.omega) * self.e_pot[i, j, k]
                
                self.updates.append(len(self.updates))
            
            x += 1

            if np.sum(np.abs(self.e_pot - old_lattice)) <= self.tolerance or x == steps:

                print(np.sum(np.abs(self.e_pot - old_lattice)))
                self.e_field_update()

                break

    @staticmethod
    def boundary_array(size: int, shift: float=0) -> array:

        base_array: array = np.zeros(size)
        boundaries: array = np.ones(size, dtype=bool)

        boundaries[base_array.ndim * (slice(1, -1),)] = False

        base_array[~boundaries] = np.random.uniform(-0.1, 0.1, (~boundaries).sum()) + shift

        return base_array

def main():

    x = Poisson(5, 'jacobi', 0.1, 1)
    x.update(200000)
    print(x.updates[-1])


if __name__ == '__main__':

    main()