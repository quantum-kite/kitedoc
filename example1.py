import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb

from pybinding.repository import graphene
from pybinding.repository.graphene import vf
from pybinding.constants import hbar

from export_lattice import Configuration, Calculation, Modification, Disorder, StructuralDisorder, \
    export_lattice, make_pybinding_model

def graphene_initial(onsite=(0, 0)):
    theta = np.pi / 3
    t= 2.8 #eV
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])
    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', - t),
        ([-1, 0], 'A', 'B', - t),
        ([-1, 1], 'A', 'B', - t)
    )

    return lat


lattice = graphene_initial()



nx = 1
ny = 1

lx = 512
ly = 512

# spectrum_range=[e_min, e_max] manually select Hamiltonian bounds.
# if spectrum_range is not selected automatic scaling of a Pybinding model with equivalent disorder strength
# is done in the background.
configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                              is_complex=False, precision=1)

calculation = Calculation(configuration)
calculation.dos(num_points=1000, num_moments=512, num_random=1, num_disorder=1)

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = Modification(magnetic_field=False)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
export_lattice(lattice, configuration, calculation, modification, 'example1.h5')
