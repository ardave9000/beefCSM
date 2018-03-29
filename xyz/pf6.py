from ase import Atoms
import ase.visualize
import ase.io

d=1.59
h=0.2
vac=1
pf6 = Atoms('PF6',positions=[(0, 0, 0), (0, 0, d), (0, 0, -d), (d, 0, 0), (-d, 0, 0), (0, d, 0), (0, -d, 0)])
ase.io.write('pf6.xyz',pf6)

#pf6.center(vacuum=vac)
#ase.visualize.view(pf6)
