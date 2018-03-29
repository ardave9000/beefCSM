from gpaw import GPAW, Mixer
from gpaw.utilities import h2gpts
from ase import Atoms
from ase.units import mol,kJ,kcal,Pascal,m
from ase.data.vdw import vdw_radii
from ase.parallel import parprint,paropen
from ase.dft.bee import BEEFEnsemble
import numpy as np
import ase.io
import ase.visualize
import pickle
import datetime
import sys
import math
from gpaw.xc.vdw import VDWFunctional
from gpaw.solvation import (
    SolvationGPAW, #solvation calculator
    EffectivePotentialCavity, #cavity form with an effective potential
    Power12Potential, #a specific effective potential (Lennard Jones)
    LinearDielectric, #rule to construct permittivity function from cavity
    GradientSurface, #rule to calculate cavity surface area
    SurfaceInteraction) #rule to calculate non-electrostatic interaction)

if len(sys.argv)!=5:
    raise ValueError("must supply sys.argv in order: solvent name | monatomic ion symbol with +/- | h-step (grid spacing) | unit cell dimension / 2")
else:
    sv = sys.argv[1]
    s = sys.argv[2]
    h = sys.argv[3]
    v = sys.argv[4]
    sv.strip()
    s.strip()
    charge_str = s[-1]
    sym = s[:-1]
    if charge_str=='+':
        charge = 1
    else:
        charge = -1

now = datetime.datetime
timestamp = now.now().strftime('%Y%-m%-d%-I%-M%-S')
shell = ase.io.read("~/xyz/{}4x.xyz".format(sv))
solvent = ase.io.read("~/xyz/{}.xyz".format(sv))
solventparams = pickle.load(open("~/beefCSM/parameters.pkl","r"))

parprint("completed")
'''
#Solvent parameters for DMSO (J. Chem. Phys. 141, 174108 (2014))
sv = 'dmso'
nd = 14.08
R = 8.6173303e-5
u0 = 0.20 #eV
epsinf = 47.2 #dielectric constant
gamma = 8.8*1e-3*Pascal*m #conversion from dyne/cm to eV/Angstrom**2
T = 298.15 #Kelvin
vdw_radii = vdw_radii.copy()
vdw_radii[1] = 1.09
atomic_radii = lambda atoms: [vdw_radii[n] for n in atoms.numbers] #create callable for GPAW calculator

#All parameters in ASE unit conventions (eV, Angstrom...)
#DFT calculation parameters
h = 0.15
vac = 6
xc = 'PBE'
beef = 'BEEF-vdW'

#Each calculation will be converged with PBE, then, using same calculator, be run with BEEF.
#This enables the PBE WFs to be used as initial guesses for BEEF convergence, hopefully enabling large box sizes.

#Data is stored in a dictionary with key: input sys.argv string "s" indicating ion-species. value : tuple of one value and one array,
#indicating Gsol and dGsol.
output = {}

#Perform solvent-only (gas phase) calculation
solvent.center(vacuum=vac) #set box size
calc = GPAW(xc=xc, txt='{}_{}_gas_{}.txt'.format(s,sv,timestamp),gpts=h2gpts(h, solvent.get_cell(), idiv=16))
solvent.calc = calc
solvent.get_potential_energy() #converge WFs with PBE

calc.set(xc=beef)
solvent.calc = calc
E_sv = solvent.get_potential_energy() #BEEF best fit energy
ense = BEEFEnsemble(calc)
dE_sv = ense.get_ensemble_energies() #N=2000 non-self-consist calculation to generate uncertainty spread

#Perform solvent-only (solvated phase) calculation
scalc = SolvationGPAW(
        xc=xc, txt='{}_{}_solvated_{}.txt'.format(s,sv,timestamp),
        gpts=h2gpts(h, solvent.get_cell(), idiv=16),
        cavity=EffectivePotentialCavity(
            effective_potential=Power12Potential(atomic_radii,u0),
            temperature=T,
            surface_calculator=GradientSurface()),
        dielectric=LinearDielectric(epsinf=epsinf),
        interactions=[SurfaceInteraction(surface_tension=gamma)])
solvent.calc = scalc
solvent.get_potential_energy()

scalc.set(xc = beef)
solvent.calc = scalc
E_sv_sol = solvent.get_potential_energy()
ense = BEEFEnsemble(scalc)
dE_sv_sol = ense.get_ensemble_energies() 

#Perform single-ion (gas phase) calculation
ion = Atoms(sym)
calc.set(xc=xc,charge = charge)
ion.calc = calc
ion.center(vacuum=vac)
ion.get_potential_energy()

calc.set(xc=beef)
ion.calc = calc
E_ion = ion.get_potential_energy() 
ense = BEEFEnsemble(calc)
dE_ion = ense.get_ensemble_energies()

#Perform ion-in-solvent-cluster (gas phase) calculation
ion = Atoms(sym)
ion.extend(shell)
calc.set(xc=xc)
ion.calc = calc
ion.center(vacuum=vac)
ion.get_potential_energy()

calc.set(xc=beef)
ion.calc = calc
E_cluster = ion.get_potential_energy()
ense = BEEFEnsemble(calc)
dE_cluster = ense.get_ensemble_energies()

#Perform ion-in-solvent-cluster (solvated phase) calculation
scalc.set(xc=xc,charge=charge)
ion.calc = scalc
ion.get_potential_energy()

scalc.set(xc=beef)
ion.calc=scalc
E_cluster_sol = ion.get_potential_energy()
ense = BEEFEnsemble(scalc)
dE_cluster_sol = ense.get_ensemble_energies()

#Calculate solvation energy according to thermodynamic cycle
Gsolsv = E_sv_sol - E_sv
Gclust = E_cluster - E_ion - 4*E_sv
Gsolcluster = E_cluster_sol - E_cluster
Gsol = Gsolcluster+Gclust - 4*Gsolsv-4*R*T*math.log(nd)


dGsolsv = dE_sv_sol[:] - dE_sv[:]
dGclust = dE_cluster[:] - dE_ion[:] - 4*dE_sv[:]
dGsolcluster = dE_cluster_sol[:] - dE_cluster[:]
dGsol = dGsolcluster[:]+dGclust[:] - 4*dGsolsv[:]

parprint ("Gsol({}) = {:.6f} + {:.6f}".format(s,Gsol,dGsol.var()))

output[s] = (Gsol,dGsol)

#Open pickle file for output dictionary and dump
f = paropen("BEEF_{}_{}.pkl".format(s,sv),mode='w')
pickle.dump(output,f)
'''
