from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

# Load the PDB file into an object
pdb = app.PDBFile('trpcage.pdb')

# Load the force field file into an object
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

# Create system object using information in the force field:
# forcefield: contains parameters of interactions
# topology: lists of atoms, residues, and bonds
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)

# Create a Langevin integrator for temperature control
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)

# Add a Monte Carlo barostat to the system for pressure control
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))

custom_bond = mm.CustomBondForce("0.5*k*(r-r0)^2")
custom_bond.addGlobalParameter("k", 100)
custom_bond.addGlobalParameter("r0", 1)
custom_bond.addBond(2, 275)

system.addForce(custom_bond)

# Use the CPU platform
platform = mm.Platform.getPlatformByName('CPU')

### If you want to add any forces to your System or modify any
### of the existing forces, you should do it here - after the
### System has been created, but before the Simulation is created.

# Implemented option B, to create a merge conflict.
def neutralizeCharge():
    atoms = list(pdb.topology.atoms())
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            print("Found nonbonded force")
            for i in range(f.getNumParticles()):
                if atoms[i].residue.name != 'HOH':
                    chg, sig, eps = f.getParticleParameters(i)
                    f.setParticleParameters(i, 0.0, sig, eps)    

# Create a Simulation object by putting together the objects above
simulation = app.Simulation(pdb.topology, system, integrator, platform)

# Set positions in the Simulation object
simulation.context.setPositions(pdb.positions)

# Minimize the energy of the system (intentionally doing a rough minimization)
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=20, tolerance=100)

# Initialize the random velocities of the system from a Maxwell-Boltzmann distribution
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# Add reporters to the simulation object, which do things at regular intervals
# while the simulation is running.
# This reporter creates a DCD trajectory file
# A: Christian Bale
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 100))

# This reporter prints information to the terminal
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, density=True, progress=True, 
    remainingTime=True, speed=True, totalSteps=10000, separator='\t'))

# Run the simulation itself
print('Running Production...')
simulation.step(10000)
print('Done!')
