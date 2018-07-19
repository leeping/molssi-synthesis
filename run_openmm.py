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

# Use the CPU platform
platform = mm.Platform.getPlatformByName('CPU')

### If you want to add any forces to your System or modify any
### of the existing forces, you should do it here - after the
### System has been created, but before the Simulation is created.

# Create a Simulation object by putting together the objects above
simulation = app.Simulation(pdb.topology, system, integrator, platform)

# Set positions in the Simulation object
simulation.context.setPositions(pdb.positions)

# Minimize the energy of the system (intentionally doing a rough minimization)
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=20, tolerance=100)

# Initialize the random velocities of the system from a Maxwell-Boltzmann distribution
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# *Incredibly Short* NPT equilibration (5 ps)
# First, add reporters to the simulation object, which do things at regular intervals
# while the simulation is running. This reporter prints information to the terminal.
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True,
    potentialEnergy=True, temperature=True, density=True, progress=True,
    remainingTime=True, speed=True, totalSteps=2500, separator='\t'))
# Next, run the equilibration simulation itself.
print('Running Equilibration...')
simulation.step(2500)

# Production (40 ps)
# Before doing anything, remember to clear the previous reporters. In the code below
# the first reporter creates a DCD trajectory file. The second reporter is otherwise
# very similar to the one above in the equilibration section. The only difference --
# other than ouput frequency -- is that the totalSteps parameter is modified to be
# the number of production steps + equilibration steps. Run the code to see why...
# A: Christian Bale
simulation.reporters.clear()
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 100))
simulation.reporters.append(app.StateDataReporter(stdout, 500, step=True, 
    potentialEnergy=True, temperature=True, density=True, progress=True, 
    remainingTime=True, speed=True, totalSteps=22500, separator='\t'))
# Run the production simulation!
print('Running Production...')
simulation.step(20000)
print('Done!')
