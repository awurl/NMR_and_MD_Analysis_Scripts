# -*- coding: utf-8 -*-
"""
Spyder Editor

Fix the position of one specified atom in all trajectory frames, and adjust remaining atoms accordingly

For now only for trajectories consisting of one molecule, tested only on .gro trajectories (however other formats should be supported)
"""

#import numpy as np
#import matplotlib.pyplot as plt
import mdtraj as md

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# CONFIGURATION:

# path to directory containing the trajectories
path='./'

trajfile=path + 'traj_run6_100fs_0to1ns_res4.gro'
#topolfile=path + 'chol0.gro' #topology file needs to be specified for .xtc trajectories

input1 ='C315' # name of atom to be fixed in space
output1 = path + 'traj_run6_100fs_0to1ns_resid4_fixedC315.gro'  # output file name + file extension (all standard gromacs formats should be supported)


#------------------------------------------------------------------------------

# load trajectory
print('Loading trajectory... \n')
traj=md.load(trajfile)

topol=traj.topology

# trajectory info
print('\nTrajectory info: %s frames, %s atoms and %s residues.\n' %(traj.n_frames,traj.n_atoms,traj.n_residues))
print('Residues: ', [residue for residue in topol.residues], '\n')
# get time step during simulation points
tstep=traj.time[1]-traj.time[0]
print('Total time simulated: %s ns (step size %s ps)\n' % (traj.time[-1]/1000.0,tstep))


# coordinate adjustments

ref_atom_idx=topol.select("name "+str(input1)) 
ref_atom=topol.atom(ref_atom_idx[0])

#frame 1
ref_pos=traj.xyz[0,ref_atom.index]

for f in range(traj.n_frames-1):
    shift=ref_pos-traj.xyz[f+1,ref_atom.index]
    traj.xyz[f+1,:]=shift+traj.xyz[f+1,:]

    
traj.save(output1)