
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytraj as pt
import pytraj as pt
import scipy
# import rdkit
from rdkit import Chem
import pickle
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def sliding_average(data, window_size=5):
    """
    Compute the sliding average (moving average) of a 1D array.

    Parameters:
    - data: Input 1D array.
    - window_size: The number of elements to consider for each sliding window.

    Returns:
    - Array of sliding averages.
    """
    # Use numpy's convolve function to compute the moving average efficiently
    # The window is defined as an array of ones, divided by the window size to average
    window = np.ones(window_size) / window_size
    return np.convolve(data, window, 'valid')



class ALG_simulations():
    def __init__(self, traj,  # Traj is a pytraj trajectory
                 acceptor_O=':549@O30', acceptor_H=':549@H66',  # Acceptors
                 donor_C=':550@C', donor_H =':550@H',       # Donors
                 cat_CG = ':29@CG', cat_O1 = ':29@OD1', cat_O2 = ':29@OD2',    # Catalytic residues
                 donor_Onames = ['O8','O','O5','O7','O6'], donor_POnames = ['O','O1','O2','O3'], 
                 AS_id = 549,
                 DS_id = 550,
                 rec_ids=[29,50,59,235,313,317], # used to calculate the hbonds
                 rec_name = ['D82','E53','R112','N288','H366','R370'],
                 resid_offset = 0):
        # Load parameters
        self.traj = traj
        self.donor_Onames = donor_Onames
        self.donor_POnames = donor_POnames
        self.rec_name = rec_name 
        self.acceptor_O = acceptor_O
        self.acceptor_H = acceptor_H
        self.donor_C = donor_C
        self.donor_H = donor_H
        self.cat_CG = cat_CG
        self.cat_O1 = cat_O1
        self.cat_O2 = cat_O2
        self.AS_id = AS_id
        self.lig_id = DS_id
        self.rec_ids = rec_ids
        self.resid_offset = resid_offset
        # Algin the trajectory
        self.traj.superpose(mask=f':1-{AS_id-1}&!(@H*)')
        # Placeholders 
        ## RMSD, RMSF calculation
        self.pro_rmsd = np.zeros(len(self.traj))
        self.AS_rmsd = np.zeros(len(self.traj))
        self.DS_rmsd = np.zeros(len(self.traj))
        
        
        ## Geometric information
        self.r_cat = np.zeros(len(self.traj))
        self.D82_donor_angle = np.zeros(len(self.traj))
        self.D82_donor_angle_2D = np.zeros(len(self.traj))
        self.acceptor_donor_angle = np.zeros(len(self.traj))
        self.acceptor_donor_angle_2D = np.zeros(len(self.traj))
        self.acceptor_donor_distance = np.zeros(len(self.traj))
        self.O_hydrogen_bonds = np.zeros((len(self.donor_Onames+self.donor_POnames),len(self.rec_ids),len(self.traj)))
        self.acceptor_donor_angle = np.zeros(len(self.traj))
        self.intra_hydrogen_bonds = np.zeros(len(self.donor_Onames))
        # Calculation        
        self.run_analyze()


    def run_analyze(self):
        # print('test')
        self.pro_rmsd, self.AS_rmsd, self.DS_rmsd, self.rmsf = self.calculate_rmsd_rmsf()
        self.DS_rmsf = pt.rmsf(self.traj,f':{self.lig_id}&!(@H*)')
        self.r_cat = self.calculate_cat_acceptor_hbond(self.acceptor_H)
        self.D82_donor_angle = self.calculate_angle(self.cat_CG, self.donor_C, self.donor_H )
        self.D82_donor_angle_2D = self.calculate_angle_2D(self.cat_CG, self.donor_C, self.donor_H )
        # self.acceptor_donor_angle = self.calculate_angle(acceptor_O, donor_C, donor_H )
        self.acceptor_donor_angle = self.calculate_angle(self.acceptor_O, self.donor_C, self.donor_H)
        self.acceptor_donor_angle_2D = self.calculate_angle_2D(self.acceptor_O, self.donor_C, self.donor_H )
        self.O_hydrogen_bonds = self.calculate_hbond_protein_donor()
        self.intra_hydrogen_bonds = self.calculate_sugar_phosphate_hbond(self.lig_id)
        self.acceptor_donor_distance = pt.distance(self.traj,self.acceptor_O + ' '+ self.donor_C)

    def calculate_cat_acceptor_hbond(self,acceptor_H):
        """
        Calculate the hbond distance between the D82 and the acceptor substrate
        return the minimal D82-acceptor distance in Angstroms
        """
        r_cat1 = pt.distance(self.traj,acceptor_H + ' '+ self.cat_O1  ) # D82
        r_cat2 = pt.distance(self.traj,acceptor_H + ' '+ self.cat_O2  ) # D82
        return  np.min([r_cat1,r_cat2],axis=0)

    def calculate_sugar_phosphate_hbond(self,lig_id):
        intra_hydrogen_bonds  = np.zeros((len(self.donor_Onames),len(self.traj)))
        POnames = ','.join(self.donor_POnames)
        for i_Oname, Oname in enumerate(self.donor_Onames):
            hbond = pt.hbond(self.traj,f':{lig_id}@{Oname},:{lig_id}@{POnames}')
            intra_hydrogen_bonds[i_Oname] = hbond.data[0].values
        return intra_hydrogen_bonds

    def calculate_angle(self, atom1, atom2, atom3): 
        return pt.angle(self.traj, atom1 + ' ' + atom2 + ' ' + atom3)
    
    def calculate_angle_2D(self, atom1, atom2, atom3):  # 
        vec1 = pt.vector.vector(self.traj,atom1+' '+atom2)
        vec3 = pt.vector.vector(self.traj,atom3+' '+atom2)
        vec1 = vec1[:,:2] / np.linalg.norm(vec1[:,:2],axis=1)[:,None]
        vec3 = vec3[:,:2] / np.linalg.norm(vec3[:,:2],axis=1)[:,None]
        cos_theta = np.sum(vec1*vec3,axis=1)
        return np.arccos(cos_theta)*180/np.pi
        
    def calculate_hbond_protein_donor(self):
        O_hydrogen_bonds  = np.zeros((len(self.donor_Onames+self.donor_POnames),len(self.rec_ids),len(self.traj)))
        for i_rec, rec_id in enumerate(self.rec_ids):
            rec_hbond = pt.hbond(self.traj,f':{rec_id}')
            for i_Oname, Oname in enumerate(self.donor_Onames+self.donor_POnames):
                hbond = pt.hbond(self.traj,f':{self.lig_id}@{Oname},:{rec_id}')
                O_hydrogen_bonds[i_Oname,i_rec] = hbond.data[0].values - rec_hbond.data[0].values
        return O_hydrogen_bonds

    def calculate_rmsd_rmsf(self):
        # Calculate the RMSD and RMSF of the trajectory
        pro_rmsd = pt.rmsd(self.traj)
        AS_rmsd  = pt.rmsd(self.traj,mask=f':{self.AS_id}',nofit=True)
        DS_rmsd = pt.rmsd(self.traj,mask=f':{self.lig_id}',nofit=True)
        rmsf = pt.rmsf(self.traj,mask='@CA') # RMSF of the enzyme
        return pro_rmsd, AS_rmsd, DS_rmsd, rmsf
        

    def plot_hbond_angle(self, ax,r,theta, cmap='Blues',log=False):
        # Plot the histogram in half polar coordinates
        binsize_r = 0.3
        binsize_theta = 180/50
        r_edges     = np.arange(3,15+0.1*binsize_r,binsize_r)
        theta_edges = np.arange(0,180+0.01*binsize_theta,binsize_theta)
        
        if log:
            H, _, _ = np.histogram2d(r,theta,bins=(r_edges, theta_edges),density=False)
            H = np.log(H+1)
        else:
            H, _, _ = np.histogram2d(r,theta,bins=(r_edges, theta_edges),density=True)
        Theta, R = np.meshgrid(theta_edges/180*np.pi, r_edges)
        im = ax.pcolormesh(Theta, R, H, cmap=cmap)
        ax.set_thetamin(0)  # Start from 0 degrees
        ax.set_thetamax(180)
        return im
    




