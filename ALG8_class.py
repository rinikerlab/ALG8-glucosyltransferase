
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytraj as pt
import pytraj as pt
import scipy
# import rdkit
# from rdkit import Chem
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
                 AS_id = 513,
                 DP_id = 514,
                 rec_ids=[26,152,154,373,369,264,56, 30, 153, 372],              # used to calculate the hbonds
                 #            0      1      2      3      4      5      6      7     8     9
                 rec_name = ['D36','H162','N164','K383','H379','R274','D66', 'H40','F163','E382'],  # name of the residues
                 A_sugar = ['C70','C71','C72','C73','C74','O58'],  # Acceptor sugar carbon
                 B_sugar = ['C52','C53','C54','C55','C56','O43'],  # Donor sugar carbon
                 C_sugar = ['C58','C59','C60','C61','C52','O48'],  # Donor sugar carbon
                 donor_Cnames = ['C25','C26','C27','C28','O7'],  # Donor carbon names of DSG
                 acceptor_O=':513@O56', acceptor_H=':513@H120',  # Acceptors
                 donor_C1=':514@C25', donor_H =':514@H47', donor_C4 =':514@C28',       # Donors of DSG
                 donor_Onames = ['O8','O','O5','O7','O6'], donor_POnames = ['O','O1','O2','O3'], 
                 cat_CG = ':26@CG', cat_O1 = ':26@OD1', cat_O2 = ':26@OD2',    # Catalytic residues
                 resid_offset = 10, Donor_substrate= True):
        # Load parameters
        self.traj = traj
        self.rec_name = rec_name 
        self.AS_id = AS_id
        self.DP_id = DP_id
        self.rec_ids = rec_ids
        self.resid_offset = resid_offset
        self.donor_Onames = donor_Onames
        self.donor_POnames = donor_POnames
        self.acceptor_O = acceptor_O
        self.acceptor_H = acceptor_H
        self.donor_C1 = donor_C1
        self.donor_C4 = donor_C4
        self.donor_H = donor_H
        self.cat_CG = cat_CG
        self.cat_O1 = cat_O1
        self.Donor_substrate = Donor_substrate
        self.A_sugar = A_sugar
        self.B_sugar = B_sugar
        self.C_sugar = C_sugar
        self.donor_Cnames = donor_Cnames
        # Algin the trajectory
        self.traj.superpose(mask=f':1-{AS_id-1}&!(@H*)')
        # Placeholders 
        ## RMSD, RMSF calculation
        self.pro_rmsd = np.zeros(len(self.traj))
        self.AS_rmsd = np.zeros(len(self.traj))
        self.DS_rmsd = np.zeros(len(self.traj))
        self.intra_hydrogen_bonds = np.zeros(len(self.donor_Onames))
        self.A_sugar_rmsd = np.zeros(len(self.traj))
        self.Bsugar_rmsd = np.zeros(len(self.traj))
        self.Csugar_rmsd = np.zeros(len(self.traj))
        self.Dsugar_rmsd = np.zeros(len(self.traj))
        ## Geometric information
        # Calculation        
        self.run_analyze()

    def run_analyze(self):
        # print('test')
        self.calculate_rmsd_rmsf()
        self.calculate_hbond_protein_substrate()
        self.calculate_hbond_protein_inter()
        if self.Donor_substrate:
            self.acceptor_donor_angle = self.calculate_angle(self.acceptor_O, self.donor_C1, self.donor_H)
            self.acceptor_donor_distance = pt.distance(self.traj,self.acceptor_O + ' '+ self.donor_C1)
        self.intra_hydrogen_bonds = self.calculate_sugar_phosphate_hbond(DP_id=self.DP_id)
        self.donor_bending_angle = self.calculate_angle(f':{self.DP_id}@P', f'{self.donor_C1}', f'{self.donor_C4}')
        self.A_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.A_sugar)}',nofit=True)
        self.B_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.B_sugar)}',nofit=True)
        self.C_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.C_sugar)}',nofit=True)
        if self.Donor_substrate:
            self.D_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.DP_id}@{",".join(self.donor_Cnames)}',nofit=True)

    def calculate_rmsd_rmsf(self):
        # Calculate the RMSD and RMSF of the trajectory
        self.pro_rmsd = pt.rmsd(self.traj)
        self.AS_rmsd  = pt.rmsd(self.traj,mask=f':{self.AS_id}',nofit=True)
        self.DP_rmsd = pt.rmsd(self.traj,mask=f':{self.DP_id}',nofit=True)
        self.DPP_rmsd = pt.rmsd(self.traj,mask=f':{self.DP_id}@P,O,O1,O2,O3',nofit=True)
        self.pro_rmsf = pt.rmsf(self.traj,mask='@CA') # RMSF of the enzyme
        self.DP_rmsf = pt.rmsf(self.traj,f':{self.DP_id}&!(@H*)')
        self.AS_rmsf = pt.rmsf(self.traj,f':{self.AS_id}&!(@H*)')
        

    def calculate_angle(self, atom1, atom2, atom3): 
        return pt.angle(self.traj, atom1 + ' ' + atom2 + ' ' + atom3)
        
    def calculate_hbond_protein_substrate(self):
        # 0: AS, 1: DP
        self.protein_substrates_bonds  = np.zeros((2,len(self.rec_ids),len(self.traj)))
        for i_rec, rec_id in enumerate(self.rec_ids):
            rec_hbond = pt.hbond(self.traj,f':{rec_id}')
            AS_intra_hbond = pt.hbond(self.traj,f':{self.AS_id}')
            DP_intra_hbond = pt.hbond(self.traj,f':{self.DP_id}')
            AS_hbond = pt.hbond(self.traj,f':{self.AS_id},{rec_id}')
            DP_hbond = pt.hbond(self.traj,f':{self.DP_id},{rec_id}')
            self.protein_substrates_bonds[0,i_rec] = AS_hbond.data[0].values - rec_hbond.data[0].values - AS_intra_hbond.data[0].values
            self.protein_substrates_bonds[1,i_rec] = DP_hbond.data[0].values - rec_hbond.data[0].values - DP_intra_hbond.data[0].values
        
    def calculate_hbond_protein_inter(self):
        self.protein_inter_hbonds  = np.zeros((len(self.rec_ids), len(self.rec_ids),len(self.traj)))
        for i_rec, rec_id in enumerate(self.rec_ids[:-1]):
            for i_rec2, rec_id2 in enumerate(self.rec_ids[i_rec:]):
                intra_hbond1 = pt.hbond(self.traj,f':{rec_id}')
                if i_rec == i_rec2:
                    self.protein_inter_hbonds[i_rec,i_rec2] = intra_hbond1.data[0].values
                    continue
                else:
                    intra_hbond2 = pt.hbond(self.traj,f':{rec_id2}')
                    inter_hbond = pt.hbond(self.traj,f':{rec_id},{rec_id2}')
                    self.protein_inter_hbonds[i_rec,i_rec2] = inter_hbond.data[0].values - intra_hbond1.data[0].values - intra_hbond2.data[0].values
        

    def calculate_sugar_phosphate_hbond(self,DP_id):
        intra_hydrogen_bonds  = np.zeros((len(self.donor_Onames),len(self.traj)))
        POnames = ','.join(self.donor_POnames)
        for i_Oname, Oname in enumerate(self.donor_Onames):
            hbond = pt.hbond(self.traj,f':{DP_id}@{Oname},:{DP_id}@{POnames}')
            intra_hydrogen_bonds[i_Oname] = hbond.data[0].values
        return intra_hydrogen_bonds



