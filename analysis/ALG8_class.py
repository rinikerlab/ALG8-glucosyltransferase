
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
                 AS_id = 513, # acceptor substrate id
                 DS_id = 514, # donor substrate id
                 rec_ids=[  26,     152,   154,   373,   369,   264,   56,    30,   153,   372,  65],              # used to calculate the hbonds
                 #            0      1      2      3      4      5      6      7     8     9     10
                 rec_name = ['D36','H162','N164','K383','H379','R274','D66', 'H40','F163','E382','E75'],  # name of the residues
                 A_sugar = ['C70','C71','C72','C73','C74','O58'],  # Acceptor sugar 
                 B_sugar = ['C52','C53','C54','C55','C56','O43'],  # Acceptor B-branch last sugar carbon
                 C_sugar = ['C58','C59','C60','C61','C52','O48'],  # Donor sugar carbon
                 D_sugar = ['C25','C26','C27','C28','C29','O7' ],  # Donor sugar carbon
                 lid_helix = [25,40],  # Residue IDs of the lid helix (T35-L50)
                 helix_loop = [150,155] ,  # Residue IDs of the helix loop (H160-Y155)
                 C_sugar_O   = ['O55','O56','O57','O59','O58'],
                 D_sugar_O   = ['O3','O4','O5','O6','O7'],
                 acceptor_O=':513@O56', acceptor_H=':513@H120',  # Acceptors
                 donor_H =':514@H47',       # Donors of DSG for bending angle and acceptor-donor angle
                 resid_offset = 10, Donor_substrate= True,  # Whether the donor is a substrate
                 ):
        # Load parameters
        self.traj = traj
        self.rec_name = rec_name 

        self.AS_id = AS_id
        self.DS_id = DS_id
        # receptor residues for hbond calculations: inter, intra
        self.rec_ids = rec_ids
        self.resid_offset = resid_offset # not used currently
         # atom selections
        self.acceptor_O = acceptor_O
        self.acceptor_H = acceptor_H
        # H is used to calculate the angle with acceptor O
        self.donor_H = donor_H
        # sugar groups for RMSD calculations
        self.Donor_substrate = Donor_substrate
        self.A_sugar = A_sugar
        self.B_sugar = B_sugar
        self.C_sugar = C_sugar
        self.C_sugar_O = C_sugar_O
        self.D_sugar = D_sugar
        self.D_sugar_O = D_sugar_O
        # lid helix and helix loop for hbond calculations
        self.lid_helix = lid_helix
        self.helix_loop = helix_loop


        # Calculation
        # self.run_analyze()



    def run_analyze(self):
        # print('test')
        # Algin the trajectory
        self.traj.superpose(mask=f':1-{self.AS_id-1}&!(@H*)')
        ## RMSD, RMSF calculation
        self.calculate_rmsd_rmsf()
        self.calculate_Hbond_onego()
        # self.calculate_hbond_protein_substrate()
        # self.calculate_hbond_protein_inter()
        # calculate other metrics
        # self.calculate_lid_helix_hbond()           
        
        
        if self.Donor_substrate:
            # self.calculate_sugar_phosphate_hbond()
            self.acceptor_donor_angle = self.calculate_angle(self.acceptor_O, f':{self.DS_id}@{self.D_sugar[0]}', self.donor_H)
            self.acceptor_donor_distance = pt.distance(self.traj,self.acceptor_O + ' '+ f':{self.DS_id}@{self.D_sugar[0]}')
            self.donor_bending_angle = self.calculate_angle(f':{self.DS_id}@P', f':{self.DS_id}@{self.D_sugar[0]}', f':{self.DS_id}@{self.D_sugar[3]}')
            self.donor_bending_angle2 = self.calculate_angle(f':{self.DS_id}@P', f':{self.DS_id}@{self.D_sugar[0]}', f':{self.DS_id}@{self.D_sugar[4]}')
        else:
            # self.calculate_sugar_phosphate_hbond()
            self.acceptor_donor_angle = None
            self.acceptor_donor_distance = None
            self.donor_bending_angle = None
        self.H40_torsion_chi1 = pt.dihedral(self.traj, f':30@N  :30@CA :30@CB :30@CG')
        self.H40_torsion_chi2 = pt.dihedral(self.traj, f':30@CA :30@CB :30@CG :30@ND1')
        
        # other interesting hbond distance
        self.Y67_O6_distance   = pt.distance(self.traj, f':57@HH :{self.DS_id}@{self.D_sugar_O[-1]}')
        self.N164_O5_distance   = pt.distance(self.traj, f':154@NE2 :{self.DS_id}@{self.D_sugar_O[3]}')
        self.N164_O6_distance   = pt.distance(self.traj, f':154@NE2 :{self.DS_id}@{self.D_sugar_O[4]}')

    def calculate_rmsd_rmsf(self):
        # Calculate the RMSD and RMSF of the trajectory
        # entities RMSD
        self.pro_rmsd = pt.rmsd(self.traj)
        self.AS_rmsd  = pt.rmsd(self.traj,mask=f':{self.AS_id}',nofit=True)
        self.DS_rmsd  = pt.rmsd(self.traj,mask=f':{self.DS_id}',nofit=True)
        # Sugar RMSD
        self.A_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.A_sugar)}',nofit=True)
        self.B_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.B_sugar)}',nofit=True)
        self.C_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.C_sugar)}',nofit=True)
        if self.Donor_substrate:
            self.D_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.DS_id}@{",".join(self.D_sugar)}',nofit=True)
        else:
            self.D_sugar_rmsd = pt.rmsd(self.traj,mask=f':{self.AS_id}@{",".join(self.D_sugar)}',nofit=True)
        self.DSP_rmsd       = pt.rmsd(self.traj,mask=f':{self.DS_id}&!(@{",".join(self.D_sugar)})',nofit=True)
        self.lid_helix_rmsd = pt.rmsd(self.traj,mask=f':{self.lid_helix[0]}-{self.lid_helix[-1]}',nofit=True)

        # RMSF 
        self.pro_rmsf = pt.rmsf(self.traj,mask='@CA') # RMSF of the enzyme
        self.DS_rmsf = pt.rmsf(self.traj,f':{self.DS_id}&!(@H*)')
        self.AS_rmsf = pt.rmsf(self.traj,f':{self.AS_id}&!(@H*)')
        
        

    def calculate_angle(self, atom1, atom2, atom3): 
        return pt.angle(self.traj, atom1 + ' ' + atom2 + ' ' + atom3)
        
    def calculate_hbond_protein_substrate(self,distance=3.0, angle=120):
        # 0: AS, 1: DP
        self.protein_substrates_bonds  = np.zeros((2,len(self.rec_ids),len(self.traj)))
        for i_rec, rec_id in enumerate(self.rec_ids):
            rec_hbond = pt.hbond(self.traj,f':{rec_id}',distance=distance, angle=angle)
            AS_intra_hbond = pt.hbond(self.traj,f':{self.AS_id}',distance=distance, angle=angle)
            DP_intra_hbond = pt.hbond(self.traj,f':{self.DS_id}',distance=distance, angle=angle)
            AS_hbond = pt.hbond(self.traj,f':{self.AS_id},{rec_id}',distance=distance, angle=angle)
            DP_hbond = pt.hbond(self.traj,f':{self.DS_id},{rec_id}',distance=distance, angle=angle)
            self.protein_substrates_bonds[0,i_rec] = AS_hbond.data[0].values - rec_hbond.data[0].values - AS_intra_hbond.data[0].values
            self.protein_substrates_bonds[1,i_rec] = DP_hbond.data[0].values - rec_hbond.data[0].values - DP_intra_hbond.data[0].values
    
    def calculate_lid_helix_hbond(self,distance=3.0, angle=120):
        lid_hbond       = pt.hbond(self.traj,f':{self.lid_helix[0]}-{self.lid_helix[-1]}',distance=distance, angle=angle)
        helix_hbond     = pt.hbond(self.traj,f':{self.helix_loop[0]}-{self.helix_loop[-1]}',distance=distance, angle=angle)
        lid_helix_hbond = pt.hbond(self.traj,f':{self.lid_helix[0]}-{self.lid_helix[-1]},{self.helix_loop[0]}-{self.helix_loop[1]}',distance=distance, angle=angle)
        self.lid_helix_hbonds = lid_helix_hbond.data[0].values - lid_hbond.data[0].values - helix_hbond.data[0].values
        
    def calculate_hbond_protein_inter(self,distance=3.0, angle=120):
        self.protein_inter_hbonds  = np.zeros((len(self.rec_ids), len(self.rec_ids),len(self.traj)))
        for i_rec, rec_id in enumerate(self.rec_ids[:-1]):
            for i_rec2, rec_id2 in enumerate(self.rec_ids[i_rec:]):
                intra_hbond1 = pt.hbond(self.traj,f':{rec_id}',distance=distance, angle=angle)
                if i_rec == i_rec2:
                    self.protein_inter_hbonds[i_rec,i_rec+i_rec2] = intra_hbond1.data[0].values
                    continue
                else:
                    intra_hbond2 = pt.hbond(self.traj,f':{rec_id2}',distance=distance, angle=angle)
                    inter_hbond = pt.hbond(self.traj,f':{rec_id},{rec_id2}',distance=distance, angle=angle)
                    self.protein_inter_hbonds[i_rec,i_rec+i_rec2] = inter_hbond.data[0].values - intra_hbond1.data[0].values - intra_hbond2.data[0].values
                    self.protein_inter_hbonds[i_rec+i_rec2,i_rec] = self.protein_inter_hbonds[i_rec,i_rec+i_rec2]  # symmetric
        

    def calculate_sugar_phosphate_hbond(self):
        if not self.Donor_substrate:
            AS_hbond = pt.hbond(self.traj,f':{self.AS_id}')
            DS_hbond = pt.hbond(self.traj,f':{self.DS_id}')
            hbond = pt.hbond(self.traj,f':{self.DS_id},{self.AS_id}')
            self.intra_hydrogen_bonds  = hbond.data[0].values - AS_hbond.data[0].values - DS_hbond.data[0].values
        else:
            DS_hbond = pt.hbond(self.traj,f':{self.DS_id}',f':{self.DS_id}')
            self.intra_hydrogen_bonds  = DS_hbond.data[0].values


    def calculate_Hbond_onego(self):
        Hbond_mask = ':'+','.join(str(id) for id in self.rec_ids)+f",{self.AS_id},{self.DS_id}"
        Hbond = pt.hbond(self.traj, Hbond_mask)
        print(Hbond.donor_acceptor[:5])
        # donor_ids = [pair.split('-')[0][3:].split('_')[0] for pair in Hbond.donor_acceptor]
        # acceptor_ids = [pair.split('-')[1][3:].split('_')[0] for pair in Hbond.donor_acceptor]

        self_hbonds = np.zeros((len(self.rec_ids)+2, len(self.rec_ids)+2, len(self.traj)))
        self.C_sugar_hbond = np.zeros((len(self.C_sugar_O), len(self.rec_ids)+2,  len(self.traj)))
        self.D_sugar_hbond = np.zeros((len(self.D_sugar_O), len(self.rec_ids)+2,  len(self.traj)))

        for i_pair, pair in enumerate(Hbond.donor_acceptor):
            donor_id = int(pair.split('-')[1][3:].split('_')[0])
            acceptor_id = int(pair.split('-')[0][3:].split('_')[0])
            donor_name = pair.split('-')[1][3:].split('_')[1]
            acceptor_name    = pair.split('-')[0][3:].split('_')[1]
            # find donor
            if donor_id in self.rec_ids:
                i_rec = self.rec_ids.index(donor_id)
            elif donor_id == self.AS_id:
                i_rec = len(self.rec_ids)
                if pair.split('-')[0][3:].split('_')[1] in self.C_sugar_O:
                    i_Csugar = self.C_sugar_O.index(pair.split('-')[0][3:].split('_')[1])
                    self.C_sugar_hbond[i_Csugar, i_rec, :] += Hbond.values[i_pair+1]

                # D sugar on acceptor
                if (not self.Donor_substrate) and (pair.split('-')[0][3:].split('_')[1] in self.D_sugar_O):
                    i_Dsugar = self.D_sugar_O.index(pair.split('-')[0][3:].split('_')[1])
                    self.D_sugar_hbond[i_Dsugar, i_rec, :] += Hbond.values[i_pair+1]

            elif donor_id == self.DS_id:
                i_rec = len(self.rec_ids)+1

            else:
                continue
            # find acceptor
            if acceptor_id in self.rec_ids:
                j_rec = self.rec_ids.index(acceptor_id)
            elif acceptor_id == self.AS_id:
                j_rec = len(self.rec_ids)
                if pair.split('-')[1][3:].split('_')[1] in self.C_sugar_O:
                    i_Csugar = self.C_sugar_O.index(pair.split('-')[1][3:].split('_')[1])
                    self.C_sugar_hbond[i_Csugar, j_rec, :] += Hbond.values[i_pair+1]

            elif acceptor_id == self.DS_id:
                j_rec = len(self.rec_ids)+1

            else:
                continue
            self_hbonds[i_rec,j_rec] += Hbond.values[i_pair+1]

            # cound C_sugar hbonds
            if (acceptor_id == self.AS_id) and (acceptor_name in self.C_sugar_O):
                i_Csugar = self.C_sugar_O.index(acceptor_name)
                self.C_sugar_hbond[i_Csugar, i_rec, :] += Hbond.values[i_pair+1]
            if (donor_id == self.AS_id) and (donor_name in self.C_sugar_O):
                i_Csugar = self.C_sugar_O.index(donor_name)
                self.C_sugar_hbond[i_Csugar, j_rec, :] += Hbond.values[i_pair+1]
            
            # count D_sugar hbonds
            if self.Donor_substrate:
                if (acceptor_id == self.DS_id) and (acceptor_name in self.D_sugar_O):
                    i_Dsugar = self.D_sugar_O.index(acceptor_name)
                    self.D_sugar_hbond[i_Dsugar, i_rec, :] += Hbond.values[i_pair+1]
                if (donor_id == self.DS_id) and (donor_name in self.D_sugar_O):
                    i_Dsugar = self.D_sugar_O.index(donor_name)
                    self.D_sugar_hbond[i_Dsugar, j_rec, :] += Hbond.values[i_pair+1]
            else: # D_sugar on acceptor
                if (acceptor_id == self.AS_id) and (acceptor_name in self.D_sugar_O):
                    i_Dsugar = self.D_sugar_O.index(acceptor_name)
                    self.D_sugar_hbond[i_Dsugar, i_rec, :] += Hbond.values[i_pair+1]
                if (donor_id == self.AS_id) and (donor_name in self.D_sugar_O):
                    i_Dsugar = self.D_sugar_O.index(donor_name)
                    self.D_sugar_hbond[i_Dsugar, j_rec, :] += Hbond.values[i_pair+1]
        

        print(self_hbonds.sum())
        # intra_hydrogen_bonds: hbonds within the protein and within the substrates, excluding protein-substrate hbonds
        self.intra_hydrogen_bonds = self_hbonds[len(self.rec_ids):,len(self.rec_ids):] + self_hbonds[len(self.rec_ids):,len(self.rec_ids):].swapaxes(0,1)
        self.intra_hydrogen_bonds[np.diag_indices(2)] /= 2  # correct double counting on diagonal
        self.protein_substrates_bonds = self_hbonds[len(self.rec_ids):,:len(self.rec_ids)] + self_hbonds[:len(self.rec_ids),len(self.rec_ids):].swapaxes(0,1)
        self.protein_inter_hbonds = self_hbonds[:len(self.rec_ids), :len(self.rec_ids)] + self_hbonds[:len(self.rec_ids), :len(self.rec_ids)].swapaxes(0,1)
        self.protein_inter_hbonds[np.diag_indices(len(self.rec_ids))] /= 2  # correct double counting on diagonal
        # print summary of hbonds
        self.Donor_intra_hbonds = pt.hbond(self.traj,f':{self.DS_id}')
        total_hbonds = self_hbonds.sum()

    def get_z_axis(self, file):
        xyz = np.loadtxt(file)
        z = xyz[:,3]
        self.D36_z = z
    
    def get_electrodensity(self, file):
        density = np.loadtxt(file)
        self.ele_density = density
    
    def get_pockets(self, file):
        pockets = np.loadtxt(file)
        self.pocket_volumes = pockets