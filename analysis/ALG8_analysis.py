# conda activate shu_pytraj
# clean

import sys
import os

# Save original path
original_sys_path = sys.path.copy()

# Remove custom entries not in original sys.path (like PYTHONPATH entries) to avoid numpy conflict
sys.path = [p for p in sys.path if "amber" not in p]
from ALG8_class import *
import argparse

"""
prefix = 'ALG8_pretransfer_DSG'
prmtop = path + f'{prefix}_stripped.prmtop'
ncs    = [path + f'{prefix}_run{run}.nc' for run in range(1,N_run+1)]
trajs  = [pt.load(nc, prmtop) for i_nc, nc in enumerate(ncs)]
state  = 'pre-DSG'  # 'pre-DSG', 'pre-DSM, 'post'
take_frame = 4000
"""

# parse command line arguments
parser = argparse.ArgumentParser(description='Run ALG8 analysis with specified parameters.')
parser.add_argument('--MD_path', type=str, default='/localhome/shuchen/projects/ALG8/MD_simulations/', help='Path to MD simulation files')
parser.add_argument('--prefix', type=str, required=True, help='Prefix for the trajectory files')
parser.add_argument('--N_run', type=int, required=True, help='Number of runs to analyze')
parser.add_argument('--state', type=str, default='pre-DSG', help='State of the system: pre-DSG, pre-DSM, post')
parser.add_argument('--take_frame', type=int, default=4000, help='Number of frames to take from each trajectory')
parser.add_argument('--ob_path', type=str, default='./objects/', help='Path to save analysis objects')

# example usage:
# python analysis/ALG8_analysis.py --prefix ALG8_pretransfer_DSG --N_run 5 --DSG_donor_C1 ':514@C25' --DSG_donor_C4 ':514@C28' 
# --DSG_donor_H ':514@H47' --donor_POnames O O1 O2 O8 --donor_Onames O3 O4 O5 O6 O7 --donor_Cnames C25 C26 C27 C28 O7

args    = parser.parse_args()
prefix  = args.prefix
MD_path = args.MD_path
N_run   = args.N_run
state   = args.state
take_frame = args.take_frame
ob_path    = args.ob_path
if not os.path.exists(ob_path): 
    os.makedirs(ob_path)

# determine donor atom selections based on the state
rec_name = ['D36','H162','N164','K383','H379','R274','D66', 'H40','F163','E382','E75'],  # name of the residues
rec_ids=[  26,     152,   154,   373,   369,   264,   56,    30,   153,   372,  65]

if 'pre' in state:
    A_sugar  = ['C70','C71','C72','C73','C74','O58']  # Acceptor A-branch last sugar carbon
    C_sugar  = ['C52','C53','C54','C55','C56','O43']  # Acceptor B-branch last sugar carbon
    B_sugar  = ['C58','C59','C60','C61','C52','O48']  # Acceptor C-branch last sugar carbon
    D_sugar  = ['C25','C26','C27','C28','C29','O7' ]  # Donor sugar carbon
    donor_H  =':514@H47' if state == 'pre-DSG' else ':514@H46'
    Donor_is_substrate = True
    C_sugar_O   = ['O55','O56','O57','O59','O58']
    D_sugar_O   = ['O3','O4','O5','O6','O7']

else:  # post
    A_sugar = ['C95','C96','C97','C98','C99','O65']  # Product A-branch last sugar carbon
    C_sugar = ['C71','C72','C73','C74','C75','O45']  # Product B-branch last sugar carbon
    B_sugar = ['C83','C84','C85','C86','C87','O55']  # Product C-branch last sugar carbon
    D_sugar = ['C101','C102','C103','C104','C105', 'O70' ]  # Product new sugar carbon
    Donor_is_substrate = False
    C_sugar_O   = ['O62','O63','O64','O66','O65']
    D_sugar_O   = ['O67','O68','69','O71','O70']
    donor_H     = None

prmtop = MD_path + f'{prefix}_stripped.prmtop'
ncs    = [MD_path + f'{prefix}_run{run}.nc' for run in range(1,N_run+1)]
print(f"prmtop: {prmtop}, ncs: {ncs}")

trajs  = [pt.load(nc, prmtop) for i_nc, nc in enumerate(ncs)]

object_list = [ALG_simulations(traj[:take_frame], 
                              donor_H = donor_H,
                              A_sugar = A_sugar,
                              B_sugar = B_sugar,
                              C_sugar = C_sugar,
                              D_sugar = D_sugar,
                                C_sugar_O = C_sugar_O, 
                                D_sugar_O = D_sugar_O,
                              Donor_substrate=Donor_is_substrate) for traj in trajs]
# state: 'pre-DSG', 'pre-DSM, 'post' 

for i_obj, obj in enumerate(object_list):
    obj_name = ob_path + f'{prefix}_run{i_obj+1}_analysis_obj.pkl'
    # Run analysis
    obj.run_analyze()
    obj.traj = None # remove traj to save space
    
    
    try:
        file_electro = MD_path + f'z_axis/{prefix}_{i_obj+1}_lipid_edens.dat'
        obj.get_electrodensity(file_electro)
        file_D36_z = MD_path + f'z_axis/{prefix}_{i_obj+1}_D36.dat'
        obj.get_z_axis(file_D36_z)
    except FileNotFoundError:
        print(f"File not found: {file_electro}")
    pickle.dump(obj, open(obj_name,'wb'))
    print(f'Object saved as {obj_name}')

# Save the analysis results




