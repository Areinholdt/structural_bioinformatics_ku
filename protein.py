"""
This is a module for a Quality Control (QC) analysis of Protein structure in the PDB database.
"""

# import package needed
from Bio.PDB import *
import matplotlib.pyplot as plt


class QC_analysis:
    """
    This class has an object of a quality control analyses of any PDB file.
    The class has different functions and returns values for evaluating
    regions in the protein that are potentially low quality.
    """

    def __init__(self, pdb_id):
        """
        This function connects the class with the PDB-id
        and set variable for; PDB file, R-free value, and Resolution.
        """
        self.pdb_id = pdb_id
        self.pdb_file = None
        self.rfree = None
        self.res = None

    def get_pdb(self):
        """
        This function retrieves the PDB file using the PDB-id
        and saves it to the variable.
        """
        pdb_get_me = PDBList()
        self.pdb_file = pdb_get_me.retrieve_pdb_file(self.pdb_id, file_format='pdb')

    def quality_control(self):
        """
        This function performs a quality control analysis for PDB files.
        That includes; detection of polypetide regions in the structure with potentially low quality
        here using the b-factor as an estimate for low quality.
        """
        print('The following lines is the result on the QC analysis performed on the protein',self.pdb_id, ':')

        parser = PDBParser(QUIET = True)
        structure = parser.get_structure(self.pdb_id, self.pdb_file)

        # create a dictionary to store low quality residues and chain-id
        low_qc_regions = {}

        # loop over the residues in the structure
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                                                 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']:
                        for atom in residue:
                            if atom.get_bfactor() > 80:
                                # get the residue number
                                residue_id = residue.get_id()
                                residue_number = residue_id[1]
                                # add the residue numbers and chain- id to the low_qc_regions dictionary
                                if chain.id not in low_qc_regions:
                                    low_qc_regions[chain.id] = [residue_number]
                                else:
                                    low_qc_regions[chain.id].append(residue_number)
        # remove duplicates
        for key in low_qc_regions:
            low_qc_regions[key] = list(set(low_qc_regions[key]))
        # print the low_qc_regions dictionary
        if low_qc_regions:
            for key in low_qc_regions:
                low_qc_regions[key].sort()
            for key, value in low_qc_regions.items():

                print('\033[031mThis is the residue position in the structure with high b-factors', f'chain {key}: {value}\033[0m')
        else:
            print('\033[92mThere is no detection of residue with atoms = b-factors > 80\033[0m')

    def other_qc(self):
        """
        This function extract values from the PDB file header
        that can provide valuable information about the overall quality
        of the structures.
        """
        # open and read file for extracting the r-free and resolution values.
        with open(self.pdb_file, 'r') as f:
            lines = f.readlines()
        # for r-free factors
        r_free_detect = False
        for line in lines:
            if 'FREE R VALUE  ' in line:
                r_value = line.split()[-1]
                try:
                    self.rfree = float(r_value)
                    print("R-free:", self.rfree)
                    r_free_detect = True
                except ValueError:
                    if r_value in ("NULL", "NONE", "NOT", "NaN"):
                        print('R-Free: NaN')
                        self.rfree = None
                        r_free_detect = True
                break
        if not r_free_detect:
            print('\033[031mThere is no R-Free factor in file\033[0m')

        # the same for resolution
        res_detect = False
        for line in lines[20:]:
            if 'RESOLUTION.' in line:
                res_val = line.split()[-2]
                try:
                    self.res = float(res_val)
                    print("Resolution:", self.res, 'Å')
                    res_detect = True
                except ValueError:
                    if res_val in ("NULL", "NOT", "NONE", "NaN"):
                        print('Resolution: NaN')
                        self.res = None
                        res_detect = True
                break
        if not res_detect:
            print('\033[031mThere is no Resolution values in file\033[0m')

    def rama(self):
        """
        This function provides a ramaschandran plot over the
        phi and psi angels in the polypeptides
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(self.pdb_id, self.pdb_file)
        pp = PPBuilder()
        total_length = 0
        for pp_chain in pp.build_peptides(structure, aa_only=True):
            total_length += len(pp_chain)
            phi_psi = pp_chain.get_phi_psi_list()
            phi_values = [x[0] for x in phi_psi if x[0] is not None]
            psi_values = [x[1] for x in phi_psi if x[1] is not None]
            # make a plot to detect outliers in structure
            plt.scatter(phi_values, psi_values, s=2)
            plt.xlabel("Phi-degrees)")
            plt.ylabel("Psi-degrees")
            plt.title('Ramachandran plot')
        plt.show()
        print('this is the total sequence length: ', total_length)


# Lets run the script and test it out :)
if __name__ == '__main__':
    pdb_id = '1F0P' # or any other pdb id
    # or any other PDB ID you want to test
    qc = QC_analysis(pdb_id)
    qc.get_pdb()
    qc.quality_control()
    qc.other_qc()
    qc.rama()
