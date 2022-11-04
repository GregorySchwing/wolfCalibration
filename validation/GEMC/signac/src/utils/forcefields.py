"""Utilities to load forcefields based on forcefield names."""
import os
import foyer

def get_ff_path(
    name: str = None,
):
    """Based on a forcefield name, it returns a path to that forcefield

    Parameters
    ----------
    name : str, default=None, optional
        Forcefield file name to load.
    """

    if name in ['oplsaa', 'trappe-ua']:
        ff_path = name
        return ff_path
    elif os.path.splitext(name)[1] == '.xml':
        ff_path = (
            str(os.path.dirname(os.path.abspath(__file__))) + "/../xmls/" + name
        )

        return ff_path
    else:
        raise ValueError("ERROR: This force field is not 'oplsaa' or 'trappe-ua', or does "
                         "not have a .xml after it. "
                         )

def get_molecule_path(mol2_or_smiles_input):
    """Based on a forcefield name, it returns a path to that forcefield

    Parameters
    ----------
    mol2_or_smiles_input : str,
        Whether to use a smiles string of mol2 file for the input.  The mol2 file must
        have the .mol2 extenstion or it will be read as a smiles string

    Returns
    ----------
    use_smiles : bool
        True if using a smiles string, and False if a mol2 file
    smiles_or_mol2_path_string : str
        The smiles string or the mol2 file with its path
    """
    if isinstance(mol2_or_smiles_input, str):
        if os.path.splitext(mol2_or_smiles_input)[1] == '':
            use_smiles = True
            smiles_or_mol2_path_string = mol2_or_smiles_input
            return [use_smiles, smiles_or_mol2_path_string]

        elif os.path.splitext(mol2_or_smiles_input)[1] == '.mol2':
            from src import xmls
            use_smiles = False
            smiles_or_mol2_path_string = (
                str(os.path.dirname(os.path.abspath(__file__))) + "/../molecules/" + mol2_or_smiles_input
            )
            return [use_smiles, smiles_or_mol2_path_string]
        else:
            raise TypeError("ERROR: For the get_molecule_path function,"
                            "a smiles string or a mol2 file that does not have a .mol2 "
                            "file extension was not found.")
    else:
        raise TypeError("ERROR: A string was not entered or the get_molecule_path function.")