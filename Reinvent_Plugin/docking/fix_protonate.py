import datamol as dm

from rdkit import Chem
from rdkit.Chem import AllChem

from openeye import oechem
from openeye import oequacpac

class Fix_and_Protonate:
    def __init__(self, in_smile):
        # Shutdown the OEChem error handler
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Warning)

        # Conver the smile to mol
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, in_smile)

        self.mol = mol

    def cansmi(self, mol, isomeric, kekule):
        """
        Generate a canonical SMILES representation of a molecule.

        Args:
            mol (oechem.OEMol): The molecule to generate SMILES for.
            isomeric (bool): Flag indicating whether to include isomeric information in the SMILES.
            kekule (bool): Flag indicating whether to kekulize the molecule before generating SMILES.

        Returns:
            str: The canonical SMILES representation of the molecule.
        """
        oechem.OEFindRingAtomsAndBonds(mol)
        oechem.OEAssignAromaticFlags(mol, oechem.OEAroModel_OpenEye)
        smiflag = oechem.OESMILESFlag_Canonical
        if isomeric:
            smiflag |= oechem.OESMILESFlag_ISOMERIC

        if kekule:
            for bond in mol.GetBonds(oechem.OEIsAromaticBond()):
                bond.SetIntType(5)
            oechem.OECanonicalOrderAtoms(mol)
            oechem.OECanonicalOrderBonds(mol)
            oechem.OEClearAromaticFlags(mol)
            oechem.OEKekulize(mol)

        # Strip Salt
        oechem.OEDeleteEverythingExceptTheFirstLargestComponent(mol)

        smi = oechem.OECreateSmiString(mol, smiflag)
        return smi

    def protomer(self, mol):
        opts = oequacpac.OEMultistatepKaModelOptions()
        multistatepka = oequacpac.OEMultistatepKaModel(mol, opts)

        # Generate all protomers and add it to a list
        if (multistatepka.GenerateMicrostates()):
            all_protomers = []

            for a in multistatepka.GetMicrostates():
                all_protomers.append(a)

            return all_protomers

    def standardize_small_mol(self, protonate=True):
        # Standardize the smile
        smi_iso = self.cansmi(self.mol, True, True)

        # Convert ISO SMILE to Mol
        mol_prot = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol_prot, smi_iso)

        # Protonate
        if protonate == True:
            all_protomers = self.protomer(mol_prot)
            all_protomers_smile = [oechem.OECreateSmiString(a) for a in all_protomers]
            return mol_prot, all_protomers_smile[0]

        else:
            return mol_prot, smi_iso


def _datamol_preprocess(smile, reionize: bool, uncharge: bool):
    mol = dm.to_mol(smile)
    mol = dm.fix_mol(mol)
    mol = dm.sanitize_mol(mol, sanifix=True, charge_neutral=False)
    mol = dm.standardize_mol(
        mol,
        disconnect_metals=False,
        normalize=True,
        reionize=reionize,
        uncharge=uncharge,
        stereo=True,
    )
    return dm.to_smiles(mol, kekulize=True)

def prepare_small_mol(smile:str, gen_3d=True, ID=None, protonate=True):
    prep = Fix_and_Protonate(smile)
    _, smile_prep = prep.standardize_small_mol(protonate=protonate)

    def generate_3d(smile):
        mol = dm.to_mol(smile)
        mol_h = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol_h)

        if ID is not None:
            mol_h.SetProp("_Name", ID)

        return mol_h

    if gen_3d == True:
        mol_prep_3d = generate_3d(smile_prep)
        return mol_prep_3d, smile_prep
    else:
        return smile_prep

def prepare_for_searching(smile):
    prep = Fix_and_Protonate(smile)
    _, smile_prep = prep.standardize_small_mol(protonate=False)
    smile_prep_prep = _datamol_preprocess(smile_prep, reionize=True, uncharge=False)
    return smile_prep_prep