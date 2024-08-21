import numpy as np

from openeye import oechem
from openeye import oeomega
from openeye import oedocking

class OpenEyeDockingPipeline:
    def __init__(self, receptor_file):
        self.receptor_file = receptor_file

    @staticmethod
    def transform_smile(smile):
        mol_openeye = oechem.OEMol()
        oechem.OESmilesToMol(mol_openeye, smile)
        return mol_openeye

    @staticmethod
    def generate_conformers(mol, max_confs):
        rms = 0.5
        strict_stereo = False
        omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
        omegaOpts.SetMaxConfs(max_confs)
        omegaOpts.SetStrictStereo(strict_stereo)
        omegaOpts.SetRMSThreshold(rms)
        omega = oeomega.OEOmega(omegaOpts)
        error_level = oechem.OEThrow.GetLevel()
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)
        status = omega(mol)
        ret_code = omega.Build(mol)
        oechem.OEThrow.SetLevel(error_level)

    @staticmethod
    def run_docking(mol, design_final, max_confs):
        rfs = oechem.oeifstream()
        if not rfs.open(design_final):
            oechem.OEThrow.Fatal("Error")
        du = oechem.OEDesignUnit()
        if not oechem.OEReadDesignUnit(rfs, du):
            oechem.OEThrow.Fatal("Failed to read design unit")
        if not du.HasReceptor():
            raise Exception("Design unit %s does not contain a receptor" % du.GetTitle())
        dockOpts = oedocking.OEDockOptions()
        dockOpts.SetScoreMethod(oedocking.OEScoreType_PLP)
        dockOpts.SetResolution(oedocking.OESearchResolution_High)
        dock = oedocking.OEDock(dockOpts)
        dock.Initialize(du)
        dockedMol = oechem.OEMol()
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol, max_confs)
        if retCode != oedocking.OEDockingReturnCode_Success:
            print(oedocking.OEDockingReturnCodeGetName(retCode))
            raise BaseException("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))
        else:
            scores = []
            for docked_conf in dockedMol.GetConfs():
                sdtag = oedocking.OEDockMethodGetName(dockOpts.GetScoreMethod())
                oedocking.OESetSDScore(docked_conf, dock, sdtag)
                dock.AnnotatePose(docked_conf)
                score = oechem.OEGetSDData(docked_conf, sdtag)
                scores.append(float(score))
        rfs.close()
        return scores

    def run_docking_pipeline(self, smile, confs=100):
        def runner(smile):
            mol_openeye = self.transform_smile(smile)
            self.generate_conformers(mol_openeye, max_confs=confs)
            try:
                scores = self.run_docking(mol_openeye, self.receptor_file, max_confs=confs)
                return np.array(scores)
            except BaseException as e:
                return None

        scores_array = runner(smile)
        return scores_array