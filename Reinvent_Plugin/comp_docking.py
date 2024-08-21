"""
OpenEye Docking
Example toml config
[component.OpenEyeDocking]
[[component.OpenEyeDocking.endpoint]]
name = "OpenEye Docking"
params.receptor_file = "path_to_receptor_file"
params.maxconfs = 100
"""

from __future__ import annotations

__all__ = ["OpenEyeDocking"]

import os
import numpy as np
from pydantic import Field
from pydantic.dataclasses import dataclass
from typing import List, Optional
import logging

from .docking.fix_protonate import prepare_small_mol
from .docking.docking_oedock import OpenEyeDockingPipeline
from ..component_results import ComponentResults
from ..add_tag import add_tag

logger = logging.getLogger(__name__)

@add_tag("__parameters")
@dataclass
class Parameters:
    """Parameters for the docking component

    Note that all parameters are always lists because components can have
    multiple endpoints and so all the parameters from each endpoint are
    collected into a list. This is also true in cases where there is only one
    endpoint.
    """

    receptor_file: List[str]
    maxconfs: List[int]

@add_tag("__component")
class OpenEyeDocking:
    def __init__(self, params: Parameters):
        self.receptor_file = params.receptor_file[0]
        self.maxconfs = params.maxconfs[0]
        self.docking = OpenEyeDockingPipeline(receptor_file=self.receptor_file)

    # def __call__(self, smilies: List[str]) -> ComponentResults:
    #     smile = smilies[0]
    #     scores = self.docking.run_docking_pipeline(smile, confs=self.maxconfs)
    #
    #     return ComponentResults([scores])

    # # FIXME: Add a Try Except block to catch errors adding a score very high (Nonetype)
    # def __call__(self, smilies: List[str]) -> ComponentResults:
    #     all_scores = []
    #     for smile in smilies:
    #         scores = self.docking.run_docking_pipeline(smile, confs=self.maxconfs)  # Numpy Array
    #         scores_mean = np.mean(scores)
    #         all_scores.append(scores_mean)
    #
    #     return ComponentResults([all_scores])

    def __call__(self, smilies: List[str]) -> ComponentResults:
        all_scores = []
        for smile in smilies:
            smile_protonated = prepare_small_mol(smile, gen_3d=False, ID=None, protonate=True)
            scores = self.docking.run_docking_pipeline(smile_protonated, confs=self.maxconfs)   # Numpy Array

            if scores is None:
                final_score = 1000
            else:
                final_score = np.mean(scores)       # Numeric

            all_scores.append(final_score)

        return ComponentResults([all_scores])