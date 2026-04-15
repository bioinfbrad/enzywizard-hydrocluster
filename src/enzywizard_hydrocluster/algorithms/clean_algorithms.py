from __future__ import annotations
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger

from typing import Dict, List, Tuple, Optional, Union
from ..utils.clean_utils import (
    standardize_resname,
    choose_atom_altloc,
    clone_atom,
    normalize_atom_name,
    is_hydrogen_atom,
)
from ..resources.aa_resources import (
    AA3_STANDARD,
    AA3_REQUIRED_HEAVY_ATOMS,
    AA3_EXPECTED_HEAVY_ATOM_SET,
    AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT,
    BACKBONE_REQUIRED_ATOMS,
)
from Bio.PDB.Atom import Atom


def check_cleaned_structure(struct: Structure, logger: Logger) -> bool:
    model_count = len(list(struct.get_models()))
    if model_count != 1:
        logger.print("[ERROR] Structure must contain exactly one model. Please run 'enzywizard clean' first.")
        return False

    model = next(struct.get_models())
    chains = list(model.get_chains())

    if len(chains) != 1:
        logger.print("[ERROR] Structure must contain exactly one chain. Please run 'enzywizard clean' first.")
        return False

    chain = chains[0]
    if chain.id != "A":
        logger.print("[ERROR] Cleaned structure must use chain ID 'A'. Please run 'enzywizard clean' first.")
        return False

    expected_resseq = 1

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        if str(hetflag).strip():
            logger.print(f"[ERROR] Non-protein or hetero residue detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        if str(icode).strip():
            logger.print(f"[ERROR] Insertion code detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        resname = res.get_resname().strip()
        resname_std = standardize_resname(resname)
        if resname != resname_std:
            logger.print(f"[ERROR] Residue name '{resname}' at residue {resseq} is not standardized. Please run 'enzywizard clean' first.")
            return False

        if resname not in AA3_STANDARD:
            logger.print(f"[ERROR] Non-standard residue '{resname}' detected at residue {resseq}. Please run 'enzywizard clean' first.")
            return False

        if int(resseq) != expected_resseq:
            logger.print(
                f"[ERROR] Residue numbering is not continuous: expected {expected_resseq}, got {resseq}. Please run 'enzywizard clean' first."
            )
            return False

        atoms_by_name: Dict[str, List[Atom]] = {}
        for atom in res.get_atoms():
            if is_hydrogen_atom(atom):
                continue

            atom_name = normalize_atom_name(atom.get_name())
            atoms_by_name.setdefault(atom_name, []).append(atom)

        for atom_name in BACKBONE_REQUIRED_ATOMS:
            if atom_name not in atoms_by_name:
                logger.print(
                    f"[ERROR] Missing backbone atom '{atom_name}' at residue {resseq}. Please run 'enzywizard clean' first.")
                return False

        for atom_name in BACKBONE_REQUIRED_ATOMS:
            chosen = choose_atom_altloc(atoms_by_name[atom_name])
            occ = chosen.get_occupancy()
            if occ is not None and occ < 0:
                logger.print(
                    f"[ERROR] Backbone atom '{atom_name}' at residue {resseq} has invalid occupancy {occ}. Please run 'enzywizard clean' first."
                )
                return False

        expected_heavy_atom_set = AA3_EXPECTED_HEAVY_ATOM_SET.get(resname)
        if expected_heavy_atom_set is None:
            logger.print(f"[ERROR] No expected heavy atom definition for residue '{resname}' at residue {resseq}.")
            return False

        actual_heavy_atom_set = set(atoms_by_name.keys())

        if not expected_heavy_atom_set.issubset(actual_heavy_atom_set):
            missing_atom_name_list = sorted(expected_heavy_atom_set - actual_heavy_atom_set)
            logger.print(
                f"[ERROR] Missing heavy atoms at residue {resseq} ({resname}): {missing_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        allowed_heavy_atom_set = AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT[resname]
        if not actual_heavy_atom_set.issubset(allowed_heavy_atom_set):
            unexpected_atom_name_list = sorted(actual_heavy_atom_set - allowed_heavy_atom_set)
            logger.print(
                f"[ERROR] Unexpected heavy atoms at residue {resseq} ({resname}): {unexpected_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        for atom_name in actual_heavy_atom_set:
            chosen = choose_atom_altloc(atoms_by_name[atom_name])
            occ = chosen.get_occupancy()
            if occ is not None and occ < 0:
                logger.print(
                    f"[ERROR] Heavy atom '{atom_name}' at residue {resseq} has invalid occupancy {occ}. Please run 'enzywizard clean' first."
                )
                return False

        expected_resseq += 1

    return True