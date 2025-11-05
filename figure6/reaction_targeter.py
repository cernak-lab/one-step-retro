#!/usr/bin/env python

import pickle
import time
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops

RDLogger.DisableLog("rdApp.*")
from typing import Any
import copy

__all__ = ["ReactionTargeter"]


class ReactionTargeter:
    # ---------- Lifecycle ----------
    def __init__(self, commercial_index: dict, experiment_id=1):
        self.commercial_index = commercial_index
        seen = set[Any]()
        target_substructures = []
        for ts1 in [
            "[c:1]",
            "[C:1]",
            "[N:1]",
            "[n:1]",
            "[O:1]",
            "[o:1]",
            "[S:1]",
            "[P:1]",
        ]:
            for ts2 in [
                "[c:2]",
                "[C:2]",
                "[N:2]",
                "[n:2]",
                "[O:2]",
                "[o:2]",
                "[S:2]",
                "[P:2]",
            ]:
                its = tuple[str, ...](sorted([ts1[1], ts2[1]]))
                if its in seen:
                    continue
                seen.add(its)
                ts = f"{ts1}-{ts2}"
                target_substructures.append(ts)
        self.target_substructures = target_substructures

        self.disconnection_map = [1, 2]
        self.substructure_atom_map = [1, 2]
        self.building_blocks_a = [
            "[O:3]",
            "[N:3]",
            "[C:3](=O)[O]",
            "[Cl:3]",
            "[Br:3]",
            "[I:3]",
            "[B:3](O)[OH]",
        ]
        self.building_blocks_b = [
            "[O:4]",
            "[N:4]",
            "[C:4](=O)[O]",
            "[Cl:4]",
            "[Br:4]",
            "[I:4]",
            "[B:4](O)[OH]",
        ]

        self.building_blocks_labels_a = [
            "hydrogen",
            "alcohol",
            "amine",
            "acid",
            "chloride",
            "bromide",
            "iodide",
            "boronate",
        ]
        self.building_blocks_labels_b = [
            "hydrogen",
            "alcohol",
            "amine",
            "acid",
            "chloride",
            "bromide",
            "iodide",
            "boronate",
        ]
        self.output_path = f"/home/bmahjour/nprt/out/{experiment_id}.xlsx"

    @classmethod
    def from_commercial_index_file(cls, path: str):
        with open(path, "rb") as f:
            commercial_index = pickle.load(f)
        return cls(commercial_index)

    # ---------- Fast availability index ----------
    @staticmethod
    def build_commercial_index(commercial_df):
        idx = {}
        for _, row in commercial_df.iterrows():
            s = row["sanitized_smiles"]
            if s not in idx:
                idx[s] = row.to_dict()
        return idx

    # ---------- Validation ----------
    @staticmethod
    def check_input_target_molecule(target_molecule_smiles):
        mol = Chem.MolFromSmiles(target_molecule_smiles)
        if mol is None:
            return None, "invalid target molecule smiles"
        return mol, ""

    @staticmethod
    def check_input_target_substructure(target_substructure_smarts):
        q = Chem.MolFromSmarts(target_substructure_smarts)
        if q is None:
            return None, "invalid target substructure smarts"
        if any(a.GetAtomMapNum() == 0 for a in q.GetAtoms()):
            return None, "target substructure missing atom number label"
        return q, ""

    @staticmethod
    def check_number_of_substructure_hits(target_mol, substructure_mol):
        hits = target_mol.GetSubstructMatches(substructure_mol)
        return (hits, "") if hits else (None, "no substructure hits in target")

    # ---------- Mapping ----------
    @staticmethod
    def copy_and_map_atoms_to_target(
        target_mol, hit_atom_numbers, substructure_atom_map
    ):
        copy = Chem.Mol(target_mol)
        for i, atom_idx in enumerate(hit_atom_numbers):
            copy.GetAtomWithIdx(atom_idx).SetAtomMapNum(int(substructure_atom_map[i]))
        return copy, Chem.MolToSmiles(copy)

    # ---------- Cutting with FragmentOnBonds ----------
    @staticmethod
    def _mapnum_to_index(mol):
        out = {}
        for a in mol.GetAtoms():
            m = a.GetAtomMapNum()
            if m:
                out[m] = a.GetIdx()
        return out

    def disconnect_synthon(self, mapped_sm, dis_labels):
        mol = Chem.MolFromSmiles(mapped_sm)
        amap = self._mapnum_to_index(mol)
        a_idx = amap[int(dis_labels[0])]
        b_idx = amap[int(dis_labels[1])]
        bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
        if bond is None:
            raise ValueError("Disconnection bond not found between mapped atoms.")
        frag = rdmolops.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
        parts = Chem.MolToSmiles(frag).split(".")
        synthon_a = max((p for p in parts if f":{dis_labels[0]}]" in p), key=len)
        synthon_b = max((p for p in parts if f":{dis_labels[1]}]" in p), key=len)
        synthon_a = synthon_a.replace("[n:1]", "[nH:1]")
        synthon_b = synthon_b.replace("[n:2]", "[nH:2]")
        return synthon_a, synthon_b

    # ---------- Map stripping ----------
    @staticmethod
    def remove_atom_labels(smiles_with_maps):
        m = Chem.MolFromSmiles(smiles_with_maps, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        return Chem.MolToSmiles(m)

    # ---------- Fast enumeration & availability ----------
    def add_fg_to_synthon(
        self,
        synthon_smiles,
        connection_label,
        building_block_smiles,
        building_block_label,
    ):
        syn_m = Chem.MolFromSmiles(synthon_smiles)
        bb_m = Chem.MolFromSmiles(building_block_smiles)
        comb = Chem.CombineMols(syn_m, bb_m)
        rw = Chem.RWMol(comb)

        def find_idx_by_map(mol, mapnum):
            for a in mol.GetAtoms():
                if a.GetAtomMapNum() == int(mapnum):
                    return a.GetIdx()
            return None

        syn_idx = find_idx_by_map(syn_m, connection_label)
        bb_idx = find_idx_by_map(bb_m, building_block_label)
        if syn_idx is None or bb_idx is None:
            raise ValueError("Could not find mapped atoms to connect.")

        offset = syn_m.GetNumAtoms()
        rw.AddBond(int(syn_idx), int(bb_idx + offset), Chem.rdchem.BondType.SINGLE)
        atom = rw.GetAtomWithIdx(int(syn_idx))
        atom.SetNumExplicitHs(0)
        out = rw.GetMol()
        return Chem.MolToSmiles(out)

    def check_synthon_for_commercial_analogs(
        self, synthon_smiles, connection_label, building_blocks, building_block_label
    ):
        candidates = [synthon_smiles]
        for bb in building_blocks:
            candidates.append(
                self.add_fg_to_synthon(
                    synthon_smiles, int(connection_label), bb, int(building_block_label)
                )
            )
        stripped = [self.remove_atom_labels(s) for s in candidates]
        return [self.commercial_index.get(s) for s in stripped]

    # ---------- Main ----------
    def target_bond_breaker(
        self,
        target_molecule_smiles,
        target_substructure_smarts,
        debug=False,
    ):
        t0 = time.time()

        building_blocks_a = self.building_blocks_a
        building_blocks_b = self.building_blocks_b
        substructure_atom_map = self.substructure_atom_map
        disconnection_map = self.disconnection_map
        bb_map = [3, 4]

        empty_cols = [
            "input_target_molecule_smiles",
            "mapped_input_target_molecule_smiles",
            "target_substructure_smarts",
            "substructure_atom_map",
            "disconnection_map",
            "synthon_a",
            "synthon_a_label",
            "synthon_b",
            "synthon_b_label",
            "synthon_a_building_block_class",
            "synthon_a_building_block_label",
            "synthon_b_building_block_class",
            "synthon_b_building_block_label",
            "synthon_a_building_block_smiles",
            "synthon_a_building_block_cas",
            "synthon_a_building_block_price",
            "synthon_a_building_block_material",
            "synthon_b_building_block_smiles",
            "synthon_b_building_block_cas",
            "synthon_b_building_block_price",
            "synthon_b_building_block_material",
        ]
        rows = {k: [] for k in empty_cols}

        tmol, msg = self.check_input_target_molecule(target_molecule_smiles)
        if not tmol:
            if debug:
                print(msg)
            return 0, time.time() - t0, pd.DataFrame(rows)

        qmol, msg = self.check_input_target_substructure(target_substructure_smarts)
        if not qmol:
            if debug:
                print(msg)
            return 0, time.time() - t0, pd.DataFrame(rows)

        hits, msg = self.check_number_of_substructure_hits(tmol, qmol)
        if not hits:
            if debug:
                print(msg)
            return 0, time.time() - t0, pd.DataFrame(rows)

        mapped_smiles_list = []
        for hit in hits:
            _, mapped_smi = self.copy_and_map_atoms_to_target(
                tmol, hit, substructure_atom_map
            )
            mapped_smiles_list.append(mapped_smi)

        for mapped_sm in mapped_smiles_list:
            try:
                syn_a, syn_b = self.disconnect_synthon(mapped_sm, disconnection_map)
            except Exception as e:
                if debug:
                    print(f"disconnect failed: {e}")
                continue

            syn_a_hits = self.check_synthon_for_commercial_analogs(
                syn_a, disconnection_map[0], building_blocks_a, bb_map[0]
            )
            if syn_a == syn_b:
                syn_b_hits = syn_a_hits
            else:
                syn_b_hits = self.check_synthon_for_commercial_analogs(
                    syn_b, disconnection_map[1], building_blocks_b, bb_map[1]
                )
            if len(syn_a_hits) == 0 or len(syn_b_hits) == 0:
                continue

            for i, sa in enumerate(syn_a_hits):
                if sa is None:
                    continue
                for j, sb in enumerate(syn_b_hits):
                    if sb is None:
                        continue
                    rows["input_target_molecule_smiles"].append(target_molecule_smiles)
                    rows["mapped_input_target_molecule_smiles"].append(mapped_sm)
                    rows["target_substructure_smarts"].append(
                        target_substructure_smarts
                    )
                    rows["substructure_atom_map"].append(substructure_atom_map)
                    rows["disconnection_map"].append(disconnection_map)
                    rows["synthon_a"].append(syn_a)
                    rows["synthon_b"].append(syn_b)
                    rows["synthon_a_label"].append(disconnection_map[0])
                    rows["synthon_b_label"].append(disconnection_map[1])
                    rows["synthon_a_building_block_smiles"].append(
                        sa["sanitized_smiles"]
                    )
                    rows["synthon_a_building_block_cas"].append(sa.get("cas"))
                    rows["synthon_a_building_block_price"].append(sa.get("price"))
                    rows["synthon_a_building_block_material"].append(sa.get("material"))
                    rows["synthon_a_building_block_class"].append(
                        self.building_blocks_labels_a[i]
                    )
                    rows["synthon_a_building_block_label"].append(bb_map[0])
                    rows["synthon_b_building_block_smiles"].append(
                        sb["sanitized_smiles"]
                    )
                    rows["synthon_b_building_block_cas"].append(sb.get("cas"))
                    rows["synthon_b_building_block_price"].append(sb.get("price"))
                    rows["synthon_b_building_block_material"].append(sb.get("material"))
                    rows["synthon_b_building_block_class"].append(
                        self.building_blocks_labels_b[j]
                    )
                    rows["synthon_b_building_block_label"].append(bb_map[1])

        df = pd.DataFrame(rows)
        dt = time.time() - t0
        return len(df), dt, df

    def run_1_ablation(self, target_molecule_smiles):
        all_hits_1 = []
        seen_rxns = set[Any]()
        target_substructures = self.target_substructures
        for target_substructure in target_substructures:
            tmol_smiles = copy.deepcopy(target_molecule_smiles)
            _, _, df_1 = self.target_bond_breaker(
                tmol_smiles,
                target_substructure,
            )
            for i, k in df_1.iterrows():
                if (
                    k["synthon_a_building_block_smiles"] is None
                    and k["synthon_b_building_block_smiles"] is None
                ):
                    continue
                trxn = f'{k["synthon_a_building_block_smiles"]}.{k["synthon_b_building_block_smiles"]}>>{k["input_target_molecule_smiles"]}'
                if trxn in seen_rxns:
                    continue
                seen_rxns.add(trxn)
                all_hits_1.append(k.to_dict())
        return pd.DataFrame(all_hits_1)

    def run(self, data):
        return self.target_bond_breaker(
            data["target_molecule"],
            data["target_substructure"],
        )


def run(data):
    rt = ReactionTargeter.from_commercial_index_file("commercial_index.pkl")
    return rt.run(data)
