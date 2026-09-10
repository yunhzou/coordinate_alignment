"""Run in the pinned competitor environment: python -m unittest discover -s tests -p test_slap_edge_sweep.py."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
import slap_edge_sweep as sweep
from slapmapper.aam import SlapAAM
from slapmapper.aam._smiles import smiles2lgp,get_numbered_rxn_smiles
from slapmapper.core import SlapMapper
from rdkit import Chem


def endpoints(reaction):
    out=[]
    for side in reaction.split('>>'):
        mol=Chem.MolFromSmiles(side)
        for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
        out.append(Chem.MolToSmiles(mol))
    return out


class EdgeSweepTests(unittest.TestCase):
    def test_uncut_matches_upstream_in_both_directions_and_modes(self):
        for reaction in ('CCO>>CC=O','CCBr.O>>CCO.Br','c1ccccc1.O>>Oc1ccccc1'):
            for reverse in (False,True):
                text='>>'.join(reaction.split('>>')[::-1]) if reverse else reaction
                for binary in (False,True):
                    with self.subTest(reaction=text,binary=binary):
                        native=SlapAAM(binary=binary);native.map_smiles(text)
                        base=smiles2lgp(text,add_Hs=True)
                        adapter=SlapMapper(binary=binary)
                        adapter.get_maps(sweep.cut_graphs(base,None),break_sym_targets=[
                            i for i,z in enumerate(base[0].props['atomic numbers']) if z>1])
                        expected=sorted((r['val'],r['smiles']) for r in native.results)
                        actual=sorted((r['val'],get_numbered_rxn_smiles(text,[g.labels for g in r['lgp']]))
                                      for r in adapter.results)
                        self.assertEqual(actual,expected)

    def test_cut_rebuilds_cache_and_does_not_change_base_or_atom_counts(self):
        base=smiles2lgp('CCO>>CC=O',add_Hs=True)
        before=pickle.dumps(base)
        cut=sweep.cut_graphs(base,(0,1))
        self.assertNotIn(1,cut[0].graph[0])
        self.assertNotIn(0,cut[0].graph[1])
        self.assertEqual(len(cut[0].labels),len(base[0].labels))
        self.assertNotEqual(cut[0]._ini_wl_labels,base[0]._ini_wl_labels)
        self.assertEqual(pickle.dumps(base),before)

    def test_case_checkpoints_all_native_outputs_and_restores_endpoints(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory)
            original='CCBr.O>>CCO.Br'
            sweep.save(run/'inputs.json',[dict(index=0,input_reaction=original)])
            tasks=[]
            for reverse in (False,True):
                text='>>'.join(original.split('>>')[::-1]) if reverse else original
                bonds=Chem.AddHs(Chem.MolFromSmiles(text.split('>>')[0])).GetNumBonds()
                tasks.append(dict(slot=len(tasks),case=0,direction='P_to_R' if reverse else 'R_to_P',
                                  mode='binary',expected_cuts=bonds+1))
            sweep.save(run/'tasks.json',tasks)
            for task in tasks:
                sweep.case(argparse.Namespace(run=run,slot=task['slot']))
                path=sweep.folder(run,task)
                rows=sweep.records(path/'records.jsonl')
                self.assertEqual(len(rows),task['expected_cuts'])
                with (path/'native.bin').open('rb') as stream:
                    for row in rows:
                        self.assertEqual(row['status'],'mapped')
                        stream.seek(row['native']['offset'])
                        blob=stream.read(row['native']['length'])
                        self.assertEqual(hashlib.sha256(blob).hexdigest(),row['native']['sha256'])
                        native=pickle.loads(gzip.decompress(blob))
                        self.assertEqual(len(native),len(row['candidates']))
                        for value in row['candidates']:
                            self.assertEqual(endpoints(value['mapped_rxn']),endpoints(original))
                with (path/'records.jsonl').open('ab') as stream:stream.write(b'{"interrupted":')
                self.assertEqual(sweep.records(path/'records.jsonl'),rows)


if __name__=='__main__':unittest.main()
