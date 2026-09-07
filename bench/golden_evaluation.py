"""Reference-blind inputs and symmetry-aware scoring of saved AAM graphs.

Coverage queries use Schreier-Sims transversal choices, not group-element or
bijection enumeration. Reference labels are only consumed after search.
"""
import ast
from collections import defaultdict
from functools import lru_cache
import time

import numpy as np
import pynauty
from rdkit import Chem
from sympy.combinatorics import Permutation, PermutationGroup

from rxn_core.domain import AAMProblem, MolecularEndpoint
from rxn_core.family_query import query_path
from rxn_core.search_graph import frozen_value


def prepare(mapped_reaction):
    endpoints, features, labels = [], [], []
    for side, smiles in zip(('R', 'P'), mapped_reaction.split('>>')):
        mol = Chem.MolFromSmiles(smiles)
        maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.MolToSmiles(mol, canonical=True)
        order = ast.literal_eval(mol.GetProp('_smilesAtomOutputOrder'))
        mol = Chem.AddHs(Chem.RenumberAtoms(mol, order))
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
        hindex = {a:i for i,a in enumerate(heavy)}
        label = {maps[old]: new for new,old in enumerate(order) if maps[old]}
        labels.append(label)
        matrix = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))
        bonds = []
        for bond in mol.GetBonds():
            a,b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            matrix[a,b] = matrix[b,a] = bond.GetBondTypeAsDouble()
            if a in hindex and b in hindex:
                bonds.append((hindex[a],hindex[b],(bond.GetBondTypeAsDouble(),str(bond.GetStereo()))))
        colors = []
        for i in heavy:
            atom = mol.GetAtomWithIdx(i)
            colors.append((atom.GetSymbol(),atom.GetFormalCharge(),atom.GetIsotope(),
                           atom.GetTotalNumHs(includeNeighbors=True),
                           atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else ''))
        features.append(dict(heavy=heavy, colors=colors, bonds=bonds))
        endpoints.append(MolecularEndpoint(tuple(a.GetSymbol() for a in mol.GetAtoms()),
            np.zeros((mol.GetNumAtoms(),3)), matrix, side,
            metadata={'unmapped_smiles':Chem.MolToSmiles(mol), 'scoring_features':features[-1]}))
    reference = {labels[0][k]:labels[1][k] for k in labels[0].keys() & labels[1].keys()}
    return AAMProblem(*endpoints), features, reference


def colored_graph(features, mapping=None):
    adjacency, colors = defaultdict(set), defaultdict(set)
    offsets, count = [], 0
    for side, feature in enumerate(features):
        offsets.append(count)
        for color in feature['colors']:
            colors[('atom',side,tuple(color))].add(count); count += 1
    def link(a,b,color):
        nonlocal count
        colors[color].add(count)
        adjacency[count].update((a,b)); adjacency[a].add(count);adjacency[b].add(count)
        count += 1
    for side,feature in enumerate(features):
        for a,b,color in feature['bonds']:
            link(offsets[side]+a,offsets[side]+b,('bond',side,tuple(color)))
    if mapping is not None:
        for a,b in sorted(mapping.items()):
            link(a,offsets[1]+b,('mapping',))
    return pynauty.Graph(count, adjacency_dict={i:sorted(adjacency[i]) for i in range(count)},
                         vertex_coloring=[colors[k] for k in sorted(colors,key=repr)])


def project(mapping, features):
    source = {a:i for i,a in enumerate(features[0]['heavy'])}
    target = {a:i for i,a in enumerate(features[1]['heavy'])}
    return {source[a]:target[b] for a,b in dict(mapping).items() if a in source and b in target}


def endpoint_generators(feature):
    n = len(feature['heavy'])
    return tuple(tuple(g[:n]) for g in pynauty.autgrp(colored_graph([feature]))[0])


def exact_action(g, feature):
    colors = feature['colors']
    if any(colors[i] != colors[g[i]] for i in range(len(g))):
        return False
    bonds = {(min(a,b),max(a,b)):tuple(c) for a,b,c in feature['bonds']}
    return all(bonds.get(tuple(sorted((g[a],g[b])))) == c for (a,b),c in bonds.items())


@lru_cache(maxsize=1024)
def transversal_factors(generators, degree):
    if not generators:
        return ()
    group = PermutationGroup([Permutation(list(g),size=degree) for g in generators])
    # An element acts as u0(u1(...uk(x))). These are subgroup transversals,
    # not all group elements; only their symbolic choice variables are stored.
    return tuple(tuple(tuple(p(i) for i in range(degree)) for p in level.values())
                 for level in reversed(group.basic_transversals) if len(level)>1)


def symbolic_path_query(mapping, groups, features, reference, timeout_ms, *, finite_domain=False,
                        independent_singleton_domains=None):
    """Query group actions, optionally with independently validated singleton choices.

    This is NOT a complete arbitrary-domain verifier. The optional domains
    must not replace correlated multi-atom fragment assignments. The diagnostic
    caller restricts them to trailing one-atom fragments and checks its witness.
    """
    import z3
    nr,np_ = len(features[0]['heavy']),len(features[1]['heavy'])
    width=max(nr,np_,1).bit_length()
    sort=z3.BitVecSort(width) if finite_domain else z3.IntSort()
    number=(lambda n:z3.BitVecVal(n,width)) if finite_domain else z3.IntVal
    variable=(lambda name:z3.BitVec(name,width)) if finite_domain else z3.Int
    solver = z3.SolverFor('QF_ABV') if finite_domain else z3.Solver()
    solver.set(timeout=max(1,int(timeout_ms)))
    choices, serial = [], 0
    def table(images, sentinel):
        array = z3.K(sort,number(sentinel))
        for i,value in enumerate(images):
            array = z3.Store(array,i,value)
        return array
    def act(values, generators, degree, label):
        nonlocal serial
        selected = []
        for factor in transversal_factors(tuple(generators),degree):
            choice = variable(f'c{serial}');serial += 1
            if finite_domain:solver.add(z3.ULT(choice,number(len(factor))))
            else:solver.add(choice>=0,choice<len(factor))
            arr = table(factor[-1],degree)
            for index in reversed(range(len(factor)-1)):
                arr = z3.If(choice==index,table(factor[index],degree),arr)
            result = []
            for value in values:
                transported = variable(f'x{serial}');serial += 1
                solver.add(transported == z3.Select(arr,value));result.append(transported)
            values = result
            choices.append((label,choice,factor)); selected.append(choice)
        return values
    values = act([number(i) for i in range(nr)], endpoint_generators(features[0]),nr,'source_equivalence')
    selected_base=dict(mapping)
    domain_variables={}
    for atom,domain in (independent_singleton_domains or {}).items():
        value=variable(f'domain_{atom}')
        solver.add(z3.Or(*(value==number(p) for p in domain)))
        selected_base[atom]=value;domain_variables[atom]=value
    if domain_variables:
        solver.add(z3.Distinct(*selected_base.values()))
    base = table([selected_base.get(i,np_) for i in range(nr)],np_)
    values = [z3.Select(base,x) for x in values]
    # Recorded chronological actions compose left-to-right; apply last first.
    for index in reversed(range(len(groups))):
        values = act(values,groups[index],np_,f'path_{index}')
    values = act(values,endpoint_generators(features[1]),np_,'target_equivalence')
    complete_reference=len(set(reference.values()))==np_
    solver.add(*(value==reference.get(i,np_) for i,value in enumerate(values)
                 if complete_reference or i in reference))
    status = solver.check()
    if status == z3.sat:
        model = solver.model()
        actions = [dict(stage=label,transversal_index=model.eval(choice).as_long(),
                                 permutation=list(factor[model.eval(choice).as_long()]))
                             for label,choice,factor in choices]
        # Independently execute the chosen finite actions, rather than trust
        # a reported SAT status without checking the transported assignment.
        source_action=list(range(nr))
        for action in actions:
            if action['stage']=='source_equivalence':
                source_action=[action['permutation'][i] for i in source_action]
        realized=dict(mapping)
        realized.update({i:model.eval(v).as_long() for i,v in domain_variables.items()})
        actual=[realized.get(i,np_) for i in range(nr)]
        for action in actions:
            if action['stage'].startswith('path_'):
                actual=[action['permutation'][i] if i<np_ else np_ for i in actual]
        normalized=[actual[i] for i in source_action]
        for action in actions:
            if action['stage']=='target_equivalence':
                normalized=[action['permutation'][i] if i<np_ else np_ for i in normalized]
        assert all(normalized[i]==p for i,p in reference.items())
        selected={i:p for i,p in enumerate(actual) if p<np_}
        if complete_reference:
            assert pynauty.certificate(colored_graph(features,selected))==pynauty.certificate(colored_graph(features,reference))
        return 'recovered', dict(actions=actions,heavy_mapping=sorted(selected.items()),
            singleton_assignments={i:realized[i] for i in domain_variables})
    return ('not_recovered',None) if status == z3.unsat else ('unknown',solver.reason_unknown())


def rank_key(mapping, problem):
    """Frozen graph-only representative rank; no reference data consumed."""
    mapping = dict(mapping)
    inverse = {p:r for r,p in mapping.items()}
    r,p = problem.reactant.wbo,problem.product.wbo
    broken=formed=changed=0
    for a,b in zip(*np.where(np.triu(r,1)>.2)):
        if a not in mapping or b not in mapping:
            broken += int(a in mapping or b in mapping)
        elif p[mapping[a],mapping[b]] <= .2: broken += 1
        elif abs(r[a,b]-p[mapping[a],mapping[b]])>.5: changed += 1
    for a,b in zip(*np.where(np.triu(p,1)>.2)):
        if a not in inverse or b not in inverse or r[inverse.get(a,0),inverse.get(b,0)]<=.2:
            formed += 1
    heavy = sum(problem.product.elements[i]!='H' for i in mapping.values())
    return (-heavy,-len(mapping),broken+formed+changed,tuple(sorted(mapping.items()))),dict(
        broken=int(broken),formed=int(formed),bond_order_changed=int(changed))


def evaluate(aam, features, reference, seconds=120, symbolic=True, query=None,
             query_timeout_ms=10000, reference_side='target', ranker=None):
    start=time.perf_counter(); deadline=start+seconds
    full_reference=dict(reference)
    reference=project(reference,features)
    side={'source':0,'target':1}[reference_side]
    complete_reference=len(reference)==len(features[side]['heavy'])
    ranker = ranker or (lambda mapping:rank_key(mapping,aam.problem))
    expected=pynauty.certificate(colored_graph(features,reference))
    terminals=list(aam.graph.terminals)
    ranked=sorted(terminals,key=lambda t:ranker(aam.graph.states[t].mapping)[0])
    candidates=[project(aam.graph.states[t].mapping,features) for t in ranked]
    # Explicit-H alternatives often have the identical heavy-atom relation.
    # Canonicalize that relation once, without discarding any AAM terminal.
    certificate_cache={}
    certificates=[]
    for mapping in candidates:
        key=tuple(sorted(mapping.items()))
        if key not in certificate_cache:
            certificate_cache[key]=pynauty.certificate(colored_graph(features,mapping))
        certificates.append(certificate_cache[key])
    hit=next((i for i,c in enumerate(certificates) if c==expected),None)
    result=dict(reference_pairs=len(reference),reference_annotation_complete=complete_reference,
        top1_correct=bool(certificates and certificates[0]==expected) if complete_reference else None,
        representative_recovery=hit is not None,reference_recovery='recovered' if hit is not None else 'not_recovered',
        top_terminal=ranked[0] if ranked else None,
        witness_terminal=ranked[hit] if hit is not None else None,
        candidate_terminals=len(terminals),unique_representative_chemistries=len(set(certificates)),
        unique_heavy_representatives=len(certificate_cache),
        best_target_heavy_coverage=max((len(m)/max(1,len(features[side]['heavy'])) for m in candidates),default=0),
        best_target_all_atom_coverage=max((len(aam.graph.states[t].mapping)/(aam.problem.source_atom_count if side==0 else aam.problem.target_atom_count) for t in terminals),default=0),
        capped=aam.graph.capped, symbolic_queries=0, unknown_queries=0)
    if hit is None and symbolic:
        seen=set(); rejected_projections=set(); target_index={a:i for i,a in enumerate(features[1]['heavy'])}
        source_heavy=set(features[0]['heavy'])
        def lifted(feature,degree):
            result=[]
            for g in endpoint_generators(feature):
                images=list(range(degree))
                for i,a in enumerate(feature['heavy']):images[a]=feature['heavy'][g[i]]
                result.append(tuple(images))
            return tuple(result)
        source_equivalence=lifted(features[0],aam.problem.source_atom_count)
        target_equivalence=lifted(features[1],aam.problem.target_atom_count)
        for path in aam.graph.paths():
            if time.perf_counter()>deadline:
                result['reference_recovery']='unknown';break
            changes_chemistry=False;projected_history=[]
            for edge in path.transitions:
                placement=aam.graph.fragment_placement(edge)
                if placement is None:continue
                projected=tuple(dict.fromkeys(tuple(target_index[g.images[a]] for a in features[1]['heavy'])
                                               for g in placement.target_generators))
                projected=tuple(g for g in projected if g!=tuple(range(len(target_index))))
                changes_chemistry |= any(not exact_action(g,features[1]) for g in projected)
                projected_domains=[]
                for domain in placement.symmetry_domains:
                    if domain.source=='exact_automorph_group':continue
                    atoms=[target_index[a] for a in domain.p_atoms if a in target_index]
                    if atoms:projected_domains.append((tuple(r for r in domain.r_atoms if r in source_heavy),tuple(atoms)))
                    for atom in atoms[1:]:
                        swap=list(range(len(target_index)));swap[atoms[0]],swap[atom]=swap[atom],swap[atoms[0]]
                        changes_chemistry |= not exact_action(swap,features[1])
                fragment=tuple(r for r in placement.r_atoms if r in source_heavy)
                if projected or projected_domains or fragment:
                    projected_history.append((fragment,tuple(projected_domains),projected,
                        tuple((r,p) for r,p in placement.representative_assignments if r in source_heavy),
                        tuple(r for r in placement.exact_fixed if r in source_heavy),
                        tuple(e for e in placement.deferred_edges if set(e)<=source_heavy),
                        tuple((r,p) for r,p in aam.graph.states[aam.graph.transitions[edge].source].mapping if r in source_heavy)))
            heavy_mapping=project(path.mapping,features)
            # Full H and domain history matter to feasibility, even when the
            # projected heavy representative is identical.
            key=(tuple(sorted(path.mapping.items())),path.context.cuts,frozen_value(path.fragments))
            if key in seen:continue
            seen.add(key)
            if complete_reference and (not changes_chemistry or len(heavy_mapping)!=len(reference)):continue
            query_function=query_path if query is None else query
            projection_key=(tuple(sorted(heavy_mapping.items())),tuple(projected_history),
                tuple(e for e in path.context.cuts if set(e)<=source_heavy))
            if projection_key in rejected_projections:continue
            if query is None:
                status,_=query_path(path,aam.problem,full_reference,
                    source_atoms=features[0]['heavy'],source_generators=source_equivalence,
                    target_generators=target_equivalence,complete_reference=complete_reference,
                    projected_atoms=features[0]['heavy'],
                    timeout_ms=min(query_timeout_ms,1000*(deadline-time.perf_counter())))
                if status=='not_recovered':rejected_projections.add(projection_key);continue
            status,witness=query_function(path,aam.problem,full_reference,
                source_atoms=features[0]['heavy'],source_generators=source_equivalence,
                target_generators=target_equivalence,complete_reference=complete_reference,
                timeout_ms=min(query_timeout_ms,1000*(deadline-time.perf_counter())))
            result['symbolic_queries']+=1
            if status=='recovered':
                realized=project(witness['mapping'],features)
                if complete_reference:
                    assert pynauty.certificate(colored_graph(features,realized))==expected
                witness['heavy_mapping']=sorted(realized.items())
                result.update(reference_recovery=status,witness_terminal=path.terminal,
                              witness_path=list(path.transitions),
                              witness_actions=witness);break
            if status=='unknown':
                result['unknown_queries']+=1;result['reference_recovery']='unknown'
    result['evaluation_seconds']=time.perf_counter()-start
    result['verifier']='full_explicit_domain_and_group_v1'
    if hit is None and symbolic:result['rejected_heavy_projections']=len(rejected_projections)
    if ranked:result['top_events']=ranker(aam.graph.states[ranked[0]].mapping)[1]
    return result


def evaluate_planned(aam, plan, features, reference, **options):
    """Score native compressed actions in search space; rank in input R/P space."""
    ordered_features = list(reversed(features)) if plan.reversed else features
    result = evaluate(aam,ordered_features,plan.to_search_mapping(reference),
        reference_side='source' if plan.reversed else 'target',
        ranker=lambda m:rank_key(plan.to_input_mapping(m),plan.input_problem),**options)
    result['search_direction']=plan.direction
    result['mapping_space']='search source/target; use orientation.json to recover input R/P indices'
    result['coverage_space']='original input product'
    if result['witness_terminal'] is not None and result['representative_recovery']:
        result['input_orientation_witness']=sorted(plan.to_input_mapping(
            aam.graph.states[result['witness_terminal']].mapping).items())
    elif result.get('witness_actions'):
        heavy=dict(result['witness_actions']['heavy_mapping'])
        concrete={ordered_features[0]['heavy'][r]:ordered_features[1]['heavy'][p]
                  for r,p in heavy.items()}
        result['input_orientation_heavy_witness']=sorted(plan.to_input_mapping(concrete).items())
        result['input_orientation_witness']=sorted(plan.to_input_mapping(result['witness_actions']['mapping']).items())
    return result
