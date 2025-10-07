# Author: ChatGPT 5
# Assistant: Cameron F. Abrams <cfa22@drexel.edu>

import networkx as nx
import logging

logger = logging.getLogger(__name__)

def elem(G: nx.Graph, n: int, element_attr: str = 'element'):
    e = G.nodes[n].get(element_attr)
    return e if e else str(n)[0]  # fallback: first char of label

def is_H(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'H'
def is_C(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'C'
def is_O(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'O'
def is_N(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'N'
def is_P(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'P'
def is_S(G: nx.Graph, n: int, element_attr: str = 'element'):  return elem(G, n, element_attr) == 'S'

def mark_cc_doubles_by_degree(G: nx.Graph, element_attr: str = 'element', set_other_cc_single: bool = True):
    """
    Set bond 'order' on C-C edges using node degrees only:
      - C-C with both endpoints degree==3 -> order=2
      - (optional) all other C-C -> order=1
    Non-C-C edges are untouched.
    Returns the list of edges marked as double.
    """
    deg = dict(G.degree())  # snapshot the DegreeView once

    # Optionally preset all C-C to single
    if set_other_cc_single:
        for u, v in G.edges():
            if is_C(G, u, element_attr) and is_C(G, v, element_attr):
                G.edges[u, v]['order'] = 1

    cc_doubles = []
    for u, v in G.edges():
        if is_C(G, u, element_attr) and is_C(G, v, element_attr):
            if deg.get(u, 0) == 3 and deg.get(v, 0) == 3:
                G.edges[u, v]['order'] = 2
                cc_doubles.append((u, v))
    return cc_doubles

import networkx as nx

HETERO_HEAVY = {'O', 'N', 'P', 'S'}  # extend if you like (F, Cl, Br, I, etc.)
ATOMIC_Z = {'H':1, 'C':6, 'N':7, 'O':8, 'P':15, 'S':16}  # extend if needed

# ---------------- helpers ----------------
def heavy_neighbors(G: nx.Graph, n: int, element_attr: str = 'element'):
    return [x for x in G.neighbors(n) if not is_H(G, x, element_attr)]

def heavy_degree(G: nx.Graph, n: int, element_attr: str = 'element'):
    return sum(1 for _ in heavy_neighbors(G, n, element_attr))

def heavy_only_nodes(G: nx.Graph, element_attr: str = 'element'):
    return [n for n in G.nodes if not is_H(G, n, element_attr)]

def cycle_atoms_heavy(G: nx.Graph, element_attr: str = 'element'):
    """Return nodes on any cycle, computed on heavy-atom induced subgraph."""
    Hsub = G.subgraph(heavy_only_nodes(G, element_attr))
    cyc_nodes = set()
    for comp in nx.connected_components(Hsub):
        sub = Hsub.subgraph(comp)
        for cyc in nx.cycle_basis(sub):
            cyc_nodes.update(cyc)
    return cyc_nodes

def seed_set(G: nx.Graph, branch_deg_threshold: int = 3, element_attr: str = 'element'):
    """
    Seeds = hetero heavy atoms, ring atoms (heavy graph), and heavy-branch carbons (heavy_degree>=3).
    Hydrogens are never seeds.
    """
    cyc = cycle_atoms_heavy(G, element_attr)
    S = set()
    for n in G.nodes:
        if is_H(G, n, element_attr):
            continue  # exclude hydrogens from seeds
        e = elem(G, n, element_attr)
        if e in HETERO_HEAVY:
            S.add(n)
            continue
        if n in cyc:
            S.add(n)
            continue
        if e == 'C' and heavy_degree(G, n, element_attr) >= branch_deg_threshold:
            S.add(n)
    return S

def carbon_neighbors_outside_S(G: nx.Graph, node: int, S: set[int], element_attr: str = 'element'):
    """Immediate carbon neighbors not in S."""
    return [v for v in G.neighbors(node)
            if v not in S and is_C(G, v, element_attr)]

def carbon_component_from_neighbor(G: nx.Graph, start: int, S: set[int], element_attr: str = 'element'):
    """
    Carbon-only connected component reachable from 'start' without entering S.
    Hydrogens are ignored implicitly since we require is_C.
    """
    frontier = [start]
    visited = set()
    while frontier:
        v = frontier.pop()
        if v in visited:
            continue
        visited.add(v)
        for w in G.neighbors(v):
            if w in S:
                continue
            if is_C(G, w, element_attr) and w not in visited:
                frontier.append(w)
    return G.subgraph(visited).copy()

def leafs_in(sub: nx.Graph, avoid: int = None, element_attr: str = 'element'):
    """Degree-1 nodes within sub; optionally ignore 'avoid'."""
    leaves = []
    for n in sub.nodes:
        if avoid is not None and n == avoid:
            continue
        if sub.degree(n) <= 1:
            leaves.append(n)
    return leaves

def anchor_paths(sub: nx.Graph, anchor_adj: int, element_attr: str = 'element'):
    """One path per terminal leaf within the connected piece containing anchor_adj."""
    if anchor_adj not in sub:
        return []
    comp_nodes = nx.node_connected_component(sub, anchor_adj)
    comp = sub.subgraph(comp_nodes).copy()
    leaves = [n for n in comp.nodes if n != anchor_adj and comp.degree(n) <= 1]
    if not leaves:
        return [[anchor_adj]]
    return [nx.shortest_path(comp, anchor_adj, leaf) for leaf in leaves]

def count_unsats_cc(G: nx.Graph, path: list[int], element_attr: str = 'element') -> int:
    """Count C=C along a carbon path using 'order' if present; otherwise 0."""
    uns = 0
    for u, v in zip(path, path[1:]):
        if is_C(G, u, element_attr) and is_C(G, v, element_attr):
            if G.edges[u, v].get('order', 1) >= 2:
                uns += 1
    return uns

# ---------------- main API ----------------
def label_lipid_chains_heavy_aware(
    G: nx.Graph,
    min_len=6,
    split_branches=True,
    include_carbonyl_C_as_seed=True,
    element_attr: str = 'element',
):
    """
    Identify lipid tails robustly when hydrogens are present.

    Seeds:
      - hetero heavy atoms (O,N,P,S),
      - heavy-graph ring atoms,
      - heavy-branch carbons (heavy_degree >= 3),
      - (optionally) carbonyl carbons are naturally included via heavy_degree rule.

    Returns:
      node_to_chain_id: dict node -> int or None
      chains: list of dicts {id, anchor, nodes, path, length, unsaturations}
    """

    node_to_chain = {n: None for n in G.nodes}
    chains = []
    used = set()  # already assigned carbon nodes (when split_branches=True)

    # Build seeds on HEAVY logic (hydrogens excluded)
    S = seed_set(G, branch_deg_threshold=3, element_attr=element_attr)

    # Optionally ensure carbonyl C are treated as seeds even if not >=3 heavy neighbors
    # (Often they already are: neighbors are O(=), O(-), and C -> heavy_degree=3.)
    if include_carbonyl_C_as_seed:
        for n in heavy_only_nodes(G, element_attr):
            if not is_C(G, n, element_attr):
                continue
            # Carbonyl-like: has at least one oxygen neighbor; tend to be anchors
            if any(is_O(G, x, element_attr) for x in heavy_neighbors(G, n, element_attr)):
                S.add(n)

    # Case A: seeds exist (typical)
    if S:
        chain_id = 0
        for s in S:
            for cn in carbon_neighbors_outside_S(G, s, S, element_attr):
                if split_branches and cn in used:
                    continue
                comp = carbon_component_from_neighbor(G, cn, S, element_attr)
                if comp.number_of_nodes() == 0:
                    continue

                if split_branches:
                    paths = anchor_paths(comp, cn, element_attr)
                    for p in paths:
                        L = sum(1 for x in p if is_C(G, x, element_attr))
                        if L < min_len:
                            continue
                        uns = count_unsats_cc(G, p, element_attr)
                        for n in p:
                            node_to_chain[n] = chain_id
                            used.add(n)
                        chains.append({
                            'id': chain_id,
                            'anchor': s,
                            'nodes': set(p),
                            'path': p,
                            'length': L,
                            'unsaturations': uns,
                        })
                        chain_id += 1
                else:
                    comp_nodes = set(comp.nodes)
                    L = sum(1 for n in comp_nodes if is_C(G, n, element_attr))
                    if L >= min_len:
                        leaves = [n for n in comp.nodes if comp.degree(n) <= 1]
                        if leaves:
                            a = leaves[0]
                            b = max(nx.single_source_shortest_path_length(comp, a),
                                    key=lambda kv: kv[1])[0]
                            c = max(nx.single_source_shortest_path_length(comp, b),
                                    key=lambda kv: kv[1])[0]
                            p = nx.shortest_path(comp, b, c)
                        else:
                            p = [cn]
                        uns = count_unsats_cc(G, p, element_attr)
                        for n in comp_nodes:
                            node_to_chain[n] = chain_id
                            used.add(n)
                        chains.append({
                            'id': chain_id,
                            'anchor': s,
                            'nodes': comp_nodes,
                            'path': p,
                            'length': L,
                            'unsaturations': uns,
                        })
                        chain_id += 1

    # Case B: no seeds (pure hydrocarbon, no rings/branches/heteroatoms)
    if not chains and not S and G.number_of_nodes() > 0:
        comp_nodes = max(nx.connected_components(G), key=len)
        sub = G.subgraph(comp_nodes)
        a = next(iter(sub.nodes))
        b = max(nx.single_source_shortest_path_length(sub, a), key=lambda kv: kv[1])[0]
        c = max(nx.single_source_shortest_path_length(sub, b), key=lambda kv: kv[1])[0]
        p = nx.shortest_path(sub, b, c)
        L = sum(1 for x in p if is_C(G, x, element_attr))
        uns = count_unsats_cc(G, p, element_attr)
        for n in p:
            node_to_chain[n] = 0
        chains.append({
            'id': 0,
            'anchor': None,
            'nodes': set(p),
            'path': p,
            'length': L,
            'unsaturations': uns,
        })

    for h in G.nodes:
        if not is_H(G, h, element_attr):
            continue
        # only consider C–H; ignore O–H, N–H, etc.
        carbon_nbrs = [n for n in G.neighbors(h) if is_C(G, n, element_attr)]
        if len(carbon_nbrs) != 1:
            continue  # skip weird/ambiguous cases
        c = carbon_nbrs[0]
        cid = node_to_chain.get(c)
        if cid is not None:
            node_to_chain[h] = cid
    return node_to_chain, chains

def _is_heavy(G, n, element_attr: str = "element"):  # exclude hydrogens
    return not is_H(G, n, element_attr)

def _is_hetero_heavy(G, n, element_attr: str = "element"):
    return not is_H(G, n, element_attr) and not is_C(G, n, element_attr)

def _Z(G, n, element_attr: str = "element"):
    return ATOMIC_Z.get(elem(G, n, element_attr), 0)

def head_tail_pairs(G: nx.Graph, chains: list[dict], element_attr: str = "element"):
    """
    For each chain in `chains`, return:
      {'chain_id', 'head', 'tail', 'head_element', 'head_distance'}
    where:
      tail = last node of the chain's path
      head = heaviest hetero atom at the minimum heavy-graph distance from the chain's anchor
    """
    # Heavy-atom-only view for distances (ignore hydrogens)
    heavy_nodes = [n for n in G.nodes if _is_heavy(G, n, element_attr)]
    H = G.subgraph(heavy_nodes)
    hetero_nodes = [n for n in H.nodes if _is_hetero_heavy(G, n, element_attr)]

    results = []
    for ch in chains:
        cid   = ch['id']
        anchor = ch.get('anchor')
        path   = ch.get('path', [])
        tail   = path[-1] if path else None

        # Fallback if anchor is None: use path start (first carbon) as anchor proxy
        anchor_proxy = anchor if anchor is not None else (path[0] if path else None)

        head = None
        head_dist = None
        tail_dist = None

        if anchor_proxy in H:
            # Distances on heavy graph only
            dist = nx.single_source_shortest_path_length(H, anchor_proxy)
            # Candidates: hetero heavy reachable nodes
            candidates = [n for n in hetero_nodes if n in dist]
            if candidates:
                # First minimize distance, then maximize atomic number, tie-break by id
                dmin = min(dist[n] for n in candidates)
                near = [n for n in candidates if dist[n] == dmin]
                head = max(near, key=lambda n: (_Z(G, n), str(n)))
                head_dist = dmin
        tail_dist = nx.shortest_path_length(H, head, tail) if (anchor_proxy in H and tail in H) else None

        results.append({
            'chain_id'      : cid,
            'head'          : head,
            'tail'          : tail,
            'tail_distance' : tail_dist,
            'head_element'  : (elem(G, head, element_attr) if head is not None else None),
            'head_distance' : head_dist,
        })

    return results
