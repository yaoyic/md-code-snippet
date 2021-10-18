import mdtraj as md

aa_map = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR"
}

def get_top_ca(seq, zero_based=True):
    if isinstance(seq, str):
        seq = [aa_map[aa] for aa in seq]
    top = md.Topology()
    chain = top.add_chain()
    elem_CA = md.element.carbon
    if zero_based:
        start = 0
    else:
        start = 1
    new_atom = None
    for i, aa in enumerate(seq):
        resi = top.add_residue(aa, chain, resSeq=start + i)
        old_atom, new_atom = new_atom, top.add_atom("CA", elem_CA, resi)
        if i > 0:
            top.add_bond(old_atom, new_atom)
    return top

if __name__ == "__main__":
    # example
    top = get_top_ca(["ALA", "GLY", "LYS", "THR"], zero_based=True)
    print(top)

    # alternatively
    top2 = get_top_ca("AGKT", zero_based=True)
    print(top)

