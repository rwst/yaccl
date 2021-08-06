import csv, os, json, argparse, sys, pronto
from rdkit import Chem

"""
Load classes from file, classify given molecule
"""
# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", help="data file",
        required=True)
parser.add_argument("-m", "--molecule", help="molecule to classify (SMILES/InChi)",
        required=True)
parser.add_argument("-g", "--go", help="Gene Ontology file",
        required=True)
parser.add_argument("-n", "--natural", help="file with root items of natural products",
        required=True)

# Read arguments from the command line
args = parser.parse_args()
if args.molecule.find('InChI=') >= 0:
    mol = Chem.MolFromInchi(args.molecule)
else:
    mol = Chem.MolFromSmiles(args.molecule)

# Check for --version or -V
script = os.path.basename(sys.argv[0])[:-3]

ff = open(args.data)
s = ff.read()
jol = json.loads(s)

print('reading data')
sitems = {}
labels = {}
smiles = {}
smarts = {}
for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    lab = dd.get('label')
    sup = d.get('super')
    go = d.get('goid')
    sm = None
    p233 = None
    p2017 = None
    p233 = d.get('p233')
    p2017 = d.get('p2017')
    if p2017 is not None and len(p2017) > 0:
        sm = p2017
    elif p233 is not None and len(p233) > 0:
        sm = p233
    if sm is not None:
        smiles[it] = sm
    p8533 = d.get('p8533')
    if p8533 is not None and len(p8533) > 0:
        smarts[it] = p8533

    i = sitems.get(it)
    if i is not None:
        i.append(sup)
    else:
        sitems[it] = [sup]
    labels[it] = lab
sitems['Q2393187'] = 'Q43460564'
labels['Q2393187'] = 'molecular entity'
labels['Q43460564'] = 'chemical entity'

# reverse links
edges = {}
for it,itsuplist in sitems.items():
    for sup in itsuplist:
        if sup in set(sitems.keys()).union(set(['Q43460564'])):
            e = edges.get(sup)
            if e is None:
                edges[sup] = set([it])
            else:
                e.add(it)

nplist = set()
with open(args.natural, 'r') as nf:
    nl = nf.readlines()
    nplist = set([line.rstrip() for line in nl])

ps = Chem.AdjustQueryParameters()
ps.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
for it,itsuplist in sitems.items():
    sm = smiles.get(it)
    sma = smarts.get(it)
    if sm is None and sma is None:
        continue
    if not sm is None and sm.find('*') < 0:
        if sm.find('@') >= 0:
            continue
    # this class has a non-wildcard SMILES
        try:
            q = Chem.MolFromSmiles(sm)
            pat = Chem.AdjustQueryProperties(q, ps)
        except Exception:
            continue
    elif not sma is None:
    # this class has a wildcard SMILES or a SMARTS
        #print('---{} {} {}'.format(it, labels.get(it), sma))
        pat = Chem.MolFromSmarts(sma)
        sm = sma # for reporting
    else:
        try:
            q = Chem.MolFromSmiles(sm)
            pat = Chem.AdjustQueryProperties(q, ps)
        except Exception:
            continue
    
    if mol.HasSubstructMatch(pat):
        print('{} {} {}'.format(it, labels.get(it), sm))
        
#print('reading GO')
