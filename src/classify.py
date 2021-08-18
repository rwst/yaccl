import csv, os, json, argparse, sys, pronto
from rdkit import Chem

"""
Load classes from directory, classify given molecule
"""
def walk_ont(sitems, hitemset, hitem, minors):
    ch = sitems.get(hitem)
    if ch is None:
        return
    for c in ch:
        if c in hitemset:
            minors.add(c)
        walk_ont(sitems, hitemset, c, minors)


# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", help="data directory",
        required=True)
parser.add_argument("-m", "--molecule", help="molecule to classify (SMILES/InChi)",
        required=True)

# Read arguments from the command line
args = parser.parse_args()
if args.molecule.find('InChI=') >= 0:
    mol = Chem.MolFromInchi(args.molecule)
else:
    mol = Chem.MolFromSmiles(args.molecule)

# Check for --version or -V
script = os.path.basename(sys.argv[0])[:-3]

# read bio process data
with open(args.data + 'data-biosyn.json', 'r') as f:
    print('reading biosyn data')
    s = f.read()
    jol = json.loads(s)

iks = {}
ik1s = {}
labels = {}
smiles = {}
smarts = {}
gos = {}
for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1s[ik[:14]] = it
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
    g = gos.get(it)
    gotup = (d.get('goid'), d.get('goLabel'))
    if g is not None:
        g.add(gotup)
        continue
    else:
        gos[it] = set(gotup)
    lab = dd.get('label')
    labels[it] = lab

# read pattern data
with open(args.data + 'data-class-pattern.json', 'r') as f:
    print('reading class pattern data')
    s = f.read()
    jol = json.loads(s)

for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    g = gos.get(it)
    if g is not None:
        continue
    lab = dd.get('label')
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1s[ik[:14]] = it
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
    labels[it] = lab

sitems = {}
with open(args.data + 'data-class-subclass.json', 'r') as f:
    print('reading superclass data')
    s = f.read()
    jol = json.loads(s)

for d in jol:
    it = d.get('item')
    sup = d.get('super')
    i = sitems.get(it)
    if i is not None:
        i.append(sup)
    else:
        sitems[it] = [sup]

sitems['Q2393187'] = 'Q43460564'
labels['Q2393187'] = 'molecular entity'
labels['Q43460564'] = 'chemical entity'
"""
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
"""

# match the InChI keys
ik = Chem.MolToInchiKey(mol)
it = iks.get(ik)
if it is not None:
    go = ''
    g = gos.get(it)
    if g is not None:
        go = g
    print('{} {} {} {}'.format(it, labels.get(it), ik, go))
else:
    it = ik1s.get(ik[:14])
    if it is not None:
        go = ''
        g = gos.get(it)
        if g is not None:
            go = g
        print('{} {} {} {}'.format(it, labels.get(it), ik[:14], go))

ps_default = Chem.AdjustQueryParameters()
ps_default.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES | Chem.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
ps_ignoreDummies = Chem.AdjustQueryParameters()
ps_ignoreDummies.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
print('collecting hits...')
hits = {}
for it,itsuplist in sitems.items():
    sm = smiles.get(it)
    sma = smarts.get(it)
    if sm is None and sma is None:
        continue
    if (not (sm is None or len(sm) == 0)) and sm.find('*') < 0:
        if sm.find('@') >= 0:
            continue
    # this class has a non-wildcard SMILES
        try:
            q = Chem.MolFromSmiles(sm)
            pat = Chem.AdjustQueryProperties(q, ps_ignoreDummies)
        except Exception:
            continue
    elif not sma is None:
    # this class has a wildcard SMILES or a SMARTS
        #print('---{} {} {}'.format(it, labels.get(it), sma))
        try:
            pat = Chem.MolFromSmarts(sma)
            #pat = Chem.AdjustQueryProperties(pat, ps_default)
            #if it == 'Q12748271':
            #    print(mol.HasSubstructMatch(pat))
            #    exit()
            #print('====={} {}'.format(it, sma))
            sm = sma # for reporting
        except Exception:
            continue
    else:
        try:
            q = Chem.MolFromSmiles(sm)
            pat = Chem.AdjustQueryProperties(q, ps_ignoreDummies)
        except Exception:
            continue
    
    if pat is None:
        print('========defect pattern: {} {}'.format(it, sma))
        continue
    if mol.HasSubstructMatch(pat):
        go = ''
        g = gos.get(it)
        if g is not None:
            go = g
        hits[it] = (it, labels.get(it), sm, go)
        print('{} {} {} {}'.format(it, labels.get(it), sm, go))

print('purging redundant hits')
hitemset = set(hits.keys())
minors = set()
for hit in hitemset:
    walk_ont(sitems, hitemset, hit, minors)
for m in minors:
    hits.pop(m)

print('------remaining hits:')
for hit in hits.keys():
    print(hits.get(hit))
