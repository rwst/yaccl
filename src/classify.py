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
#parser.add_argument("-n", "--natural", help="file with root items of natural products",
#        required=True)
parser.add_argument("-s", "--super", help="superclass database file default data-class-subclass.json",
        default='data-class-subclass.json')

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
iks = {}
ik1s = {}
labels = {}
smiles = {}
smarts = {}
gos = {}
for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    lab = dd.get('label')
    go = d.get('goid')
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
    labels[it] = lab

sitems = {}
with open(args.super, 'r') as f:
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
"""
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
        print('{} {} {} {}'.format(it, labels.get(it), sm, go))
        
