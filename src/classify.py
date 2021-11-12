import csv, os, json, argparse, sys, random
from rdkit import Chem

"""
Load data from directory, classify given molecule or test file 
"""

# Walk reverse ontology subgraph (root = hitem) upwards, depth first. If any member of hitemset
# is encountered, add it to minor set. These are later discarded, leaving only the most
# specialized hits
def walk_ont(sitems, hitemset, hitem, minors):
    ch = sitems.get(hitem)
    if ch is None:
        return
    for c in ch:
        if c in hitemset:
            minors.add(c)
        walk_ont(sitems, hitemset, c, minors)

# Walk reverse ontology subgraph (root = node) upwards, breadth first. The first encountered
# GO process is returned.
def walk_ont_go(sitems, node, gos):
    visited = set()    # List to keep track of visited nodes.
    queue = []      #Initialize a queue
    visited.add(node)
    queue.append(node)

    while queue:
        s = queue.pop(0)
        go = gos.get(s)
        if go is not None:
            return go
        ch = sitems.get(s)
        if ch is None:
            continue
        for c in ch:
            if c not in visited:
                visited.add(c)
                queue.append(c)
    return None

# Walk reverse ontology subgraph (root = hitem) upwards, depth first.
# Return True if citem is encountered (i.e. citem is superclass of hitem)
def walk_find(sitems, hitem, citem):
    ch = sitems.get(hitem)
    # print('{} {} {}'.format(hitem, citem, ch))
    if ch is None:
        return False
    if citem in ch:
        return True
    for c in ch:
        if walk_find(sitems, c, citem):
            return True
    return False

# Walk reverse ontology subgraph (root = hitem) upwards, depth first.
# If one of items in list is encountered return full path to it 
path = []
def walk_find_list(sitems, hitem, item_set):
    global path
    ch = sitems.get(hitem)
    # print('{} {} {}'.format(hitem, citem, ch))
    if ch is None:
        return False
    for c in ch:
        if c is None:
            continue
        if c in item_set:
            path.append(str(c))
            return True
        if walk_find_list(sitems, c, item_set):
            path.append(c)
            return True
    return False

# Find path to topmost natural product, return as array of names
def path_to_nproot(hit, nplist):
    global path
    if hit in nplist:
        return [hit]
    path = []
    success = walk_find_list(sitems, hit, set(nplist))
    if not success:
        return []
    path.append(hit)
    return path

ctr = 0
def is_elem(e, s):
    global ctr
    ctr = ctr + 1
    return e in s

# Walk ontology downwards, breadth first, starting from given root, add children to set
def walk_collect_npclasses(edges, node, collection):
    visited = set()    # List to keep track of visited nodes.
    queue = []      #Initialize a queue
    visited.add(node)
    queue.append(node)

    while queue:
        s = queue.pop(0)
        e = edges.get(s)
        if e is None:
            continue
        for c in e:
            if c not in visited:
                visited.add(c)
                queue.append(c)
    collection |= visited

#rule A-1
def check_unspec_alkaloid(mol):
    #ring-N
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[#7R]')):
        return False
    #!peptide
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([NR0]),$([Nr;!r3;!r4;!r5;!r6;!r7;!r8])]~C(~[OX1H0,OX2H1])')):
        return False
    #!tetrapyrroles
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]1~[#6]~[#6]2~[#6]~[#6]3~[#6]~[#6]~[#6](~[#6]~[#6]4~[#6]~[#6]~[#6](~[#6]~[#6]5~[#6]~[#6]~[#6](~[#6]~[#6]~1~[#7]~2)~[#7]~5)~[#7]~4)~[#7]~3')):
        return False
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]1~[#6]~[#7]~[#6](~[#6R0]~[#6]2~[#6]~[#6]~[#6](~[#6R0]~[#6]3~[#6]~[#6]~[#6](~[#6R0]~[#6]4~[#6]~[#6]~[#6]~[#7]~4)~[#7]~3)~[#7]~2)~[#6]~1')):
        return False
    #!nucleotides
    if mol.HasSubstructMatch(Chem.MolFromSmarts('POCC1CCC(n2cnc3cncnc32)O1')):
        return False
    if mol.HasSubstructMatch(Chem.MolFromSmarts('POCC2CCC([#7R1]1~[#6R1]~[#7R1]~[#6R1]~[#6R1]~[#6R1]~1)O2')):
        return False
    #!diketopiperazine
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]1~[#7]~[#6](~O)~[#6]~[#7]~[#6](~O)~1')):
        return False
    return True

# Try to match ALL patterns. Remove redundant. Return remaining.
def get_hits(mol, silent=False):
    if not silent:
        print('collecting hits...')
    hits = {}
    # match the InChI keys
    ik = Chem.MolToInchiKey(mol)
    it = iks.get(ik)
    if it is not None:
        go = ''
        g = gos.get(it)
        if g is not None:
            go = g
        hits[it] = (it, labels.get(it), ik)
        if args.rawhits:
            print('{} {} {} {}'.format(it, labels.get(it), ik, go))
    else:
        it = ik1s.get(ik[:14])
        if it is not None:
            go = ''
            g = gos.get(it)
            if g is not None:
                go = g
            hits[it] = (it, labels.get(it), ik[:14])
            if args.rawhits:
                print('{} {} {} {}'.format(it, labels.get(it), ik[:14], go))

    ps_default = Chem.AdjustQueryParameters()
    ps_default.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES | Chem.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
    ps_ignoreDummies = Chem.AdjustQueryParameters()
    ps_ignoreDummies.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
    for it,itsuplist in sitems.items():
        sm = smiles.get(it)
        smas = smarts.get(it)
        if sm is None and smas is None:
            continue
        if (not (sm is None or len(sm) == 0)) and sm.find('*') < 0:
            if sm.find('@') >= 0:
                continue
        # this class has a non-wildcard SMILES
            try:
                q = Chem.MolFromSmiles(sm)
                pats = [Chem.AdjustQueryProperties(q, ps_ignoreDummies)]
            except Exception:
                continue
        elif not smas is None:
        # this class has a wildcard SMILES or a SMARTS
            #print('---{} {} {}'.format(it, labels.get(it), sma))
            try:
                pats = [Chem.MolFromSmarts(sma) for sma in smas]
                #pat = Chem.AdjustQueryProperties(pat, ps_default)
                #if it == 'Q12748271':
                #    print(mol.HasSubstructMatch(pat))
                #    exit()
                #print('====={} {}'.format(it, sma))
            except Exception:
                continue
        else:
            try:
                q = Chem.MolFromSmiles(sm)
                pats = [Chem.AdjustQueryProperties(q, ps_ignoreDummies)]
            except Exception:
                continue
        
        if pats is None:
            print('========defect pattern: {} {}'.format(it, sma))
            continue
        for p in pats:
            if mol.HasSubstructMatch(p):
                go = ''
                g = gos.get(it)
                if g is not None:
                    go = g
                hits[it] = (it, labels.get(it), sm)
                if args.rawhits:
                    print('{} {} {} {}'.format(it, labels.get(it), sm, go))
                break

    if check_unspec_alkaloid(mol):
        hits['Q70702'] = ('Q70702', 'unspecified alkaloid', 'Rule A-1')
    
    if not silent:
        print('purging redundant hits')
    hitemset = set(hits.keys())
    minors = set()
    for hit in hitemset:
        walk_ont(sitems, hitemset, hit, minors)
    for m in minors:
        hits.pop(m)

    return hits


# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", help="data directory", required=True)
parser.add_argument("-m", "--molecule", help="molecule to classify (SMILES/InChi)")
parser.add_argument("-a", "--rawhits", help="output all hits before hit processing starts",
        action="store_true")
parser.add_argument("-t", "--testfile", help="classify all InChis from file, first line may contain check item", type=str)
parser.add_argument("-n", "--nptest", help="Load N random NP Inchis from WD and classify", action="store_true")
parser.add_argument("-j", "--json", help="together with -m outputs JSON formatted result", action="store_true")
parser.add_argument("-N", "--natural", help="prune ontology to only include natural products", action="store_true")
parser.add_argument("--list_nproots", help="list of topmost natural products classes", action="store_true")

# Read arguments from the command line
args = parser.parse_args()
is_normal_run = not (args.molecule is None) or not(args.testfile is None) or args.nptest
is_extra_service = args.list_nproots

if not is_normal_run and not is_extra_service:
    print('One of -m or -t or -n is needed')
    parser.print_usage()
    exit()

silent = args.json

if not silent:
    print(args)

if is_normal_run and args.testfile is None and args.nptest is False:
    if args.molecule.find('InChI=') >= 0:
        mol = Chem.MolFromInchi(args.molecule)
    else:
        mol = Chem.MolFromSmiles(args.molecule)


iks = {}
ik1s = {}
labels = {}
smiles = {}
smarts = {}

# read bio process data
with open(args.data + 'data-biosyn-classes.json', 'r') as f:
    if not silent:
        print('reading biosyn class data')
    s = f.read()
    jol = json.loads(s)

gos = {}
for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1s[ik[:14]] = it
    else:
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
            t = smarts.get(it)
            if t is None:
                smarts[it] = [p8533]
            else:
                t.append(p8533)
    g = gos.get(it)
    gotup = (d.get('goid'), d.get('goLabel'))
    if g is not None:
        g.add(gotup)
        continue
    else:
        gos[it] = set(gotup)
    lab = dd.get('label')
    labels[it] = lab

sitems = {}
with open(args.data + 'data-biosyn-compounds.json', 'r') as f:
    if not silent:
        print('reading biosyn compound data')
    s = f.read()
    jol = json.loads(s)

for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1s[ik[:14]] = it
    else:
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
            t = smarts.get(it)
            if t is None:
                smarts[it] = [p8533]
            else:
                t.append(p8533)
    g = gos.get(it)
    gotup = (d.get('goid'), d.get('goLabel'))
    if g is not None:
        g.add(gotup)
        continue
    else:
        gos[it] = set(gotup)
    lab = dd.get('label')
    labels[it] = lab
    p31 = d.get('p31')
    if p31 is not None:
        sitems[it] = p31

# read pattern data
with open(args.data + 'data-class-pattern.json', 'r') as f:
    if not silent:
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
    else:
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
            t = smarts.get(it)
            if t is None:
                smarts[it] = [p8533]
            else:
                t.append(p8533)
    labels[it] = lab

childs = {}
# filling subclass structure
with open(args.data + 'data-class-subclass.json', 'r') as f:
    if not silent:
        print('reading superclass data')
    s = f.read()
    jol = json.loads(s)

for d in jol:
    it = d.get('item')
    sup = d.get('super')
    if not args.natural:
        i = sitems.get(it)
        if i is not None:
            i.append(sup)
        else:
            sitems[it] = [sup]
    c = childs.get(sup)
    if c is not None:
        c.append(it)
    else:
        childs[sup] = [it]

with open(args.data + 'data-class-subclass-names.json', 'r') as ff:
    s = ff.read()
    jol = json.loads(s)

for d in jol:
    it = d.get('item')
    lab = d.get('lab')
    if lab is None:
        continue
    labels[it] = lab

# root items missing from query data
sitems['Q2393187'] = ['Q43460564']
labels['Q2393187'] = 'molecular entity'
labels['Q43460564'] = 'chemical entity'

#load list of NP roots
nplist = set()
with open('natural.txt', 'r') as nf:
    nl = nf.readlines()
    nplist = set([line.rstrip() for line in nl])

if args.list_nproots:
    j = []
    for i in nplist:
        d = {}
        d['item'] = i
        d['name'] = labels.get(i)
        j.append(d)
    print(json.dumps(j))
    exit()

#create specialized ontology
if args.natural:
    nps = set()
    for nproot in nplist:
        walk_collect_npclasses(childs, nproot, nps)
#        print((nproot,len(nps),ctr))
        ctr = 0
    for np in nps:
        clist = childs.get(np)
        if clist is None:
            continue
        for c in clist:
            i = sitems.get(c)
            if i is not None:
                i.append(np)
            else:
                sitems[c] = [np]

if args.testfile is not None:
    print('classifying testfile {}'.format(args.testfile))
    with open(args.testfile, 'r') as tf:
        tfs = [line.strip() for line in tf.readlines()]
    count = 0
    fcount = 0
    for t_inchi in tfs:
        if t_inchi.startswith('Q'):
            check_item = t_inchi
            print('test class: {}'.format(check_item))
            continue
        mol = Chem.MolFromInchi(t_inchi)
        if mol is None:
            print('{} {}'.format(count, t_inchi))
            exit()
        hits = get_hits(mol, silent=True)
        count = count + 1
        if hits.get(check_item) is None and all([walk_find(sitems, hit, check_item) == False for hit in hits.keys()]):
            print('{} FAIL: {}'.format(count, t_inchi))
            fcount = fcount + 1
        else:
            print('{} OK'.format(count, hits))
    print('FAIL: {}/{}'.format(fcount, count))
elif args.nptest:
    N = 100
    print('classifying {} random small biomolecules'.format(N))
    with open(args.data + 'test/npitems.txt', 'r') as npf:
        qs = [line.strip() for line in npf.readlines()]
    sample = random.choices(qs, k=N)
    query="""
    SELECT ?inchi
    WHERE
    {{
      VALUES ?item {{ {} }}
      ?item wdt:P234 ?inchi.
    }}
    """.format('wd:' + ' wd:'.join(sample))
    f = open('npsample.rq', 'w')
    f.write(query)
    f.close()
    ret = os.popen('wd sparql npsample.rq >npsample.json')
    if ret.close() is not None:
        print('query failed')
        exit()
    with open('npsample.json', 'r') as tf:
        tfs = [line.strip() for line in tf.readlines()]
    count = 0
    fcount = 0
    check_item = 'Q206229'
    for t_inchi in tfs:
        mol = Chem.MolFromInchi(t_inchi)
        if mol is None:
            print('{} {}'.format(count, t_inchi))
            exit()
        hits = get_hits(mol, silent=True)
        count = count + 1
        found = False
        for hit in hits.keys():
            if walk_find(sitems, hit, check_item):
                print('{} OK {}'.format(count, hits.get(hit)))
                found = True
                break
        if not found:
            print('{} FAIL: {}'.format(count, t_inchi))
            fcount = fcount + 1
    print('FAIL: {}/{}'.format(fcount, count))
else:
    silent = False
    if args.json:
        silent = True
    hits = get_hits(mol, silent)
    if args.rawhits:
        print('----------------------------------')
    j = []
    for hit in hits.keys():
        if args.json:
            p = path_to_nproot(hit, nplist)
            if len(p) == 0:
                continue
            d = {}
            d['classification_names'] = [labels.get(it) for it in p]
        #finding closest biosynthetic processes
        cgo = walk_ont_go(sitems, hit, gos)
        if args.json:
            if cgo is not None:
                d['biological_process'] = list(cgo)
            j.append(d)
        else:
            print('{} {}'.format(hits.get(hit), cgo))
    if args.json:
        print(json.dumps(j))


