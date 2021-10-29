
import os, json, argparse, sys

"""
Load compounds from Wikidata that link to biosynthetic processes, save in canonical form (Q numbers sorted numerically), output statistics
"""
# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="perform SPARQL query",
        action="store_true")
parser.add_argument("-p", "--patternless", help="after stats, output list of compounds lacking any SMILES/SMARTS/InChI key",
        action="store_true")

# Read arguments from the command line
args = parser.parse_args()

# Check for --version or -V
dontquery = not args.query
script = os.path.basename(sys.argv[0])[:-3]
OFFSET = 1000000000

if dontquery is False:
    print('performing query ...', file=sys.stderr)
    ret = os.popen('wd sparql biosyn-compounds.rq >biosyn-compounds.json'.format(script, script))
    if ret.close() is not None:
        raise
    ff = open('biosyn-compounds.json'.format(script))
    s = ff.read()
    jol = json.loads(s)
    ff.close()
    items = set([d.get('item').get('value') for d in jol])
    query="""
    SELECT DISTINCT ?item ?p31
    WHERE
    {{
      VALUES ?item {{ {} }}
      ?item wdt:P31 ?p31.
    }}
    """.format("wd:" + " wd:".join(items))
    f = open('{}-1.rq'.format(script), 'w')
    f.write(query)
    f.close()

    print('querying P31 for {} items...'.format(len(items)))
    ret = os.popen('wd sparql {}-1.rq >{}.P31'.format(script, script))
    if ret.close() is not None:
        exit()

    itp31 = {}
    with open('{}.P31'.format(script), 'r') as ff:
        s = ff.read()
        jj = json.loads(s)
        for d in jj:
            it = d.get('item')
            p31 = d.get('p31')
            if p31 == 'Q11173':
                continue
            i = itp31.get(it)
            if i is None:
                itp31[it] = [p31]
            else:
                i.append(p31)

    for dd in jol:
        it = dd.get('item').get('value')
        p31 = itp31.get(it)
        if p31 is not None:
            dd['p31'] = p31

    with open('data-biosyn-compounds.json', 'w+') as f:
        f.write(json.dumps(sorted(jol, key=lambda data: OFFSET*int(data.get('item').get('value')[1:]) +\
            int(data.get('goid')[3:])), indent=0, ensure_ascii=False))
else:
    ff = open('data-biosyn-compounds.json'.format(script))
    s = ff.read()
    jol = json.loads(s)

items = set()
gos = {}
iks = {}
ik1 = {}
smarts = {}
smiles = {}
for d in jol:
    it = d.get('item').get('value')
    items.add(it)
    g = gos.get(it)
    gotup = (d.get('goid'), d.get('goLabel'))
    if g is not None:
        g.add(gotup)
        continue
    else:
        gos[it] = set(gotup)
    sma = d.get('p8533')
    if sma is not None:
        smarts[sma] = it
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1[ik[:14]] = it
    csm = d.get('p233')
    ism = d.get('p2017')
    smi = None
    if ism is not None:
        smi = ism
    elif csm is not None:
        smi = csm
    if smi is not None:
        smiles[smi] = it

pless = items.difference(set(iks.values()).union(set(smiles.values()).union(set(smarts.values()))))

print('#items: {}'.format(len(items)))
print('#InChI keys: {}'.format(len(iks.keys())))
print('#InChI1 keys: {}'.format(len(ik1.keys())))
print('#smarts: {}'.format(len(smarts.keys())))
print('#smiles: {}'.format(len(smiles.keys())))
print('# w/o pattern: {}'.format(len(pless)))

if args.patternless:
    for i in pless:
        print(gos.get(i))
