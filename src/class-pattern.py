
import os, json, argparse, sys

"""
Load chemclasses from Wikidata together with any pattern, save in canonical form (Q numbers sorted numerically), output statistics
"""
# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="perform SPARQL query",
        action="store_true")

# Read arguments from the command line
args = parser.parse_args()

# Check for --version or -V
dontquery = not args.query
script = os.path.basename(sys.argv[0])[:-3]
OFFSET = 1000000000

if dontquery is False:
    print('performing query...', file=sys.stderr)
    ret = os.popen('wd sparql class-pattern.rq >class-pattern.json'.format(script, script))
    if ret.close() is not None:
        raise
    ff = open('class-pattern.json'.format(script))
    s = ff.read()
    jol = json.loads(s)
    with open('data-class-pattern.json', 'w+') as f:
        f.write(json.dumps(sorted(jol, key=lambda data: int(data.get('item').get('value')[1:])), indent=0, ensure_ascii=False))
else:
    ff = open('data-class-pattern.json'.format(script))
    s = ff.read()
    jol = json.loads(s)

items = set()
iks = {}
ik1 = {}
smarts = {}
smiles = {}
for d in jol:
    it = d.get('item').get('value')
    items.add(it)
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

print('#items: {}'.format(len(items)))
print('#InChI keys: {}'.format(len(iks.keys())))
print('#InChI1 keys: {}'.format(len(ik1.keys())))
print('#smarts: {}'.format(len(smarts.keys())))
print('#smiles: {}'.format(len(smiles.keys())))
print('#wildcard smiles: {}'.format(sum(sm.find('*') >= 0 for sm in smiles.keys())))


