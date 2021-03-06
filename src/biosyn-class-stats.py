
import os, json, argparse, sys

"""
Load chemclasses from Wikidata that link to biosynthetic processes, save in canonical form (Q numbers sorted numerically), output statistics
"""
# see https://stackoverflow.com/questions/793761/built-in-python-hash-function
def c_mul(a, b):
  return eval(hex((int(a) * b) & (2**64 - 1))[:-1])

def py25hash(self):
  if not self:
    return 0 # empty
  value = ord(self[0]) << 7
  for char in self:
    value = c_mul(1000003, value) ^ ord(char)
  value = value ^ len(self)
  if value == -1:
    value = -2
  if value >= 2**63:
    value -= 2**64
  return value

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="perform SPARQL query",
        action="store_true")
parser.add_argument("-p", "--patternless", help="after stats, output list of processes where conected chemical items lack any SMILES/SMARTS/InChI key",
        action="store_true")

# Read arguments from the command line
args = parser.parse_args()

# Check for --version or -V
dontquery = not args.query
script = os.path.basename(sys.argv[0])[:-3]
OFFSET1 = 1000000000*2**64

if dontquery is False:
    print('performing query...', file=sys.stderr)
    ret = os.popen('wd sparql biosyn-classes.rq >biosyn-classes.json'.format(script, script))
    if ret.close() is not None:
        raise
    ff = open('biosyn-classes.json'.format(script))
    s = ff.read()
    jol = json.loads(s)

if dontquery:
    ff = open('data-biosyn-classes.json'.format(script))
    s = ff.read()
    jol = json.loads(s)
    ff.close()

with open('data-biosyn-classes.json', 'w+') as f:
    f.write(json.dumps(sorted(jol, key=lambda data: OFFSET1*int(data.get('item').get('value')[1:]) +\
        abs(py25hash(data.get('p8533'))) + int(data.get('goid')[3:])), indent=0, ensure_ascii=False))

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
