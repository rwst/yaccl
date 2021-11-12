
import os, json, argparse, sys

"""
Load chemclasses from Wikidata together with any pattern, save in canonical form (Q numbers sorted numerically), output statistics
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

# Read arguments from the command line
args = parser.parse_args()

# Check for --version or -V
dontquery = not args.query
script = os.path.basename(sys.argv[0])[:-3]
OFFSET = 1000000000

if dontquery is False:
    print('performing query 1...', file=sys.stderr)
    ret = os.popen('wd sparql class-pattern1.rq >class-pattern1.json'.format(script, script))
    if ret.close() is not None:
        raise
ff = open('class-pattern1.json'.format(script))
s = ff.read()
jol = json.loads(s)

items = set()
iks = {}
ik1 = {}
smarts = {}
smiles = {}
outj = []
for d in jol:
    outj.append(d)
    it = d.get('item').get('value')
    items.add(it)
    ik = d.get('p235')
    if ik is not None:
        iks[ik] = it
        if ik[15:25] == 'UHFFFAOYSA':
            ik1[ik[:14]] = it

if dontquery is False:
    print('performing query 2...', file=sys.stderr)
    ret = os.popen('wd sparql class-pattern2.rq >class-pattern2.json'.format(script, script))
    if ret.close() is not None:
        raise
ff = open('class-pattern2.json'.format(script))
s = ff.read()
jol = json.loads(s)

for d in jol:
    outj.append(d)
    it = d.get('item').get('value')
    items.add(it)
    sma = d.get('p8533')
    if sma is not None:
        smarts[sma] = it
    csm = d.get('p233')
    ism = d.get('p2017')
    smi = None
    if ism is not None:
        smi = ism
    elif csm is not None:
        smi = csm
    if smi is not None:
        smiles[smi] = it

with open('data-class-pattern.json', 'w+') as f:
    f.write(json.dumps(sorted(outj, key=lambda data: 2**64*int(data.get('item').get('value')[1:]) + abs(py25hash(data.get('p8533')))), indent=0, ensure_ascii=False))

print('#items: {}'.format(len(items)))
print('#InChI keys: {}'.format(len(iks.keys())))
print('#InChI1 keys: {}'.format(len(ik1.keys())))
print('#smarts: {}'.format(len(smarts.keys())))
print('#smiles: {}'.format(len(smiles.keys())))
print('#wildcard smiles: {}'.format(sum(sm.find('*') >= 0 for sm in smiles.keys())))


