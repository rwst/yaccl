
import os, json, argparse, sys

"""
Load chemclasses from Wikidata and their superclasses, save in canonical form (Q numbers sorted numerically), output statistics
"""

def qnone2zero(s):
    if s is None:
        return 0
    else:
        return int(s[1:])

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

itsups = {}
if dontquery is False:
    print('performing query...', file=sys.stderr)
    ret = os.popen('wd sparql class-subclass.rq >class-subclass.json'.format(script, script))
    if ret.close() is not None:
        raise
    ff = open('class-subclass.json'.format(script))
    s = ff.read()
    jol = json.loads(s)
    with open('data-class-subclass.json', 'w+') as f:
        f.write(json.dumps(sorted(jol, key=lambda data: OFFSET*(int(data.get('item')[1:]))+qnone2zero(data.get('super'))), indent=0, ensure_ascii=False))
else:
    ff = open('data-class-subclass.json'.format(script))
    s = ff.read()
    jol = json.loads(s)

sups = set()
for d in jol:
    it = d.get('item')
    sup = d.get('super')
    if sup is None:
        exit()
    sups.add(sup)
    i = itsups.get(it)
    if i is None:
        itsups[it] = set(sup)
    else:
        i.add(sup)

print('#items with superclass: {}'.format(len(itsups.keys())))
print('#superclasses: {}'.format(len(sups)))
print('#superclasses not in items: {}'.format(len(sups.difference(set(itsups.keys())))))


