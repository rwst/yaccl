
import csv, os, json, argparse, sys

"""
Load classes from Wikidata, save in canonical form (Q numbers sorted numerically)
"""
# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="perform SPARQL query",
        action="store_true")

# Read arguments from the command line
args = parser.parse_args()
#print(args)

# Check for --version or -V
dontquery = not args.query
script = os.path.basename(sys.argv[0])[:-3]
OFFSET = 1000000000

if dontquery is False:
    print('performing query...', file=sys.stderr)
    ret = os.popen('wd sparql {}.rq >{}.json'.format(script, script))
    if ret.close() is not None:
        raise
ff = open('{}.json'.format(script))
s = ff.read()
jol = json.loads(s)
with open('data.json', 'w+') as f:
    f.write(json.dumps(sorted(jol, key=lambda data: OFFSET*int(data.get('item').get('value')[1:]) +\
            int(data.get('super')[3:])), indent=0, ensure_ascii=False))


exit()

sitems = {}
labels = {}
for d in jol:
    dd = d.get('item')
    it = dd.get('value')
    lab = dd.get('label')
    sup = d.get('super')
    go = d.get('goid')
    p233 = d.get('p233')
    p2017 = d.get('p2017')
    p8533 = d.get('p8533')
    i = sitems.get(it)
    if i is not None:
        i.append(sup)
    else:
        sitems[it] = [sup]
    labels[it] = lab
sitems['Q2393187'] = 'Q43460564'
labels['Q2393187'] = 'molecular entity'
labels['Q43460564'] = 'chemical entity'

#print('Q415812' in set(items.keys()).union(set(['Q43460564'])))
#print([it for it,itsuplist in items.items() if 'Q415812' in itsuplist ])

edges = {}
for it,itsuplist in items.items():
    for sup in itsuplist:
        if sup in set(items.keys()).union(set(['Q43460564'])):
            e = edges.get(sup)
            if e is None:
                edges[sup] = set([it])
            else:
                e.add(it)
#print(edges.get('Q415812'))

seen = set()
def walk(E, edges, prefix):
    pfix = '├──'
    if E in seen:
        pfix = '╞══'
    print('{}[[{}]] {}'.format(prefix + pfix, E, labels.get(E)))
    if E in seen and edges.get(E) is not None:
        print(prefix + '... see above')
        return
    seen.add(E)
    children = edges.get(E)
    prefix = ' │   ' + prefix
    if len(prefix) > 70 or children is None:
        return
    for c in sorted(children, key=lambda c: labels.get(c)):
        walk(c, edges, prefix)

walk('Q43460564', edges, ' ')
print()
walk('Q36496', edges, ' ')

