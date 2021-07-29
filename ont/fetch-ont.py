
import os, json, argparse, sys

"""
Load classes from Wikidata, save in canonical form (Q numbers sorted numerically)
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
    ret = os.popen('wd sparql {}.rq >{}.json'.format(script, script))
    if ret.close() is not None:
        raise
ff = open('{}.json'.format(script))
s = ff.read()
jol = json.loads(s)
with open('data.json', 'w+') as f:
    f.write(json.dumps(sorted(jol, key=lambda data: OFFSET*int(data.get('item').get('value')[1:]) +\
            int(data.get('super')[3:])), indent=0, ensure_ascii=False))
