===2101===
* initial version

===2102===
* patterns: added >30 SMARTS, replacing unsuitable wildcard SMILES
* biosyn: added 24 SMARTS
* classify: add option to load superclass dataset; does nothing atm
* add `class-subclass.py` for updating superclass dataset

===2103===
* unified dataset handling (all datasets are loaded from given directory)
* data update

===2104===
* classify: remove redundant hits
* classify: get nearest parent GO process if none on class
* classify: add option to get redundant hits as before
* superclass data: more class types, bigger ontology, dead-end diagnostic

===2105===
* all data: monoterpenoids
* test data: monoterpenoids, 14/197 FAIL
* classify: test capability
* multiple SMARTS (denoting sets of patterns) now handled

===2106===
* all data: sesquiterpenoids
* test data: sesquiterpenoids, 7/240 FAIL

===2107===
* all data: diterpenoids
* test data: diterpenoids, 0/177 FAIL
* don't load SMILES/SMARTS if InChI present, 10x speed up of classification

===2108===
* remove bogus pattern leading to Boost warnings
* revised/added to acyclic monoterpenoids, FAIL: 0/197
* all data: sesterterpenoids
* test data: sesterterpenoids, 1/14 FAIL
* split pattern query, makes it faster, more resilient, one third reduction in pattern db size

===2109===
* some new files were missing in last release
* reduce test run output
* all data: triterpenoids
* test data: triterpenoids, 1/85 FAIL

===2110===
* pattern fixes (not all possible double bonds were allowed)
* all data: tetraterpenoids
* test data: tetraterpenoids, 7/199 FAIL

===2111===
* all data: steroids
* test data: steroids, 3/1561 FAIL

===2112===
* add -n option to classify small sample of WD natural products
* all data: flavonoids pt. 1 (est. 60% of all)
* test data: flavonoids

===2113===
* all data: flavonoids pt. 2 (est. 85% of all)

===2114===
* add -j option for JSON output (NP only)
* all data: flavonoids pt. 3 (est. 95% of all)

==2115==
* all data: alkaloids pt. 1 (est. 25% of all)
* test data: alkaloids

==2116==
* fixed bug where compounds with biosyn process were not seen with -j option
* all data: alkaloids pt. 2 (est. 45% of all)

==2117==
* sort pattern data on SMARTS hash to avoid most reshuffling
* all data: alkaloids pt. 3 (est. 62% of all)

==2118==
* fix bug where hits to NP root classes were ignored
* again fix data write order; Py hash() was salted
* add rule to recognize unspecified alkaloids
* all data: alkaloids pt. 4 (est. 35% unspecified)
* add option ouputting NP roots in JSON format
* revamp top NP hierarchy

==2119==
* all data: polyketides pt. 1 (est. 90% of all)
* test data: polyketides
* prettify `--list_nproots`
* further NP root optimizations

==2120==
* test data polyketides was missing
* add rule to recognize unspecified macrolides
* all data: polyketides pt. 2 (est. 97% of all)

==2121==
* classify: revamp -n code
* classify: include mol data with JSON output
* improve macrolide rule
* reduced unspec. alkaloids to 25% by adding more small groups
* unspec. alkaloids and makrolides tagged in WD

==2122==
* all data: fatty acyl pt. 1 (est. 77% of all)
* test data: fatty acyl
* classify: fatty acyl rules


