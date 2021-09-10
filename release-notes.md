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

===2107==
* all data: diterpenoids
* test data: diterpenoids, 0/177 FAIL
* don't load SMILES/SMARTS if InChI present, 10x speed up of classification

