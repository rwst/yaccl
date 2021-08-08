# yaccl
yet another ChemClassifier (Python based on wikibase-cli and rdk)

This is a proof of the concept of a simple compound classifier that relies completely on knowledge from Wikidata. At this stage it uses InChI keys, SMILES and SMARTS strings, downloaded from Wikidata, to get hits on compounds or classes. About 650 biosynthetic processes from Gene Ontology are directly associated with classes and compound, so they are potential hits. Matching is done by going through the list of classes---fast enough to find all hits in the dataset with >110k classes within a few seconds.

## Version

The current version is 2101.

## Prerequisites

* Python 3
* rdkit python package
* wikibase-cli (optional, for dataset updates)

Cloning the repo will get you all datasets, so no need to download from Wikidata if you just want to play with yaccl.

## Installation

Clone the repo.

## Usage

Change into the repo `src` directory. There are three scripts. One is for classification, and two are for updating / statistics of datasets.

### classify

Takes a SMILES/InChI and outputs a classification. Examples for each dataset:

```
$ python3 classify.py -d data-biosyn.json -m 'CCO'
Q153 ethanol LFQSCWFLJHTTHZ-UHFFFAOYSA-N {'ethanol biosynthetic process', 'GO:0006115', ('GO:0019431', 'acetyl-CoA biosynthetic process from ethanol'), ('GO:0044576', 'pentose catabolic process to ethanol'), ('GO:0044577', 'xylose catabolic process to ethanol')}
Q153 ethanol CCO {'ethanol biosynthetic process', 'GO:0006115', ('GO:0019431', 'acetyl-CoA biosynthetic process from ethanol'), ('GO:0044576', 'pentose catabolic process to ethanol'), ('GO:0044577', 'xylose catabolic process to ethanol')}
Q156 alcohols O[*] {'GO:0046165', 'alcohol biosynthetic process'}
Q2832210 primary alcohol [OX2H1][CX4H2] {'primary alcohol biosynthetic process', 'GO:0034309'}
```
```
$ python3 classify.py -d data-class-pattern.json -m 'InChI=1S/C15H22O10/c16-3-6-9(19)10(20)11(21)14(23-6)24-13-7-5(1-2-22-13)8(18)12-15(7,4-17)25-12/h1-2,5-14,16-21H,3-4H2/t5-,6-,7-,8+,9-,10+,11-,12+,13+,14+,15-/m1/s1'
reading data
Q105151716 catalpol LHDWRKICQLTVDL {None}
Q156 alcohols O[*] {None}
Q103230 ethers [#6!$([#6]=[!#6])][OX2][#6!$([#6]=[!#6])] {None}
Q408028 epoxide C1OC1 {None}
Q416840 iridoid [$(C);!$(C(~[#6])~[#6])]~[$(C);!$(C(~C)(~C)(~C)~C)]1[$(C);!$(C(~C)(~C)~C)]~[$(C);!$(C(~C)(~C)~C)]~[$(C);!$(C(~C)(~C)(~C)~C)]2~[$(C);!$(C(~C)(~C)(~C)~C)]~1~[$(C);!$(C(~C)~C)]O[$(C);!$(C(~C)~C)]~C~2 {None}
Q2832210 primary alcohol [OX2H1][CX4H2] {None}
Q2832211 secondary alcohol [CX4H1][OX2H1] {None}
Q74173050 D-glucoside [C@@H]1(OC([C@H](O)[C@H]([C@@H]1O)O)O*)CO {None}
Q74568405 alpha-D-mannoside OC[C@@H]1[C@H]([C@@H]([C@H](O)[C@H](O1)O*)O)O {None}
Q75078202 pyranoside O1C(C(C(C(C1O*)O)O)O)CO {None}
Q75082510 beta-D-galactoside OC[C@H]1O[C@@H](O[*])[C@H](O)[C@@H](O)[C@H]1O {None}
Q105151716 catalpol OCC1OC(OC2OC=CC3C(O)C4OC4(CO)C23)C(O)C(O)C1O {None}
```
### datasets
The difference between the `biosyn` and the `class-pattern` datasets is that the former contains classes AND compounds that have associated GO processes in Wikidata (possibly lacking any pattern), while the latter contains all classes (no compounds) that have associated patterns (SMILES, SMARTS or InChI). Note that the catalpol (Q105151716) hit in the second example is actually a class (group of stereoisomers), not a compound.

The design of datasets will be subject to change in future versions.

### biosyn-class-stats
Output stats of `biosyn` dataset.

```
$ python3 biosyn-class-stats.py 
#items: 1168
#InChI keys: 634
#InChI1 keys: 241
#smarts: 13
#smiles: 695
# w/o pattern: 463
```

### class-pattern
Output stats of `biosyn` dataset.

```
$ python3 class-pattern.py 
#items: 113486
#InChI keys: 112388
#InChI1 keys: 81085
#smarts: 35
#smiles: 104011
#wildcard smiles: 994
```
Yes, Wikidata has >110k compound classes, most of which (>81k) have unspecified stereochemistry, i.e. they are groups of stereoisomers.

### Contribute
Most errors or inconsistencies are due to data, so please contribute to Wikidata. For example catalpol (Q105151716) should be renamed 2-[[5-hydroxy-2-(hydroxymethyl)-3,9-dioxatricyclo[4.4.0.02,4]dec-7-en-10-yl]oxy]-6-(hydroxymethyl)oxane-3,4,5-triol because the stereochemistry is unspecified. More involved and pressing is to add patterns to the items of the `biosyn` dataset, or even add missing classes to expand Wikidata's class ontology.

Of course, any comment or PR is highly appreciated.
