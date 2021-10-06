# yaccl
yet another ChemClassifier (Python based on wikibase-cli and rdk)

This is a proof of the concept of a simple compound classifier that relies completely on knowledge from Wikidata. At this stage it uses InChI keys, SMILES and SMARTS strings, downloaded from Wikidata, to get hits on compounds or classes. About 650 biosynthetic processes from Gene Ontology are directly associated with classes and compounds, so they are potential hits. Matching is done by going through the list of classes---fast enough to find all hits in the dataset with >110k classes within a few seconds.

While yaccl is, in principle, a general classifier, development of patterns focuses on biomolecules.

## Version

The current version is 2112.

## Prerequisites

* Python 3
* rdkit python package
* wikibase-cli (optional, for dataset updates)

Cloning the repo will get you all datasets, so no need to download from Wikidata if you just want to play with yaccl.

## Installation

Clone the repo.

## Usage

Change into the repo `src` directory. There are four scripts. One is for classification, and three are for updating / statistics of datasets.

### classify

Takes a SMILES/InChI and outputs a classification. Examples:

```
$ python3 classify.py -d ./ -m 'CCO'
reading biosyn data
reading class pattern data
reading superclass data
collecting hits...
purging redundant hits
('Q153', 'ethanol', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N') {('GO:0044576', 'pentose catabolic process to ethanol'), ('GO:0044577', 'xylose catabolic process to ethanol'), 'GO:0006115', 'ethanol biosynthetic process', ('GO:0019431', 'acetyl-CoA biosynthetic process from ethanol')}
('Q2832210', 'primary alcohol', '[OX2H1][CX4H2]') {'primary alcohol biosynthetic process', 'GO:0034309'}
```
```
$ python3 classify.py -d ./ -m 'InChI=1S/C15H22O10/c16-3-6-9(19)10(20)11(21)14(23-6)24-13-7-5(1-2-22-13)8(18)12-15(7,4-17)25-12/h1-2,5-14,16-21H,3-4H2/t5-,6-,7-,8+,9-,10+,11-,12+,13+,14+,15-/m1/s1'
reading biosyn data
reading class pattern data
reading superclass data
collecting hits...
purging redundant hits
('Q105151716', '2-[[5-hydroxy-2-(hydroxymethyl)-3,9-dioxatricyclo[4.4.0.02,4]dec-7-en-10-yl]oxy]-6-(hydroxymethyl)oxane-3,4,5-triol', 'OCC1OC(OC2OC=CC3C(O)C4OC4(CO)C23)C(O)C(O)C1O') {'GO:0016138', 'glycoside biosynthetic process'}
('Q408028', 'epoxide', 'C1OC1') {'GO:1901503', 'ether biosynthetic process'}
('Q2832210', 'primary alcohol', '[OX2H1][CX4H2]') {'primary alcohol biosynthetic process', 'GO:0034309'}
('Q2832211', 'secondary alcohol', '[CX4H1][OX2H1]') {'GO:1902653', 'secondary alcohol biosynthetic process'}
('Q74173050', 'D-glucoside', '[C@@H]1(OC([C@H](O)[C@H]([C@@H]1O)O)O*)CO') {'GO:0016138', 'glycoside biosynthetic process'}
('Q74568405', 'alpha-D-mannoside', 'OC[C@@H]1[C@H]([C@@H]([C@H](O)[C@H](O1)O*)O)O') {'GO:0016138', 'glycoside biosynthetic process'}
('Q75078202', 'pyranoside', 'O1C(C(C(C(C1O*)O)O)O)CO') {'GO:0016138', 'glycoside biosynthetic process'}
('Q75082510', 'beta-D-galactoside', 'OC[C@H]1O[C@@H](O[*])[C@H](O)[C@@H](O)[C@H]1O') {'GO:0016138', 'glycoside biosynthetic process'}
('Q105151716', '2-[[5-hydroxy-2-(hydroxymethyl)-3,9-dioxatricyclo[4.4.0.02,4]dec-7-en-10-yl]oxy]-6-(hydroxymethyl)oxane-3,4,5-triol', 'OCC1OC(OC2OC=CC3C(O)C4OC4(CO)C23)C(O)C(O)C1O') {'GO:0016138', 'glycoside biosynthetic process'}
```
Given a file of InChI strings the `-t` option runs a test checking that at least one of the hits is a subclass of the test class item.
```
ralf@ark:~/wikidata/yaccl/src> python3 classify.py -t test/sesterterpenoids.txt -d ./ 
reading biosyn data
reading class pattern data
reading superclass data
classifying testfile test/sesterterpenoids.txt
test class: Q107363222
1 OK
2 OK
3 OK
4 OK
5 OK
6 OK
7 OK
8 OK
9 FAIL: InChI=1S/C25H44/c1-9-21(6)13-11-17-25(8,16-10-12-19(2)3)24-18-22(7)14-15-23(24)20(4)5/h9,18-19,21,23-24H,1,4,10-17H2,2-3,5-8H3
...
```

### datasets
The difference between the `biosyn` and the `class-pattern` datasets is that the former contains classes AND compounds that have associated GO processes in Wikidata (possibly lacking any pattern), while the latter contains all classes (no compounds) that have associated patterns (SMILES, SMARTS or InChI).

The design of datasets will be subject to change in future versions.

### biosyn-class-stats
Output stats of `biosyn` dataset.

```
$ python3 biosyn-class-stats.py 
#items: 1167
#InChI keys: 635
#InChI1 keys: 242
#smarts: 35
#smiles: 695
# w/o pattern: 440
```

### class-pattern
Output stats of `class-pattern` dataset.

```
$ python3 class-pattern.py 
#items: 113788
#InChI keys: 112432
#InChI1 keys: 81120
#smarts: 374
#smiles: 1024
#wildcard smiles: 979
```
Yes, Wikidata has >110k compound classes, most of which (>81k) have unspecified stereochemistry, i.e. they are groups of stereoisomers.

### Contribute
Most errors or inconsistencies are due to data, so please contribute to Wikidata. More involved and pressing is to add patterns to the items of the `biosyn` dataset, or even add missing classes to expand Wikidata's class ontology.

If you want to tackle a specific superclass we recommend usage of `wdtaxonomy` to get an overview of the Wikidata situation.

Of course, any comment or PR is highly appreciated.

### TODO
* option to work with NP subgraph only (should be much faster)
* output JSON
* if there are superclass patterns match them first, narrowing following matches
* write toolserver client
* use RHEA to list enzymatic reactions
* handle iridoid + glycoside --> iridoid glycoside
