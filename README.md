# yaccl
yet another ChemClassifier (Python based on wikibase-cli and rdk)

This is a proof of the concept of a simple compound classifier relying completely on knowledge from Wikidata. At this stage it uses InChI keys, SMILES and SMARTS strings, downloaded from Wikidata, to get hits on compounds or classes, and about 650 associated biosynthetic processes from Gene Ontology. Matching is done by going through the list of classes---fast enough to find all hits in the dataset with >110k classes within a few seconds.
