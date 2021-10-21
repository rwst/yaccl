In place of a detailed outline of the strategy for creating SMARTS patterns for classification (which doesn't exist), the following at least lists common features that we have used in constructing the patterns that constitute the core of `yaccl`.
This page can be considered an extension of such introductory pages by Daylight as:
* https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
* https://www.daylight.com/dayhtml_tutorials/languages/smarts/
* https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

The patterns are all working in `rdkit`.

SMILES | comment
--- | ---
`C` | aliphatic C
`c` | aromatic C
`[#6]` | any C
`[CR0]` | aliphatic C not in ring
`[CR1]` | aliphatic C in exactly one ring
`[cr5]` | aromatic C in 5-member ring
`[Cr;!r3;!r4;!r5;!r6;!r7]` | aliphatic C in macrocycle (>7)
`[CX4H1]` | C bound to 4 atoms, of which one is H
--- | ---
`[$(C);$(CO)]` | C bound to O---this only describes and matches the C!
`[$(C);$([CX4H1](C)(C)O),$([CX4H2](C)C)]` | C bound either to C,C,O,H or to C,C,H,H
`[$(C);!$(C(~C)(~C)(~C)~C)]` | C not bound to 4 other C
