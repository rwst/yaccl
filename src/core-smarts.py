import csv, os, json, argparse, sys
from rdkit import Chem

"""
"""

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument('core', type=str, help='an integer for the accumulator')

# Read arguments from the command line
args = parser.parse_args()
#print(args.core.replace('<1>', '[CX4H3,CX3H2]').replace('<2>', '[CX4H2,CX3H1]').replace('<3>', '[CX4H1,CX3H0]').replace('<1a>', '[CX4H3]').replace('<2r>', '[CX4H2,CX3H1]').replace('<3r>', '[CX4H1,CX3H0]').replace('<4r>', '[CR0]').replace('<c>', '[cX3H1]'))
print(args.core.replace('<1>', '[$(C);!$(C(~C)~C)]').replace('<2>', '[$(C);!$(C(~C)(~C)~C)]').replace('<3>', '[$(C);!$(C(~C)(~C)(~C)~C)]').replace('<1a>', '[$(C);!$(C(~C)~C)]').replace('<2r>', '[$([CR0]);!$(C(~C)(~C)~C)]').replace('<3r>', '[$([CR0]);!$(C(~C)(~C)(~C)~C)]').replace('<4r>', '[CR0]').replace('<2c>', '[$(c);!$(c(c)(c)C)]').replace('<3c>', '[$(c);!$(c(c)(c)(C)C)]').replace('<1q>', '[$([#6]);!$([#6](~[#6])~C)]').replace('<2q>', '[$([#6]);!$([#6](~[#6])(~[#6])~C)]').replace('<3q>', '[$([#6]);!$([#6](~[#6])(~[#6])(~[#6])~C)]'))

