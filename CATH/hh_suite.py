#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
import pdb

#Other script imports
from hh_reader import read_result


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs hhblitz to create HMMs of domain sequences, then runs hhalign on all
						pairs and checks that over 75 % has been aligned using hh_reader.py''')
 
parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to directory with input files.')








read_result


