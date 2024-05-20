#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:01:49 2024

@author: frulcino
"""

from gurobi_optimods import datasets, opf


case = datasets.load_opf_example("case9")
result = opf.solve_opf(case, opftype="AC")
