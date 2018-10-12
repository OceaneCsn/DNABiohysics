#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 08:52:51 2018

@author: ocassan
"""

import pandas as pd
import os

class Genome():
    
    def __init__(self, f, first_gen = False):
        
        if first_gen:
            pathToFiles = '/home/ocassan/ProjetSimADN/'
            tmp = pd.read_csv(os.path.join(pathToFiles,'tousgenesidentiques.gff'), header = None, names =  ['name'])
            self.genes = pd.DataFrame(tmp.name.str.split('\t',8).tolist(), columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame' , 'Attribute'])
            self.f = f
            
            self.TSS = pd.read_csv(os.path.join(pathToFiles,'TSS.dat'), sep =  '\t')
            self.TTS = pd.read_csv(os.path.join(pathToFiles,'TTS.dat'), sep =  '\t')
            self.barrier = pd.read_csv(os.path.join(pathToFiles,'prot.dat'), sep =  '\t')
            self.env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        
    def copy(self, gen):
        self.genes = gen.genes
        self.f = gen.f
        self.TSS = gen.TSS
        self.TTS = gen.TTS
        self.barrier = gen.barrier
        self.env = gen.env
        
    #def inversion(self):
        
        
    #def insertion(self):
        
    #def deletion(self):
        
    #def fitness(self):
        
g0 = Genome(0.5, first_gen=True)        
g = Genome(0.5)
g.copy(g0)
g.genes
g.env
g.barrier
g.TSS
