#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 08:52:51 2018

@author: ocassan
"""

import pandas as pd
import os

class Genome():
    
    
    def __init__(self, f=None, first_gen = False):
        ## Creates a genome. If this is the first genome of an evolutive process, 
        # the genome is initialized with the files in the folder pathToFiles
        # else, the genome will be initialized by copy later
        if first_gen:
            pathToFiles = '/home/ocassan/ProjetSimADN/'
            tmp = pd.read_csv(os.path.join(pathToFiles,'tousgenesidentiques.gff'), header = None, names =  ['name'])
            self.data = pd.DataFrame(tmp.name.str.split('\t',8).tolist(), columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame' , 'Attribute'])
            self.genes = pd.DataFrame(self.data.iloc[4:len(self.data),:])
            self.genes = self.genes.reset_index(drop = True)
            self.size = int(self.genes['End'][0])
            self.header = self.data.iloc[0:4,:]
            self.f = f
            self.TSS = pd.read_csv(os.path.join(pathToFiles,'TSS.dat'), sep =  '\t')
            self.TTS = pd.read_csv(os.path.join(pathToFiles,'TTS.dat'), sep =  '\t')
            self.barrier = pd.read_csv(os.path.join(pathToFiles,'prot.dat'), sep =  '\t')
            self.env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
    
    
    def copy(self, gen):
        # sets the genome atributes with the atributes of the genome given in
        #argument
        self.genes = gen.genes
        self.f = gen.f
        self.TSS = gen.TSS
        self.TTS = gen.TTS
        self.barrier = gen.barrier
        self.env = gen.env
        self.size = gen.size
        self.header = gen.header
    
    def create_genome(self):
        #create a linear genome
        genome = ['i' for i in range(self.size)]        
        for i in range(1, len(self.genes)):
            start = int(self.genes["Start"][i])
            end = int(self.genes["End"][i])
            gene_name = self.genes["Attribute"][i].split('=')[-1]
            genome[start:end+1] = [gene_name] * (end-start+1)
        self.gen = list(genome)
        return genome
            

    #def inversion(self):
        
        
    #def insertion(self):
        
    #def deletion(self):
        
    #def fitness(self):
        
g0 = Genome(f = 0.5, first_gen=True)        
g = Genome()
g.copy(g0)
l = g0.create_genome()
