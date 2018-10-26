#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 08:52:51 2018

@author: ocassan
"""

import pandas as pd
import os
import numpy as np
import random as rd

class Genome():


    def __init__(self, f=None, first_gen = False, indel_size = 60):
        ## Creates a genome. If this is the first genome of an evolutive process,
        # the genome is initialized with the files in the folder pathToFiles
        # else, the genome will be initialized by copy later
        if first_gen:
            #pathToFiles = '/home/ocassan/ProjetSimADN/'
            pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/'
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
            self.indel_size = indel_size


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
        self.genes_list = list()
        self.genes_TU = {}
        for i in range(1, len(self.genes)):
            start = int(self.genes["Start"][i])
            end = int(self.genes["End"][i])
            gene_name = self.genes["Attribute"][i].split('=')[-1]
            self.genes_list.append(gene_name)
            self.genes_TU[gene_name] = i-1
            genome[start:end+1] = [gene_name] * (end-start+1)
        #add barriers
        for b in self.barrier['prot_pos']:
            genome[b] = 'b'

        self.gen = np.array(list(genome))
        return self.gen


    #def inversion(self):


    def insertion(self):
        intergenes = np.where(self.gen == 'i')[0]
        ins_pos = rd.choice(intergenes)
        print("Insertion at : ", ins_pos)
        self.gen = np.insert(self.gen, ins_pos, np.array(['i' for j in range(self.indel_size)]))

    def update_files(self):
        #update barriers
        self.barrier['prot_pos'] = np.where(self.gen == 'b')[0]

        #update TSS
        for g in self.genes_list:
            tmp = np.where(self.gen==g)[0]
            tss = tmp[0]
            tts = tmp[-1]
            self.TSS.iloc[self.genes_TU[g], 2] = tss
            self.TTS.iloc[self.genes_TU[g], 2] = tts

        #update genes
        self.genes.iloc[0,4] = self.gen.size
        self.genes.iloc[1:,3] = list(self.TSS["TSS_pos"])
        self.genes.iloc[1:,4] = list(self.TTS["TTS_pos"])


    def deletion(self):
        intergenes = np.where(self.gen == 'i')[0]
        for i in list(g0.barrier["prot_pos"])+list(g0.TSS["TSS_pos"]):
            print(i)
            if i<self.indel_size:
                new_i = self.gen.size - (self.indel_size - i)
                print("new i ", new_i)
                intergenes = np.delete(intergenes, range(new_i, self.gen.size))
#A mon avis, ce qu'on fait c'est pas okay : il faudrait prendre le génome copié et ensuite prendre les 'i' et delete ça car dans intergenes l'association index-position est faussée
                intergenes = np.delete(intergenes, range(0, self.indel_size-i))
            else:
                intergenes = np.delete(intergenes, range(i-self.indel_size, i))
        del_pos = rd.choice(intergenes)
        print("Deletion at : ", del_pos)

        if(del_pos > self.gen.size - self.indel_size):
            print("Positions deleted : ", range(del_pos, self.gen.size), " and ", range(0, self.indel_size - (self.gen.size - del_pos)))
            self.gen = np.delete(self.gen, range(del_pos, self.gen.size))
            self.gen = np.delete(self.gen, range(0, self.indel_size - (self.gen.size - del_pos)))
        else:
            print("Positions deleted : ", range(del_pos, del_pos + self.indel_size))
            self.gen = np.delete(self.gen, range(del_pos, del_pos + self.indel_size))

	#def fitness(self):

g0 = Genome(f = 0.5, first_gen=True)
l = g0.create_genome()
print(g0.barrier)
#g0.insertion()
g0.update_files()
print(g0.barrier)
print(g0.gen.shape)

print(g0.genes)
g0.deletion()
g0.update_files()
print(g0.genes)
