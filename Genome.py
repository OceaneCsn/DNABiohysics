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
import math

class Genome():


    def __init__(self, f=None, indel_size = 60, indelInvRatio = True, T0 = 0.1):
        '''
        Creates a genome.
        The genome is initialized with the files in the folder pathToFiles
        Mettre la signification de tous les attributs aussi
        '''
        pathToFiles = 'D:/ProjetSimADN/Init_files'
        #pathToFiles = '/home/ocassan/ProjetSimADN/Init_files'
        #pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/'
        tmp = pd.read_csv(os.path.join(pathToFiles,'tousgenesidentiques.gff'), header = None, names =  ['name'])
        self.data = pd.DataFrame(tmp.name.str.split('\t',8).tolist(), columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame' , 'Attribute'])
        self.genes = pd.DataFrame(self.data.iloc[4:len(self.data),:])
        self.genes = self.genes.reset_index(drop = True)
        self.size = int(self.genes['End'][0])
        self.header = self.data.iloc[0:4,:]
        #rapport indels/inversions
        self.f = f
        self.T0 = T0
        self.TSS = pd.read_csv(os.path.join(pathToFiles,'TSS.dat'), sep =  '\t')
        self.TTS = pd.read_csv(os.path.join(pathToFiles,'TTS.dat'), sep =  '\t')
        self.barrier = pd.read_csv(os.path.join(pathToFiles,'prot.dat'), sep =  '\t')
        self.env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        self.indel_size = indel_size
        self.indelInvRatio = indelInvRatio

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
        #defines the current genome
        self.gen_ancetre = np.array(list(genome))
        self.fitness = self.compute_fitness(init = True)
        return self.gen_ancetre


    def inversion(self):
        self.gen = np.array(self.gen_ancetre)
        print('inversion pas encore codee')


    def insertion(self):
        '''
        randomly inserts a gene portion of size 60 bp, in the intergenic region
        '''
        self.gen = np.array(self.gen_ancetre)
        intergenes = np.where(self.gen == 'i')[0]
        ins_pos = rd.choice(intergenes)
        print("Insertion at : ", ins_pos)
        self.gen = np.insert(self.gen, ins_pos, np.array(['i' for j in range(self.indel_size)]))

    
    def deletion(self):
        '''
        randomly deletes a gene portion of size 60 bp, in the intergenic region
        and avoiding barrier proteins
        '''
        self.gen = np.array(self.gen_ancetre)
        intergenes = np.where(self.gen == 'i')[0]
        for i in list(g0.barrier["prot_pos"])+list(g0.TSS["TSS_pos"]):
            print(i)
            if i<self.indel_size:
                new_i = self.gen.size - (self.indel_size - i)
                print("new i ", new_i)
                
                to_delete = []
                to_delete2 = []
                for inter in intergenes:
                    if inter in range(new_i, self.gen.size):
                        to_delete.append(inter)
                    
                intergenes = np.delete(intergenes, to_delete)
                
                for inter in intergenes:
                    if inter in range(new_i, self.gen.size):
                        to_delete2.append(inter)

                intergenes = np.delete(intergenes, to_delete2)
            else:
                to_delete = []
                for inter in intergenes:
                    if inter in range(i-self.indel_size, i):
                        to_delete.append(inter)
                intergenes = np.delete(intergenes, to_delete)
        del_pos = rd.choice(intergenes)
        print("Deletion at : ", del_pos)

        if(del_pos > self.gen.size - self.indel_size):
            print("Positions deleted : ", range(del_pos, self.gen.size), " and ", range(0, self.indel_size - (self.gen.size - del_pos)))
            self.gen = np.delete(self.gen, range(del_pos, self.gen.size))
            self.gen = np.delete(self.gen, range(0, self.indel_size - (self.gen.size - del_pos)))
        else:
            print("Positions deleted : ", range(del_pos, del_pos + self.indel_size))
            self.gen = np.delete(self.gen, range(del_pos, del_pos + self.indel_size))
            
            
    def compute_fitness(self, init = False):
        #pathToFiles = '/home/ocassan/ProjetSimADN/'
        pathToFiles = 'D:/ProjetSimADN'
        if init:
            pathToFiles = os.path.join(pathToFiles, 'Init_files')
        new_env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        return np.exp(-sum((new_env.iloc[:,1]-self.env.iloc[:,1])/self.env.iloc[:,1]))
            
    def update_files(self):
        #update barriers positions
        self.barrier['prot_pos'] = np.where(self.genome_ancetre == 'b')[0]

        #update TSS and TTS
        for g in self.genes_list:
            tmp = np.where(self.genome_ancetre==g)[0]
            tss = tmp[0]
            tts = tmp[-1]
            self.TSS.iloc[self.genes_TU[g], 2] = tss
            self.TTS.iloc[self.genes_TU[g], 2] = tts

        #update genes
        self.genes.iloc[0,4] = self.genome_ancetre.size
        self.genes.iloc[1:,3] = list(self.TSS["TSS_pos"])
        self.genes.iloc[1:,4] = list(self.TTS["TTS_pos"])
        
        self.data.iloc[4,4] = self.genome_ancetre.size
        self.data.iloc[5:,3] = list(self.TSS["TSS_pos"])
        self.data.iloc[5:,4] = list(self.TTS["TTS_pos"])


    def dataframes_to_text(self):
        '''
        converts the dataframe attributes to .dat files that can be
        taken as input by the transcription simulator
        '''
        self.TSS.to_csv('TSS.dat', header=True, index=False, sep='\t', mode='w')
        self.TTS.to_csv('TTS.dat', header=True, index=False, sep='\t', mode='w')
        self.barrier.to_csv('prot.dat', header=True, index=False, sep='\t', mode='w')
        self.data.to_csv('tousgenesindentiques.gff', header=False, index=False, sep='\t', mode='w')
        
        
        
    def evolution_step(self, t):
        #creation du genome
        self.create_genome()
        
        #mutation de ce genome suivant la probabilite relative
        r = rd.random()
        if(self.indelInvRatio):
            if  r < self.f:
                if rd.random() < 0.5:
                    self.insertion()
                else:
                    self.deletion()
            else:
                self.inversion()
        else:
            if  r < self.f:
                self.inversion()
            else:
                if rd.random() < 0.5:
                    self.insertion()
                else:
                    self.deletion()
                    
        #evaluation de la fitness du nouveau genome par simulation
        '''utiliser le code du prof de github'''
        new_fitness = self.compute_fitness()
        keep = False
        #garder le nouveau genome ou non
        if new_fitness >= self.fitness:
            keep = True
        else:
            #on tire la probabilite dans une loi exponentielle
            proba = self.T0*math.exp(-self.T0*t)
            if rd.random() < proba:
                keep = True
        if keep:
            #mise a jour des attributs en consequent
            self.genome_ancetre = np.array(self.gen)
            self.fitness = new_fitness
            self.update_files()
            self.dataframes_to_text()
                
    def evolution(self, T):
        '''
        simulates a genome evolution by a Monte Carlo algorithm
        '''
        for t in range(T):
            self.evolution_step(t)


g0 = Genome(f = 0.5)
l = g0.create_genome()
a = g0.evolution_step(0)
a = g0.data
#g0.deletion()
