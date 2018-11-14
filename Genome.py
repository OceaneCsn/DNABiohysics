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
        #pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/Init_files/'
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

        #selection of two inversion bounds in intergenic region
        intergenes = np.where(self.gen == 'i')[0]
        pos1 = rd.choice(intergenes)
        pos2 = rd.choice(intergenes)
        print('inversion de ', min(pos1,pos2), ' a ', max(pos1,pos2))

        #inversion des nucleotides
        inverted = np.flip(self.gen[min(pos1,pos2) : max(pos1, pos2)])
        print(inverted.size, inverted)
        list_gen = self.gen.tolist()
        list_gen[min(pos1,pos2) : max(pos1, pos2)] = inverted.tolist()
        self.gen = np.array(list_gen)

        #changement du sens des genes dans le dataframe self.genes
        for g in self.genes_list:
            if g in self.gen[min(pos1,pos2) : max(pos1, pos2)]:
                print(g, 'inverted')
                #on regarde ou il faut changer le signe du strand dans self.genes
                #print('index ',i, ' tu ', self.genes_TU[g])
                if self.genes['Strand'].iloc[self.genes_TU[g]+1,] == '+':
                    self.genes['Strand'].iloc[self.genes_TU[g]+1,] = '-'
                else:
                    self.genes['Strand'].iloc[self.genes_TU[g]+1,] = '+'

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
        #Get all the non intergenes indices
        not_intergenes = np.where(self.gen != 'i')[0]
        genes_name = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10']
        gene_length = 1000

		#Only keep the beginning of the genes
        get_genes_starts = []
        for i in genes_name:
            get_genes_starts.append(list(self.gen).index(i))
        for i in get_genes_starts:
            to_delete = [j for j in range(i + 1, i + gene_length)]
            not_intergenes = not_intergenes[np.invert(np.isin(not_intergenes, to_delete))]
        print(not_intergenes)

		#Apply a protection field before each non-intergene element : between 
        #this number and the correpsonding one of not_intergenes, 
        #none can be deleted
        to_protect_bounds = not_intergenes - self.indel_size
        for i in to_protect_bounds:
            if i < 0:
                to_protect_bounds = np.append(to_protect_bounds, [self.gen.size + i, 0])
                to_protect_bounds = to_protect_bounds[to_protect_bounds >= 0]
        to_protect_bounds.sort()
        print(to_protect_bounds)
        not_possible = []
        for i in zip(to_protect_bounds, not_intergenes):
            not_possible.extend([j for j in range(i[0], i[1])])
        for i in get_genes_starts:
            not_possible.extend([j for j in range(i, i + gene_length)])
        not_possible.sort()

		#Set the indices where the del_pos can be
        possible = np.array([i for i in range(0, self.gen.size)])
        possible = possible[np.invert(np.isin(possible, not_possible))]
        del_pos = rd.choice(possible)
        print("del_pos = ", del_pos, " Is del_pos in not_possible ?", del_pos in not_possible)

        #Remove
        if(del_pos > self.gen.size - self.indel_size):
            print("Positions deleted : ", range(del_pos, self.gen.size), " and ", range(0, self.indel_size - (self.gen.size - del_pos)))
            self.gen = np.delete(self.gen, range(del_pos, self.gen.size))
            self.gen = np.delete(self.gen, range(0, self.indel_size - (self.gen.size - del_pos)))
        else:
            print("Positions deleted : ", range(del_pos, del_pos + self.indel_size))
            self.gen = np.delete(self.gen, range(del_pos, del_pos + self.indel_size))

        '''
        intergenes = np.where(self.gen == 'i')[0]
        for i in list(g0.barrier["prot_pos"])+list(g0.TSS["TSS_pos"]):
            #print(i)
            if i<self.indel_size:
                new_i = self.gen.size - (self.indel_size - i)

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
        '''

        return self.gen

    def compute_fitness(self, init = False):
        #pathToFiles = '/home/ocassan/ProjetSimADN/'
        #pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/'
        pathToFiles = 'D:/ProjetSimADN'
        if init:
            pathToFiles = os.path.join(pathToFiles, 'Init_files')
        new_env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        return np.exp(-sum((new_env.iloc[:,1]-self.env.iloc[:,1])/self.env.iloc[:,1]))

    def update_files(self):
        #update barriers positions
        self.barrier['prot_pos'] = np.where(self.gen_ancetre == 'b')[0]

        #update TSS and TTS
        for g in self.genes_list:
            tmp = np.where(self.gen_ancetre==g)[0]
            tss = tmp[0]
            tts = tmp[-1]
            if self.genes['Strand'].iloc[self.genes_TU[g]+1,] == '-':
                temp = tss
                tss = tts
                tts = temp
            self.TSS.iloc[self.genes_TU[g], 2] = tss
            self.TTS.iloc[self.genes_TU[g], 2] = tts

        #update genes size, tss and tts
        self.genes.iloc[0,4] = self.gen_ancetre.size
        self.genes.iloc[1:,3] = list(self.TSS["TSS_pos"])
        self.genes.iloc[1:,4] = list(self.TTS["TTS_pos"])

        #orientation dans TSS et TTS
        self.TSS['TUorient'] = list(self.genes['Strand'].iloc[1:,])
        self.TTS['TUorient'] = list(self.genes['Strand'].iloc[1:,])

        #tout mis a jour dans data
        self.data.iloc[4,4] = self.gen_ancetre.size
        self.data.iloc[5:,3] = list(self.TSS["TSS_pos"])
        self.data.iloc[5:,4] = list(self.TTS["TTS_pos"])
        self.data['Strand'].iloc[5:,] = list(self.genes['Strand'].iloc[1:,])

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

        #mutation du genome suivant la probabilite relative
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
        '''
        ICI : utiliser le code du prof de github qui va generer le nouvel
        environnement.dat utilise par compute_fitness
        '''
        new_fitness = self.compute_fitness()
        keep = False
        #garder le nouveau genome?
        if new_fitness > self.fitness:
            keep = True
        else:
            #on tire la probabilite dans une loi exponentielle
            proba = self.T0*math.exp(-self.T0*t)
            if rd.random() < proba:
                keep = True
        if True:
            #mise a jour des attributs en consequent
            self.gen_ancetre = np.array(self.gen)

            self.fitness = new_fitness
            self.update_files()
            self.dataframes_to_text()


    def evolution(self, T):
        '''
        simulates a genome evolution using a Monte Carlo algorithm
        '''
        self.create_genome()
        for t in range(T):
            self.evolution_step(t)


g0 = Genome(f = 0.8)

g0.evolution(20)

a = g0.genes

