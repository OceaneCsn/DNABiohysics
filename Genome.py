#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 08:52:51 2018

@authors: Julie Etienne, Amaury Prin, Océane Cassan

Projet de biologie computationnelle

L’objectif de projet est de développer un code de simulation d’évolution sur un 
génome dans lequel nous allons effectuer des modifications évolutives dans 
l’objectif de tester l’adaptation de réseaux de régulation transcriptionnels.

"""
import pandas as pd
import os
import numpy as np
import random as rd
import math
import simulation2 as sim
import pickle

class Genome():

    def __init__(self, pathToFiles = '/home/ocassan/ProjetSimADN', nb_pol = None, f=0.8, indel_size = 60, indelInvRatio = True, T0 = 0.1):
        '''
        Creates a genome.
        The genome is initialized with the files in the folder pathToFiles
        Mettre la signification de tous les attributs aussi
        '''
        self.pathToFiles = pathToFiles
        pathToInitFiles = os.path.join(self.pathToFiles, 'Init_files')
        pathToInitFiles = 'Init_files'
        tmp = pd.read_csv(os.path.join(pathToInitFiles,'tousgenesidentiques.gff'), header = None, names =  ['name'])
        #dataframe attributes of the genome, and parameters
        self.data = pd.DataFrame(tmp.name.str.split('\t',8).tolist(), columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame' , 'Attribute'])
        self.genes = pd.DataFrame(self.data.iloc[4:len(self.data),:])
        self.genes = self.genes.reset_index(drop = True)
        self.size = int(self.genes['End'][0])
        self.header = self.data.iloc[0:4,:]
        #rapport indels/inversions
        self.f = f
        self.T0 = T0
        self.TSS = pd.read_csv(os.path.join(pathToInitFiles,'TSS.dat'), sep =  '\t')
        self.TTS = pd.read_csv(os.path.join(pathToInitFiles,'TTS.dat'), sep =  '\t')
        self.barrier = pd.read_csv(os.path.join(pathToInitFiles,'prot.dat'), sep =  '\t')
        self.env = pd.read_csv(os.path.join(pathToInitFiles,'environment.dat'), header = None, sep =  '\t') #environment of reference : the goal bacteria wants to achieve to survive
        self.indel_size = indel_size
        self.indelInvRatio = indelInvRatio
        self.nb_pol = nb_pol

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
        self.gen = np.array(self.gen_ancetre)
        self.update_files(self.genes)
        self.dataframes_to_text()
        if self.nb_pol == None:
            sim.start_transcribing(os.path.join(self.pathToFiles,'paramsInit.ini'),
                                   os.path.join(self.pathToFiles, 'testRes'))
        else:
            sim.start_transcribing(os.path.join(self.pathToFiles,'params.ini'),
                                   os.path.join(self.pathToFiles, 'testRes'), nb_pol = int(self.nb_pol))
        self.fitness = self.compute_fitness()
        #self.fitness = 1
        self.fitnesses.append(self.fitness)
        self.tmp_genes = self.genes.copy()
        return self.gen_ancetre


    def inversion(self):
        self.gen = np.array(self.gen_ancetre)

        #selection of two inversion bounds in intergenic region
        intergenes = np.where(self.gen == 'i')[0]
        pos1 = rd.choice(intergenes)
        pos2 = rd.choice(intergenes)
        print('Inversion de ', min(pos1,pos2), ' a ', max(pos1,pos2))

        #inversion des nucleotides
        inverted = np.flip(self.gen[min(pos1,pos2) : max(pos1, pos2)], axis=0)
        list_gen = self.gen.tolist()
        list_gen[min(pos1,pos2) : max(pos1, pos2)] = inverted.tolist()
        self.gen = np.array(list_gen)

        #changement du sens des genes dans le dataframe self.genes

        for g in self.genes_list:
            if g in self.gen[min(pos1,pos2) : max(pos1, pos2)]:
                #print(g, 'inverted')
                #on regarde ou il faut changer le signe du strand dans self.genes
                #print('index ',i, ' tu ', self.genes_TU[g])
                if self.genes['Strand'].iloc[self.genes_TU[g]+1,] == '+':
                    self.tmp_genes['Strand'].iloc[self.genes_TU[g]+1,] = '-'
                else:
                    self.tmp_genes['Strand'].iloc[self.genes_TU[g]+1,] = '+'
        return max(pos1, pos2) - min(pos1, pos2)

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

		#Apply a protection field before each non-intergene element : between
        #this number and the correpsonding one of not_intergenes,
        #none can be deleted
        to_protect_bounds = not_intergenes - self.indel_size*2
        for i in to_protect_bounds:
            if i < 0:
                to_protect_bounds = np.append(to_protect_bounds, [self.gen.size + i, 0])
                to_protect_bounds = to_protect_bounds[to_protect_bounds >= 0]
        to_protect_bounds.sort()
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
        #print("del_pos = ", del_pos, " Is del_pos in not_possible ?", del_pos in not_possible)

        #Remove
         
        if(del_pos > self.gen.size - self.indel_size):
            maxi = self.indel_size - (self.gen.size - del_pos)
            for b in np.where(self.gen == 'b')[0]:
                if b in range(0, self.indel_size - (self.gen.size - del_pos)):
                    maxi = b-1
            #print("Positions deleted : ", range(del_pos, self.gen.size), " and ", range(0, self.indel_size - (self.gen.size - del_pos)))
            self.gen = np.delete(self.gen, range(del_pos, self.gen.size))
            self.gen = np.delete(self.gen, range(0, maxi))
            
                    
            print("deletion at ", del_pos, 'entre ', del_pos, ' et ', self.indel_size - (self.gen.size - del_pos) )
        else:
            #print("Positions deleted : ", range(del_pos, del_pos + self.indel_size))
            mini = del_pos
            for b in np.where(self.gen == 'b')[0]:
                if b in range(del_pos, del_pos + self.indel_size+1):
                    print('bariere supprimee!!!!!!!!!!!!!!!!!!!!!!', b, 'entre ', del_pos, ' et ', del_pos+self.indel_size )
                    if b == mini:
                        mini = b+1
            self.gen = np.delete(self.gen, range(mini, del_pos + self.indel_size))
            print("deletion at ", del_pos, 'entre ', del_pos, ' et ', del_pos+self.indel_size )
        return self.gen

    def compute_fitness(self, init = False):

        if init:
            self.pathToFiles = os.path.join(self.pathToFiles, 'Init_files')
            
        else:
            self.pathToFiles = self.pathToFiles

        #new_env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        total_nb_transcrits = self.nb_transcrits().sum()[0]
        new_env = self.nb_transcrits() * 1.0/total_nb_transcrits
        #new_env.iloc[:,1] = self.nb_transcrits() * 1.0/total_nb_transcrits
        #new_env.to_csv('environment.dat', header = None, index = False, sep = '\t', mode = 'w')
        #new_env = pd.read_csv(os.path.join(pathToFiles,'environment.dat'), header = None, sep =  '\t')
        #new_env = self.nb_transcrits() * 1.0/total_nb_transcrits
        #return np.exp(-sum((new_env.iloc[:,1]-self.env.iloc[:,1])/self.env.iloc[:,1]))
        #print("The env of the bacteria :\n", new_env)
        #print("The reference :\n", self.env)
        return np.exp(-sum((new_env.iloc[:,0]-self.env.iloc[:,1])/self.env.iloc[:,1]))

    def update_files(self, genesDf):
        #update barriers positions
        #print(len(np.where(self.gen == 'b')[0]), ' barrières : ', np.where(self.gen == 'b')[0])
        self.barrier['prot_pos'] = np.where(self.gen == 'b')[0]

        #update TSS and TTS
        for g in self.genes_list:
            tmp = np.where(self.gen==g)[0]
            tss = tmp[0]
            tts = tmp[-1]
            if genesDf['Strand'].iloc[self.genes_TU[g]+1,] == '-':
                temp = tss
                tss = tts
                tts = temp
            self.TSS.iloc[self.genes_TU[g], 2] = tss
            self.TTS.iloc[self.genes_TU[g], 2] = tts

        #update genes size, tss and tts
        genesDf.iloc[0,4] = self.gen.size
        genesDf.iloc[1:,3] = list(self.TSS["TSS_pos"])
        genesDf.iloc[1:,4] = list(self.TTS["TTS_pos"])

        #orientation dans TSS et TTS
        self.TSS['TUorient'] = list(genesDf['Strand'].iloc[1:,])
        self.TTS['TUorient'] = list(genesDf['Strand'].iloc[1:,])

        genesDf.iloc[0,3] = int(genesDf.iloc[0,3])
        genesDf = genesDf.sort_values(by='Start')
        self.TSS = self.TSS.sort_values(by=['TSS_pos'])
        self.TTS = self.TTS.sort_values(by=['TTS_pos'])
        #print(self.TSS['TSS_pos'], self.TTS['TTS_pos'])
        #self.barrier.sort_index(by = ['prot_pos'])
        #tout mis a jour dans data
        self.data.iloc[4,4] = self.gen.size
        self.data.iloc[5:,3] = list(self.TSS["TSS_pos"])
        self.data.iloc[5:,4] = list(self.TTS["TTS_pos"])
        self.data['Strand'].iloc[5:,] = list(genesDf['Strand'].iloc[1:,])

    def dataframes_to_text(self):
        '''
        converts the dataframe attributes to .dat files that can be
        taken as input by the transcription simulator
        '''
        self.TSS.to_csv('TSS.dat', header=True, index=False, sep='\t', mode='w')
        self.TTS.to_csv('TTS.dat', header=True, index=False, sep='\t', mode='w')
        self.barrier.to_csv('prot.dat', header=True, index=False, sep='\t', mode='w')
        self.data.to_csv('tousgenesidentiques.gff', header=False, index=False, sep='\t', mode='w')


    def evolution_step(self, t):
        #mutation du genome suivant la probabilite relative
        r = rd.random()
        inv_size = 0
        if(self.indelInvRatio):
            if  r < self.f:
                if rd.random() < 0.5:
                    event = 'insertion'
                    self.insertion()
                else:
                    self.deletion()
                    event = 'deletion'
            else:
                inv_size = self.inversion()
                event = 'inversion'
        else:
            if  r < self.f:
                inv_size = self.inversion()
                event = 'inversion'
            else:
                if rd.random() < 0.5:
                    event = 'insertion'
                    self.insertion()
                else:
                    self.deletion()
                    event = 'deletion'

        #evaluation de la fitness du nouveau genome par simulation        
        self.update_files(self.tmp_genes)
        self.dataframes_to_text()
        try:
            keep = False
            if self.nb_pol == None:
                sim.start_transcribing(os.path.join(self.pathToFiles,'params.ini'),
                                       os.path.join(self.pathToFiles, 'testRes'))
            else:
                sim.start_transcribing(os.path.join(self.pathToFiles,'params.ini'),
                                       os.path.join(self.pathToFiles, 'testRes'), nb_pol = int(self.nb_pol))
            new_fitness = self.compute_fitness()
            #new_fitness = 1
            
            #garder le nouveau genome?
            if new_fitness > self.fitness:
                keep = True
            else:
                #on tire la probabilite dans une loi exponentielle
                proba = self.T0*math.exp(-self.T0*t)
                if rd.random() < proba:
                    keep = True
            #print("old fit : ", self.fitness, " new fi : ", new_fitness, " keep = ", keep)
            if event == 'inversion':
                self.delta_fitness[event].append((self.fitness-new_fitness, inv_size))
            else:
                self.delta_fitness[event].append(self.fitness-new_fitness)
            if keep:
                #print("keep = ", keep)
                #mise a jour des attributs en consequent
                self.genes = self.tmp_genes.copy()
                self.gen_ancetre = np.array(self.gen)
                self.fitness = new_fitness
                self.events.append((t,event))
            self.fitnesses.append(self.fitness)
        
        except ValueError:
            print('Genes cotes a cotes dans simulation....')
            #plt.figure()
            #plt.plot(self.gen)

    #Function to get the number of transcrits per gene
    def nb_transcrits(self):
        pathToResFiles = os.path.join(self.pathToFiles, 'testRes')
        transcrits = pd.read_csv(os.path.join(pathToResFiles, 'save_tr_nbr.csv'), header = None)
        return transcrits

    def evolution(self, T, dump = False):
        '''
        simulates a genome evolution using a Monte Carlo algorithm
        '''
        self.fitnesses = []
        self.events = []
        self.delta_fitness = {ev:[] for ev in ['insertion', 'deletion', 'inversion']}
        self.create_genome()
        for t in range(T):
            self.evolution_step(t)
        if dump:
            path = os.path.join(self.pathToFiles, 'Binary_files')
            with open(os.path.join(path, 'genome.file'), "wb") as f:
                pickle.dump(g0, f, pickle.HIGHEST_PROTOCOL)

def heatmap(X, x_min, x_max, n_x, Y, y_min, y_max, n_y, nRep, t_sim = 50, filename = 'heatmap_files.file'):
    
    xs = np.linspace(x_min, x_max, num = n_x)
    ys = np.linspace(y_min, y_max, num = n_y)
    res = np.zeros((len(xs), len(ys)))
    cpt = 0
    for i, x in enumerate(xs):
        f = 0.5
        t0 = 0.1
        indel_size = 60
        nb_pol = None
        if X == 'f':
            f = x
        if X == 'T0':
            t0 = x
        if X == 'indel_size':
            indel_size = x
        if X == 'Nb_pol':
            nb_pol = x
        for j,y in enumerate(ys):
            mean_fitness = 0
            if Y == 'f':
                f = y
            if Y == 'T0':
                t0 = y
            if Y == 'indel_size':
                indel_size = y
            if Y == 'Nb_pol':
                nb_pol = y
            for r in range(nRep):
                print(cpt, '/', len(ys)*len(xs)*nRep ,', ',X,' : ', x, ', ', Y,' : ', y)
                g0 = Genome(pathToFiles=pathToFiles, nb_pol=nb_pol, f=f, T0=t0, indel_size=int(indel_size))
                g0.evolution(t_sim)
                mean_fitness += g0.fitness
                cpt+=1
            res[len(xs)-i-1,j] = 1.0*mean_fitness/nRep
    path = os.path.join(pathToFiles, 'Binary_files')
    with open(os.path.join(path, 'heatmap_files_'+X+'_'+Y+'.file'), "wb") as f:
        pickle.dump((res, xs, ys, X, Y), f, pickle.HIGHEST_PROTOCOL)
    return (res, xs, ys, X, Y)

#pathToFiles = 'D:/ProjetSimADN'
#pathToFiles = '/home/ocassan/ProjetSimADN'
#pathToFiles = '/home/amaury/ProjetSimADN'
#pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/'


pathToFiles = '/home/biosciences/users/OceaneAmauryJulie/ProjetSimADN'
#g0 = Genome(pathToFiles = pathToFiles, f = 0.5)
#g0.evolution(10, dump = True)


#nohup command &

heatmap('Nb_pol', 1, 10, 5, 'f', 0, 1, 4, nRep = 1, t_sim = 1)



