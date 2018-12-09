# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 08:41:27 2018

@author: Océane
"""
import matplotlib.pyplot as plt
from collections import OrderedDict
import pickle
import numpy as np
import os

pathToFiles = 'D:/ProjetSimADN'
#pathToFiles = '/home/ocassan/ProjetSimADN'
#pathToFiles = '/home/amaury/ProjetSimADN'
#pathToFiles = '/home/julie/Documents/5BIM/BacteriaEvolution/ProjetSimADN/'

def plot_fitness(genome, fig_name = "fitness.png"): 
        #shows the fitness evolution in time
        #with the events as coloured circles on the curve
        plt.figure()
        plt.plot(genome.fitnesses, color = 'k')
        for ev in list(genome.events):
            if ev[1] == "insertion":
                col = 'green'
            elif ev[1] == "deletion":
                col = 'red'
            else:
                col = 'blue'
            plt.plot(ev[0],genome.fitnesses[ev[0]], 'o', color = col, label = ev[1])
        plt.title("Fitness of the genome in time")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        plt.ylabel('Fitness')
        plt.xlabel('Time')
        plt.savefig(fig_name)


def plot_hist(genome, event, fig_name = "hist.png"):
    #plots the histogram of the fitness change induced by the 
    #mutational events
    plt.figure()
    if event != 'inversion':
        plt.hist(genome.delta_fitness[event], bins = 20)
    else:
        plt.hist( [ev[0] for ev in genome.delta_fitness[event]], bins = 20)
    plt.title('Histogram of the fitness changes induced by '+event+'s')
    plt.xlabel('Fitness difference')
    plt.ylabel('Counts')
    
def plot_hists(genome, fig_name = "hist.png"):
    #plots the histogram of the fitness change induced by the 
    #mutational events
    plt.figure()
    for i, event in enumerate(['inversion', 'deletion', 'inversion']):
        fc=[0, 0, 0, 0.2]
        fc[i] = 1
        if event != 'inversion':
            plt.hist(genome.delta_fitness[event], bins = 20, fc = fc )
        else:
            plt.hist( [ev[0] for ev in genome.delta_fitness[event]], bins = 20, fc = fc)
    plt.title('Histogram of the fitness changes induced by '+event+'s')
    plt.xlabel('Fitness difference')
    plt.ylabel('Counts')
        
    
def plot_inv_size_fitness(genome):
    plt.figure()
    delta = genome.delta_fitness['inversion']
    plt.plot([d[1] for d in delta], [d[0] for d in delta], 'o')

def heatmap(X, Y):
    path = os.path.join(pathToFiles, 'Binary_files')
    with open(os.path.join(path, 'heatmap_files_'+X+'_'+Y+'.file'), "rb") as f:
        heatmap_files = pickle.load(f)
    res = heatmap_files[0]
    xs = heatmap_files[1]
    ys = heatmap_files[2]
    X = heatmap_files[3]
    Y = heatmap_files[4]
    
    print(res)
    
    fig, ax = plt.subplots()
    im = ax.imshow(res)
    #We want to show all ticks...
    ax.set_yticks(np.arange(len(xs)))
    ax.set_xticks(np.arange(len(ys)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(list(reversed([round(f, 1) for f in xs])))
    ax.set_ylabel(X)
    ax.set_xticklabels([round(f, 1) for f in ys])
    ax.set_xlabel(Y)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    ax.set_title("Heatmap of the last fitness value for different "+X+" and "+Y)
    fig.tight_layout()
    plt.savefig('heatmap'+X+'_'+Y+'.png')

path = os.path.join(pathToFiles, 'Binary_files')
with open(os.path.join(path, 'genome.file'), "rb") as f:
    genome = pickle.load(f)
#plot_hist(genome, 'insertion')

#heatmap()
plot_fitness(genome)
#plot_hist(genome, 'inversion')
#plot_hists(genome)
#plot_inv_size_fitness(genome)
