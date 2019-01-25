# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 16:31:36 2019

@author: Julie
"""

def encode_genome(genome):
    '''
    To visualize a genome
    Intergene region encoded by 0
    Barrier encoded by 1 
    Gene encoded by their height (their number)
    '''
    gen = list(genome)
    new_gen = []
    for index,value in enumerate(gen):
        if value == 'i':
            new_gen.append(0)
        elif value == 'b':
            new_gen.append(1)
        else:
            new_gen.append(int(value.split('g')[-1]))
    return new_gen

def visualize_genome_fitness(delta_fitness, stock_genomes, threshold):
    '''
    To visualize the genome structure
    If a delta of fitness goes above the threshold at step t
    This function will give the genome structure at time t-1 and t
    '''
    indices = []
    for i,v in enumerate(delta_fitness):
        if v > threshold:
            indices.append(i)
    
    for i in indices:
        plt.figure()
        new_genome1 = encode_genome(stock_genomes[i-1])
        new_genome2 = encode_genome(stock_genomes[i])
        
        plt.subplot(2,1,1)
        plt.plot(range(len(new_genome1)), new_genome1)
        title = "Genome at step " + str(i) + " and " + str(i+1)
        plt.title(title)
        plt.xlabel('Positions')
        
        plt.subplot(2, 1, 2)
        plt.plot(range(len(new_genome2)), new_genome2)
        plt.xlabel('Positions')

def visualize_nb_transcripts(genomes):
    '''
    Argument : genomes saved during an entire simulation given in argument
    Return different plots associated to these genomes :
        Number of transcripts per time
        Fitness per time
        Delta of fitness per time
        Structure of the genes in the genome when there is a big delta of fitness
    '''
    
    #Number of transcripts per gene as a function of time
    plt.figure()
    ax = plt.subplot(111)
    for gene in range(0,10):
        ax.plot(g.transcripts[gene], label = gene)
        ax.set_xlabel("Time t")
        ax.set_ylabel("Nombre de transcrits")
        plt.title("Nombre de transcrits par gène en fonction du temps")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1, chartBox.height])
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)

    #Fitness as a function of time
    plot_fitness(g)

    #Plot the genome sutructure depending on the threshold and the delta of fitness
    deltas = [g.fitnesses[i+1]-g.fitnesses[i-1] for i in range(0,len(g.fitnesses)-1)]
    plt.figure()
    #Define threshold
    sorted_deltas = sorted(deltas)
    thr = sorted_deltas[-9]
    plt.plot(deltas)
    print(thr, sorted_deltas)
    plt.hlines(thr, color = 'r', xmin = 0, xmax = len(deltas))
    visualize_genome_fitness(deltas, genomes, thr)

    plt.show()

#g = Genome(pathToFiles = pathToFiles, f = 0.8)
#genomes = g.evolution(10)
#visualize_nb_transcripts(genomes)
