import os,sys
import configparser
import pandas as pd
import numpy as np
import collections as col
from pylab import *
import errno
import csv
from shutil import copy
from scipy.optimize import fsolve

###########################################################
#                       Functions                         #
###########################################################

def create_config_file(config, config_file_path, TSS_file, SIGMA_0, RNAPS_NB): # DELTA_X, D, J_0, SIGMA_0, RNAPS_NB,
    # Create the config file
    config.set('INPUTS','TSS', str(TSS_file))
    config.set('SIMULATION','SIGMA_0', str(SIGMA_0))
    config.set('SIMULATION','RNAPS_NB', str(RNAPS_NB))
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)

# Read the config files
def read_config_file(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    # Update INI file without removing comments
    config = configparser.ConfigParser(allow_no_value=True)
    if not os.path.exists(path):
        print("Input file was not found at path %s"%path)
        sys.exit(1)
    config.read(path)
    return config


"""
# Read the config files and return the values of each variable
# this function will be useful when we are in another script
def read_config_file_v2(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    config.read(path)
    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # get values from the config file
    m = config.getfloat('GLOBAL', 'm')
    sigma_t = config.getfloat('GLOBAL', 'sigma_t')
    epsilon = config.getfloat('GLOBAL', 'epsilon')

    RNAPs_genSC = config.getfloat('SIMULATION', 'RNAPs_genSC')
    DELTA_X = config.getfloat('SIMULATION', 'DELTA_X')
    DELTA_T = config.getfloat('SIMULATION', 'DELTA_T')
    RNAPS_NB = config.getint('SIMULATION', 'RNAPS_NB')
    SIM_TIME = config.getfloat('SIMULATION', 'SIM_TIME')
    OUTPUT_STEP = config.getfloat('SIMULATION', 'OUTPUT_STEP')

    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')
    TOPO_CTE = config.getfloat('SIMULATION', 'TOPO_CTE')
    GYRASE_CTE = config.getfloat('SIMULATION', 'GYRASE_CTE')
    #TOPO_EFFICIENCY = config.getfloat('SIMULATION', 'TOPO_EFFICIENCY')
    k_GYRASE = config.getfloat('SIMULATION', 'k_GYRASE')
    x0_GYRASE = config.getfloat('SIMULATION', 'x0_GYRASE')
    k_TOPO = config.getfloat('SIMULATION', 'k_TOPO')
    x0_TOPO = config.getfloat('SIMULATION', 'x0_TOPO')

    # Calculate SIGMA_0 based on Topoisomerases concentration.
    try:
        SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    except:
        # compute sigma_0 from concentrations of topoisomerases
        func = lambda sig0 : -GYRASE_CONC*1/(1+np.exp(-k_GYRASE*(sig0-x0_GYRASE)))*GYRASE_CTE + TOPO_CONC*1/(1+np.exp(k_TOPO*(sig0-x0_TOPO)))*TOPO_CTE
        sig0_initial_guess = -0.03
        SIGMA_0 = fsolve(func, sig0_initial_guess)[0]
        if not isinstance(SIGMA_0,float):
            print("Error computing SIGMA_0")

    return GFF_file, TSS_file, TTS_file, Prot_file, m, sigma_t, epsilon, SIGMA_0, DELTA_X, DELTA_T, RNAPS_NB, SIM_TIME, OUTPUT_STEP, GYRASE_CONC, TOPO_CONC, TOPO_CTE, GYRASE_CTE, k_GYRASE, x0_GYRASE, k_TOPO, x0_TOPO
"""






###################### Reading files ######################

# you can combine those two functions
def load_gff(filename):
    gff_df_raw = pd.read_table(filename, sep='\t', comment='#', header=0)
    return gff_df_raw

def load_tab_file(filename):
    data = pd.read_table(filename, sep='\t', header=0)
    return data

def str2num(s):
    s[s == '+'] = 1 #True
    s[s == '-'] = -1 #False
    return s

def get_tr_nbr_csv(csv_file):
    csv_tr_nbr = pd.read_csv(csv_file, sep=';', header=None)
    tr_nbr = csv_tr_nbr.values
    return tr_nbr.flatten()

######################## Others ###########################

# Get the genome size from the header of gff file (befor renaming it)
def get_genome_size(gff_df):
    genome_size = int(gff_df.columns[4]) - int(gff_df.columns[3])
    return genome_size

# Rename the header (columns names)
def rename_gff_cols(gff_df):
    names=["seqid", "source", "type","start","end","score","strand","phase","attributes"]
    gff_df.columns = names
    return gff_df

# # Whether the gene is on the + strand or - strand
# def in_forward(tr_id):
#     if strands[tr_id] == 1. :
#         return True
#     else:
#         return False


# # Get the transciption unit with the list of transcripts
# def get_TU(TUindex_list):
#     TU_dict = col.defaultdict(list)
#     for index, TUindex in enumerate(TUindex_list):
#         TU_dict[TUindex].append(index)
#     return TU_dict

# calculate the initiation rate
def f_init_rate(tr_prob, sig, sigma_t, epsilon, m):
    tr_prob_sig = tr_prob * np.exp((1/(1+np.exp((sig-sigma_t)/epsilon)))*m)
    return tr_prob_sig



# Get the list of all possible transcripts
def get_tr_info(tss, tts, TU_tts, Kon, Poff):
    this_TU_tts = []
    tr_id = []
    tr_start = []
    tr_end = []
    tr_strand = []
    tr_size = []
    tr_rate = []
    sum_Kon = np.sum(Kon)

    j = 0 # trancript id indice
    for i in tss.index.values: # All TSSs
        # get the TU of this tss
        TU_id = tss['TUindex'][i]
        # the list of TTS that are in the same TU of this tss_id (i)
        # TU_tts ex : defaultdict(list, {0: [1150, 2350], 1: [6250]})
        # On prend tt les tts qui existent dans la meme UT du tss choisi
        this_TU_tts = TU_tts[TU_id] # pour 0 => [1150, 2350]
        # + or -
        if tss['TUorient'][i] == '+' :
            # go right
            k = 0 # TTS id index : k start from the first position of each TU
            proba_rest = 1
            while proba_rest > 0 :
                if tss['TSS_pos'][i] < this_TU_tts[k]:
                    tr_id.append(j)
                    tr_strand.append(1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting the TSSs, we shall (in every loop) generate a new tr_end
                    tr_end.append(this_TU_tts[k])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
        else:
            # go left
            k = 0
            proba_rest = 1
            while proba_rest > 0 and k < len(this_TU_tts) :
                if this_TU_tts[k] < tss['TSS_pos'][i] :
                    tr_id.append(j)
                    tr_strand.append(-1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting them, we shall (in every loop) generate a new tr_end
                    tr_end.append(this_TU_tts[i])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
    tr_size = np.abs(np.array(tr_start) - np.array(tr_end))
    ts_beg_all_trs = np.zeros(len(tr_id), dtype=int64)
    ts_remain_all = np.around(tr_size)
    return (tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all)








# Get the list of all possible transcripts
def get_tr_info_1(tss, tts, TU_tts, Kon, Poff):
    this_TU_tts = []
    tr_id = []
    tr_start = []
    tr_end = []
    tr_strand = []
    tr_size = []
    tr_rate = []
    sum_Kon = np.sum(Kon)

    j = 0 # trancript id indice
    for i in tss.index.values: # All TSSs
        # get the TU of this tss
        TU_id = tss['TUindex'][i]
        # the list of TTS that are in the same TU of this tss_id (i)
        # TU_tts ex : defaultdict(list, {0: [1150, 2350], 1: [6250]})
        # On prend tt les tts qui existent dans la meme UT du tss choisi
        this_TU_tts = TU_tts[TU_id] # pour 0 => [1150, 2350]
        # + or -
        if tss['TUorient'][i] == '+' :
            # go right
            k = TU_id # TTS id index : k start from the first position of each TU
            proba_rest = 1
            while proba_rest > 0 :
                if tss['TSS_pos'][i] < tts['TTS_pos'][k]:
                    tr_id.append(j)
                    tr_strand.append(1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting the TSSs, we shall (in every loop) generate a new tr_end
                    tr_end.append(tts['TTS_pos'][i])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
        else:
            # go left
            k = 0
            proba_rest = 1
            while proba_rest > 0 and k < len(this_TU_tts) :
                if tts['TTS_pos'][k] < tss['TSS_pos'][i] :
                    tr_id.append(j)
                    tr_strand.append(-1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting them, we shall (in every loop) generate a new tr_end
                    tr_end.append(tts["TTS_pos"][i])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
    tr_size = np.abs(np.array(tr_start) - np.array(tr_end))
    ts_beg_all_trs = np.zeros(len(tr_id), dtype=int64)
    ts_remain_all = np.around(tr_size)
    return (tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all)






# Get the list of all possible transcripts
def get_tr_info_old(tss, tts, TU_tts, Kon, Poff):
    this_TU_tts = []
    tr_id = []
    tr_start = []
    tr_end = []
    tr_strand = []
    tr_size = []
    tr_rate = []
    sum_Kon = np.sum(Kon)

    j = 0 # trancript id indice
    for i in tss.index.values: # All TSSs
        # get the TU of this tss
        TU_id = tss['TUindex'][i]
        # the list of TTS that are in the same TU of this tss_id (i)
        # TU_tts ex : defaultdict(list, {0: [1150, 2350], 1: [6250]})
        # On prend tt les tts qui existent dans la meme UT du tss choisi
        this_TU_tts = TU_tts[TU_id] # pour 0 => [1150, 2350]
        # + or -
        if tss['TUorient'][i] == '+' :
            # go right
            k = TU_id # TTS id index : k start from the first position of each TU
            proba_rest = 1
            while proba_rest > 0 :
                if tss['TSS_pos'][i] < tts['TTS_pos'][k]:
                    tr_id.append(j)
                    tr_strand.append(1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting the TSSs, we shall (in every loop) generate a new tr_end
                    tr_end.append(tts['TTS_pos'][k])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
        else:
            # go left
            k = 0
            proba_rest = 1
            while proba_rest > 0 and k < len(this_TU_tts) :
                if tts['TTS_pos'][k] < tss['TSS_pos'][i] :
                    tr_id.append(j)
                    tr_strand.append(-1)
                    tr_start.append(tss['TSS_pos'][i])
                    # after getting them, we shall (in every loop) generate a new tr_end
                    tr_end.append(this_TU_tts[k])
                    # the probability to choose a specific transcript
                    tr_rate.append(Kon[i] * (Poff[k] * proba_rest))
                    proba_rest = (1 - Poff[k]) * proba_rest
                    j += 1
                k += 1
    tr_size = np.abs(np.array(tr_start) - np.array(tr_end))
    ts_beg_all_trs = np.zeros(len(tr_id), dtype=int64)
    ts_remain_all = np.around(tr_size)
    return (tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all)

def f_prob_init_rate(init_rate, sum_init_rate, DELTA_T):
    return (1-np.exp(-sum_init_rate*DELTA_T)) * (init_rate/sum_init_rate)

def f_prob_unhooked_rate(sum_Kon, DELTA_T, RNAPs_unhooked_nbr):
    return np.exp(-sum_Kon*DELTA_T)/RNAPs_unhooked_nbr

# Get the transciption unit with the list of tts belonging to TU.
def get_TU_tts(tss, tts):
    TU_tts = col.defaultdict(list)
    for index, TUindex in enumerate(tss['TUindex'].values):
        TU_tts[TUindex].append(tts['TTS_pos'][index])
    return TU_tts

def calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO, TOPO_CTE, DELTA_T):
    d_sigma = (-GYRASE_CONC*1/(1+np.exp(-k_GYRASE*(Barr_sigma-x0_GYRASE)))*GYRASE_CTE + TOPO_CONC*1/(1+np.exp(k_TOPO*(Barr_sigma-x0_TOPO)))*TOPO_CTE) * DELTA_T
    Barr_sigma += d_sigma

    return Barr_sigma

###################### Saving files #######################

def save_files(output_path,
                Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma,
                tr_nbr, tr_times, save_RNAPs_info, save_tr_info,
                save_Dom_sigma, save_Barr_pos, save_mean_sig_wholeGenome, save_Dom_size,
                DELTA_X, RNAPs_genSC,
                RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id, RNAPs_hooked_id,
                RNAPs_strand, ts_beg, ts_remain, save_nbr_RNAPs_hooked,
                init_rate, Kon, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC):

    # make sure that the output direcory exists, and create one if it doesn't
    os.makedirs("%s/resume_sim" %output_path, exist_ok=True)
    os.makedirs("%s/all_res" %output_path, exist_ok=True)

    # save tr_nbr
    tr_nbr = pd.Series(tr_nbr)
    tr_nbr.to_csv("%s/save_tr_nbr.csv" %output_path, sep=';', index=False)

    # convert tr_times dict to pandas serie
    tr_times = pd.DataFrame.from_dict(tr_times, orient='index')
    # save the tr_times to csv file
    tr_times.to_csv("%s/save_tr_times.csv" %output_path, sep=';', index=True, header = False)
    # save the number of RNAPs hooked
    np.savez("%s/save_nbr_RNAPs_hooked.npz" %output_path, nbr_RNAPs_hooked = save_nbr_RNAPs_hooked)

    # Save last info
    np.savez("%s/resume_sim/resume_sim_RNAPs.npz" %output_path, RNAPs_tr = RNAPs_tr,
                                                               RNAPs_pos = RNAPs_pos,
                                                               RNAPs_unhooked_id = RNAPs_unhooked_id,
                                                               RNAPs_strand = RNAPs_strand,
                                                               ts_beg = ts_beg,
                                                               ts_remain = ts_remain,
                                                               RNAPs_hooked_id = RNAPs_hooked_id)

    np.savez("%s/resume_sim/resume_sim_tr.npz" %output_path, tr_nbr = tr_nbr,
                                                            init_rate = init_rate)

    np.savez("%s/resume_sim/resume_sim_Barr.npz" %output_path, Barr_pos = Barr_pos,
                                                              Barr_type = Barr_type,
                                                              Dom_size = Dom_size,
                                                              Barr_ts_remain = Barr_ts_remain,
                                                              Barr_sigma = Barr_sigma)

    # Save all info
    np.savez("%s/all_res/save_RNAPs_info" %output_path, RNAPs_info = save_RNAPs_info)
    np.savez("%s/all_res/save_tr_info" %output_path, tr_info = save_tr_info)
    np.savez("%s/all_res/save_sigma_info" %output_path, dom_sigma_info = save_Dom_sigma, save_Barr_pos = save_Barr_pos, mean_sig_wholeGenome = save_mean_sig_wholeGenome, Dom_size = save_Dom_size)

###########################################################
#         Transcription Process (Simulation)              #
###########################################################

def start_transcribing(INI_file, first_output_path=None, resume_output_path=None, resume=False):

    """Example function with types documented in the docstring.

    `PEP 484`_ type annotations are supported. If attribute, parameter, and
    return types are annotated according to `PEP 484`_, they do not need to be
    included in the docstring:

    Args:
        INI_file (str): The path to the parameter file.
        first_output_path (str, optional): The path in which all the generated files (when) will be saved.
        resume_output_path (str, optional): The path in which all the generated files (when resuming) will be saved.
        resume (bool, optional): Whether this function is used to resume an already done simulation or not.

    Returns:
        GFF_file : Relative path to the GFF file
        TSS_file : Relative path to the TSS file
        TTS_file : Relative path to the TTS file
        SIM_TIME : Simulation time
        RNAPS_NB : Numbers of RNAPol
        tr_nbr : Number of transcripts
        tr_times : Transcription time
        init_rate : Initiation rate
        RNAPs_tr : Contains the id of the picked transcript
        RNAPs_pos : The position of RNAPols
        RNAPs_unhooked_id : id of Unhoocked RNAPols
        save_RNAPs_info : Contains the RNAPs_tr and RNAPs_pos
        save_tr_info : Contains the tr_nbr and init_rate
        save_Dom_sigma : Contains the SC density (sigma) value of each domaine
        save_Barr_pos : Contains the barriers positions (whether Barr_fix or RNAPol)
        cov_bp : Coverage
        tr_end : The end (position) of transcripts

    """

    ###########################################################
    #                 initiation of variables                 #
    ###########################################################

    ####################### Params info ###################

    config = read_config_file(INI_file)

    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # get values from the config file
    m = config.getfloat('PROMOTER', 'm')
    sigma_t = config.getfloat('PROMOTER', 'sigma_t')
    epsilon = config.getfloat('PROMOTER', 'epsilon')

    DELTA_X = config.getfloat('GLOBAL', 'DELTA_X')
    DELTA_T = config.getfloat('GLOBAL', 'DELTA_T')

    RNAPs_genSC = config.getfloat('SIMULATION', 'RNAPs_genSC')
    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    RNAPS_NB = config.getint('SIMULATION', 'RNAPS_NB')
    SIM_TIME = config.getfloat('SIMULATION', 'SIM_TIME')
    OUTPUT_STEP = config.getfloat('SIMULATION', 'OUTPUT_STEP')
    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')

    TOPO_CTE = config.getfloat('TOPOISOMERASES', 'TOPO_CTE')
    GYRASE_CTE = config.getfloat('TOPOISOMERASES', 'GYRASE_CTE')
    k_GYRASE = config.getfloat('TOPOISOMERASES', 'k_GYRASE')
    x0_GYRASE = config.getfloat('TOPOISOMERASES', 'x0_GYRASE')
    k_TOPO = config.getfloat('TOPOISOMERASES', 'k_TOPO')
    x0_TOPO = config.getfloat('TOPOISOMERASES', 'x0_TOPO')
    # Calculate SIGMA_0 based on Topoisomerases concentration.
    #SIGMA_0 = 0 #((-np.log(((GYRASE_CONC*GYRASE_CTE)/TOPO_CONC*TOPO_CTE)-1))/k)+x_0

    # path to the input files (remove the "params.ini" from the path)
    pth = INI_file.rpartition("/")[0] #+ "/"
    # or you can use : pth = os.path.split(INI_file)[0] + "/"
    #if pth=="/":
        #pth="./"
    if pth=="":
        pth="."
    if pth[-1]!="/":
        pth+="/"

    gff_df_raw = load_gff(pth+GFF_file)
    tss = load_tab_file(pth+TSS_file)
    tts = load_tab_file(pth+TTS_file)
    # we will load the prot_file later

    # get the TSS position
    TSS_pos = (tss['TSS_pos'].values/DELTA_X).astype(int)

    # get the initiation rate (Kon)
    Kon = tss['TSS_strength'].values

    # get the Poff
    Poff = tts['TTS_proba_off'].values

    # get the genome size
    genome_size = get_genome_size(gff_df_raw)
    gff_df = rename_gff_cols(gff_df_raw)

    # Dict of transciption units with the list of tts belonging to TU.
    # One TU starts from a single TSS but can have several TTS...
    TU_tts = get_TU_tts(tss, tts)

    # The RNAPs id
    RNAPs_id = np.full(RNAPS_NB, range(0, RNAPS_NB), dtype=int)

    # RNAPs_last_pos
    RNAPs_last_pos = np.full(RNAPS_NB, NaN)

    ## get the strands orientation
    # strands = str2num(gff_df['strand'].values)

    # list of all possible transcripts
    tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all = get_tr_info_1(tss, tts, TU_tts, Kon, Poff)
    print(tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all)
    # convert all variables to numpy array
    tr_id = np.array(tr_id)
    tr_strand = np.array(tr_strand)

    tr_start = np.array(tr_start)/DELTA_X
    tr_start = tr_start.astype(int64)

    tr_end = np.array(tr_end)/DELTA_X
    tr_end = tr_end.astype(int64)

    tr_rate = np.array(tr_rate)

    tr_size = np.array(tr_size)/DELTA_X
    tr_size = tr_size.astype(int64)

    ts_beg_all_trs = np.array(ts_beg_all_trs)

    ts_remain_all = np.array(ts_remain_all)/DELTA_X
    ts_remain_all = ts_remain_all.astype(int64)

    genome = int(genome_size/DELTA_X)

    if resume == False:
        # The position of RNAPs
        RNAPs_pos = np.full(RNAPS_NB, NaN)

        # The number of times transcripts has been transcribed
        tr_nbr = np.zeros(len(tr_id), dtype=int)


        # try to read the prot_file which contains the fixed barriers positions
        try:
            prot = load_tab_file(pth+Prot_file)
            Barr_fix = (prot['prot_pos'].values/DELTA_X).astype(int)

            # just for the echo we can assign it directely
            Barr_pos = np.copy(Barr_fix)

            # update in case where no fixed barriers !!!
            # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
            Dom_size = np.abs(np.ediff1d(Barr_pos)) #, dtype=int

            Dom_size = np.append(Dom_size, genome-(np.max(Barr_pos)-np.min(Barr_pos)))

            # Barr_type contains the barrier type
            # 0 : fixed barrrier (e.g Protein)
            # -1 : -RNAPol (direction : <--)
            # 1 : +RNAPol (direction : -->)
            Barr_type = np.full(len(Barr_fix), 0, dtype=int)
            Barr_sigma = np.full(len(Barr_fix), SIGMA_0)
            # here we need to make an Barr_ts_remain
            # to track the position of each RNAPol
            # each position in Barr_ts_remain is associated with the same position in Barr_pos
            Barr_ts_remain = np.full(len(Barr_fix), NaN) # The Barr_ts_remain of fixed barr is NaN

        # if prot_file is empty or doesn't exist then:
        except (pd.io.common.EmptyDataError, OSError, ValueError):
            # we'll have one Dom_size which is the whole genome
            # There is no Barr_fix
            Dom_size = np.array([genome], dtype=int)
            # create the other variables
            Barr_type = np.array([], dtype=int)
            Barr_sigma = np.array([SIGMA_0]) # on the whole genome
            Barr_pos = np.array([], dtype=int)
            Barr_ts_remain = np.array([])


        RNAPs_unhooked_id = np.copy(RNAPs_id)
        RNAPs_strand = np.full(RNAPS_NB, NaN)
        ts_beg = np.full(RNAPS_NB, NaN)
        ts_remain = np.full(RNAPS_NB, NaN)
        # RNAPs_tr will contain the id of the picked transcript
        RNAPs_tr = np.full(RNAPS_NB, -1, dtype=(int64))
        # get the TSSs ids
        tss_id = tss.index.values

        # in the case of RNAP_NBR = 0
        RNAPs_hooked_id = []

        # no resume ==> set the 'output_path' to 'first_output_path'
        output_path = first_output_path
        # 'output_path' variable is used to save the output files
        # in 'save_files()' function

    else:
        #RNAPs_info contains : ['RNAPs_unhooked_id', 'RNAPs_pos', 'RNAPs_tr']
        RNAPs_info = np.load("%s/resume_sim/resume_sim_RNAPs.npz" %first_output_path)

        # get the RNAPs position
        RNAPs_pos = RNAPs_info['RNAPs_pos']

        # The number of times transcripts has been transcribed
        csv_path = "%s/save_tr_nbr.csv" %first_output_path
        tr_nbr = get_tr_nbr_csv(csv_path)

        # when we resume, we won't need info from 'prot' file because we already
        # have all what we need in 'resume_sim_Barr.npz' file ;)

        # Get info from NPZ file
        Barr_info = np.load("%s/resume_sim/resume_sim_Barr.npz" %first_output_path)

        # aaand here we go !
        Barr_pos = Barr_info['Barr_pos']
        Dom_size = Barr_info['Dom_size']

        Barr_type = Barr_info['Barr_type']
        Barr_sigma = Barr_info['Barr_sigma']

        # here we need to make an Barr_ts_remain
        # so we can track the position of each RNAPol
        # each position in Barr_ts_remain is associated with the same position in Barr_pos
        Barr_ts_remain = Barr_info['Barr_ts_remain']

        ## do the same for RNAPs_info
        # get the RNAPs_info
        RNAPs_info = np.load("%s/resume_sim/resume_sim_RNAPs.npz" %first_output_path)

        # get the RNAPs_hooked_id and RNAPs_pos
        RNAPs_hooked_id = RNAPs_info["RNAPs_hooked_id"]
        RNAPs_pos = RNAPs_info["RNAPs_pos"]
        # deduce the RNAPs_hooked_pos from the extracted info ;)
        #!!! check out this one (it's not used)
        RNAPs_hooked_pos = RNAPs_pos[RNAPs_hooked_id].astype(int)

        # since we continue the simulation, we shall retrieve the RNAPs_unhooked_id from the npz file
        RNAPs_unhooked_id = RNAPs_info["RNAPs_unhooked_id"]

        RNAPs_strand = RNAPs_info["RNAPs_strand"]
        ts_beg = RNAPs_info["ts_beg"]
        ts_remain = RNAPs_info["ts_remain"]
        # RNAPs_tr contains the id of the picked transcript
        RNAPs_tr = RNAPs_info["RNAPs_tr"]
        # get the TSSs ids
        tss_id = tss.index.values

        # will do the same for RNAPs_hooked_id
        RNAPs_hooked_id = RNAPs_info["RNAPs_hooked_id"]

        # resume ==> set the 'output_path' to 'resume_output_path'
        output_path = resume_output_path
        # 'output_path' variable is used to save the output files
        # in 'save_files()' function

    ######### Variables used to get the coverage ##########

    id_shift_fwd = list(range(1, genome))
    id_shift_fwd.append(0)
    id_shift_fwd = np.array(id_shift_fwd)
    id_shift_bwd = list(range(0, genome-1))
    id_shift_bwd.insert(0, genome-1)
    id_shift_bwd = np.array(id_shift_bwd)

    cov_bp = np.arange(0, genome_size, DELTA_X)
    cov_bp = np.resize(cov_bp, genome)

    # save the time when RNApoly is starting trasncribing a specific transcript
    #tr_times = col.defaultdict(list)
    tr_times={}
    for transcript in tr_id:
        tr_times[transcript]=[]

    # numpy array where all RNAPs info will be saved except the nbr_RNAPs_hooked
    save_RNAPs_info = np.full([RNAPS_NB, 2, int(SIM_TIME/(DELTA_T*OUTPUT_STEP))], np.nan) # nbr d'ele (cols)

    # this array will contain the number of RNAPs hooked at each time step
    # the length of the array is equivalent to the simulation length
    save_nbr_RNAPs_hooked = np.full(int(SIM_TIME/(DELTA_T*OUTPUT_STEP)), np.nan)

    # the same for transcripts info
    save_tr_info = np.full([len(tr_id), 2, int(SIM_TIME/(DELTA_T*OUTPUT_STEP))], np.nan)

    # # in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    save_Dom_sigma = list()




    # in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    save_Dom_sigma = list()
    save_Dom_size = list()
    save_Barr_pos = list()
    save_mean_sig_wholeGenome = list()
    save_Barr_pos = list()

    ########### Go !

    for t in range(0,int(SIM_TIME/DELTA_T)):
        # we need to know each TSS belong to which Domaine
        # if we have only one domaine (e.g no Barr_fix) then
        # all the TSS belong to the only domaine that exists (e.g TSS_pos_idx will be 0)
        TSS_pos_idx = np.searchsorted(Barr_pos, TSS_pos)

        # after knowing the domaine of each TSS we can get sigma
        try:
            sigma_tr_start = Barr_sigma[TSS_pos_idx-1]
        except IndexError:
            sigma_tr_start = np.array([SIGMA_0])

        # get the initiation rates
        # calculate the initiation rate of each transcript/gene
        init_rate = f_init_rate(tr_rate, sigma_tr_start, sigma_t, epsilon, m)
        # use the calculated init_rate to get the probability
        # of RNAPol's binding to each TSS (e.g prob_init_rate)
        sum_init_rate = np.sum(init_rate)
        prob_init_rate = f_prob_init_rate(init_rate, sum_init_rate, DELTA_T)

        if np.size(RNAPs_unhooked_id)!=0:
            # get the probability of RNAPol's stay unhooked
            prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, DELTA_T, len(RNAPs_unhooked_id))
            # craete the numpy array
            prob_unhooked_rate = np.full(len(RNAPs_unhooked_id), prob_unhooked_rate)
            all_prob = np.concatenate([prob_init_rate, prob_unhooked_rate])
            # create the numpy array that will contains [ nTSS , Unhooked RNAPS ]
            # e.g if we have 2 TSSs and 3 unhooked RNAPols then
            # tss_and_unhooked_RNAPs = [0, 1, -1, -1, -1]
            tss_and_unhooked_RNAPs = np.concatenate([tss_id, np.full(len(RNAPs_unhooked_id), -1, dtype=int)])
            # pick up a random transcipt
            picked_tr = np.random.choice(tss_and_unhooked_RNAPs, len(RNAPs_unhooked_id), replace=False, p=all_prob) #RNAPs_unhooked_id
            # This is the KEY !
            picked_tr_hooked_id = picked_tr[np.where(picked_tr!=-1)[0]]
            picked_tr_unhooked_id = picked_tr[np.where(picked_tr==-1)[0]]

            new_RNAPs_hooked_id = RNAPs_unhooked_id[np.where(picked_tr!=-1)[0]]

            RNAPs_tr[new_RNAPs_hooked_id] = picked_tr[picked_tr!=-1]
            RNAPs_strand[new_RNAPs_hooked_id] = tr_strand[picked_tr[np.where(picked_tr!=-1)]]

            # The new position of each polymerase
            # if there is no RNAP already at this position
            RNAPs_pos[new_RNAPs_hooked_id] = tr_start[picked_tr[np.where(picked_tr!=-1)]].astype(int)

            # Bug #1 Fix
            # We use new_RNAPs_hooked_id and RNAPs_pos[new_RNAPs_hooked_id] to create an ordered dictionary
            # new_hooked_RNAPs_pos_dict contains [new hooked RNAP position as a key: new RNAP hooked id as value]
            new_hooked_RNAPs_pos_dict = dict(zip(RNAPs_pos[new_RNAPs_hooked_id], new_RNAPs_hooked_id))
            # We sort the dict by Keys (e.g position)
            # At the end we'll get the sorted positions with their corresponding ids
            new_hooked_RNAPs_pos_ordered_dict = col.OrderedDict(sorted(new_hooked_RNAPs_pos_dict.items()))
            # Here idx are sorted based on the positions
            new_hooked_RNAPs_idx_sorted = [idx for idx in new_hooked_RNAPs_pos_ordered_dict.values()]
            new_hooked_RNAPs_pos_sorted = [pos for pos in new_hooked_RNAPs_pos_ordered_dict.keys()]

            # take the positions and use them to get the index in which we will insert the new recruited RNAPs in Barr_pos array
            Barr_pos_RNAPs_idx = np.searchsorted(Barr_pos, new_hooked_RNAPs_pos_sorted)
            # after getting the idx, we start inserting
            # Barr_pos = np.insert(in Barr pos, in the follwing positions, the new hooked RNAPs)
            Barr_pos = np.insert(Barr_pos, Barr_pos_RNAPs_idx, new_hooked_RNAPs_pos_sorted)
            # if we have at least two barrier
            try:
                # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
                Dom_size = np.abs(np.ediff1d(Barr_pos))
                Dom_size = np.append(Dom_size, genome-(np.max(Barr_pos)-np.min(Barr_pos)))
                # if at least two RNAPols are hooked at the same time
                # and we have at least one Barr_pos (that already exists)
                # then we insert a new Barr_sigma value
                # otherwise we won't insert anything cuz we still have one domaine
                if len(new_hooked_RNAPs_idx_sorted) >= 1 and len(Barr_pos) > 1:
                    Barr_sigma = np.insert(Barr_sigma, Barr_pos_RNAPs_idx, Barr_sigma[Barr_pos_RNAPs_idx-1])
                    # if no RNAPol is hoooked
                    # and at least two RNAPols are hooked at the same time we remove SIGMA_0
                    # since we have already Barr_sigma=([SIGMA_0])
                    if len(Dom_size) < len(Barr_sigma):
                        # remove Barr_sigma=([SIGMA_0]) which still there
                        Barr_sigma = np.delete(Barr_sigma, -1)
            # in case we have one or zero barrier
            except Exception:
                # in case we have less than 2 Barriers
                Dom_size = np.array([genome])
                Barr_sigma = np.array([SIGMA_0])

            Barr_type = np.insert(Barr_type, Barr_pos_RNAPs_idx, RNAPs_strand[new_hooked_RNAPs_idx_sorted])

            # RNAPs_last_pos
            RNAPs_last_pos[new_hooked_RNAPs_idx_sorted] = tr_end[picked_tr_hooked_id]
            ts_beg[new_hooked_RNAPs_idx_sorted] = 0
            ts_remain[new_hooked_RNAPs_idx_sorted] = ts_remain_all[picked_tr_hooked_id] # NOT picked_tr
            Barr_ts_remain = np.insert(Barr_ts_remain, Barr_pos_RNAPs_idx, ts_remain[new_hooked_RNAPs_idx_sorted])
            RNAPs_hooked_id = np.where(RNAPs_tr!=-1)[0]

        ts_beg[RNAPs_hooked_id] += 1
        ts_remain[RNAPs_hooked_id] -= 1

        # save the time when RNAPol FINISHS trasncribing a specific transcript
        for x in RNAPs_tr[np.where(ts_remain==0)] :
            tr_times[x].append(t*DELTA_T) # + 0.5

        tr_nbr[RNAPs_tr[np.where(ts_remain==0)]]+=1

        # look in the net : numpy where two conditions
        Barr_ts_remain[np.where(Barr_type == -1)]-=1
        Barr_ts_remain[np.where(Barr_type == 1)]-=1

        # Get the index of RNAPs to remove
        rm_RNAPs_idx = np.where(Barr_ts_remain == 0)[0]

        # recover sigma value of the removed position
        removed_sigma = Barr_sigma[rm_RNAPs_idx]
        removed_dom_size = Dom_size[rm_RNAPs_idx]

        # recover the old_dom_size : the size of the previous domaine before combination/merging
        old_dom_size = Dom_size[rm_RNAPs_idx-1]
        old_sigma = Barr_sigma[rm_RNAPs_idx-1]

        # update Dom_size
        Dom_size[rm_RNAPs_idx-1] += removed_dom_size
        # or
        # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
        #Dom_size = np.abs(np.ediff1d(Barr_pos))
        #Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0])

        Barr_sigma[rm_RNAPs_idx-1] = (old_dom_size*old_sigma+removed_dom_size*removed_sigma)/(old_dom_size+removed_dom_size)

        # and reomve them
        Barr_pos = np.delete(Barr_pos, rm_RNAPs_idx)
        Barr_type = np.delete(Barr_type, rm_RNAPs_idx)
        Barr_ts_remain = np.delete(Barr_ts_remain, rm_RNAPs_idx)
        Barr_sigma = np.delete(Barr_sigma, rm_RNAPs_idx)
        Dom_size = np.delete(Dom_size, rm_RNAPs_idx)

        # update the RNAPs_tr array
        RNAPs_tr[np.where(ts_remain==0)] = -1
        # update the RNAPs_unhooked_id based on RNAPs_tr
        RNAPs_unhooked_id = np.where(RNAPs_tr==-1)[0]

        # reset the arrays
        RNAPs_strand[RNAPs_unhooked_id] = NaN
        RNAPs_pos[RNAPs_unhooked_id] = NaN
        RNAPs_last_pos[RNAPs_unhooked_id] = NaN
        ts_beg[RNAPs_unhooked_id] = NaN
        ts_remain[RNAPs_unhooked_id] = NaN

        Barr_pos[np.where(Barr_type == -1)]-=1
        Barr_pos[np.where(Barr_type == 1)]+=1

        # Update the position of polymerases still transcribing
        RNAPs_pos[np.where(RNAPs_strand == 1)]+=1
        RNAPs_pos[np.where(RNAPs_strand == -1)]-=1

        # Update the Dom_size (+1 or -1)
        # if we have at least two barrier
        try:
            # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
            Dom_size = np.abs(np.ediff1d(Barr_pos))
            Dom_size = np.append(Dom_size, genome-(np.max(Barr_pos)-np.min(Barr_pos)))
        # in case we have one or zero barrier
        except (IndexError, ValueError):
            Dom_size = np.array([genome])
            Barr_sigma = np.array([SIGMA_0])

        # UPDATE SIGMA
        # R_plus_pos : the ids of RNA pol in the + strand
        R_plus_pos = np.where(Barr_type == 1)[0].astype(int)
        # R_minus_pos : the ids of RNA pol in the - strand
        R_minus_pos = np.where(Barr_type == -1)[0].astype(int)
        #### Extract all types of domaines (Those are ids of domaines)
        # Barr_type_ahead to make the extraction circular ;)
        Barr_type_ahead = np.roll(Barr_type, -1)
        # __|__________O+____
        Barr_Dom_RPlus = np.where((Barr_type==0) & (Barr_type_ahead==1))
        # __|__________O-____
        Barr_Dom_RMinus = np.where((Barr_type==0) & (Barr_type_ahead==-1))
        # __|__________|_____
        Barr_Dom_Barr = np.where((Barr_type==0) & (Barr_type_ahead==0))
        # ___O+_________O+___
        RPlus_Dom_RPlus = np.where((Barr_type==1) & (Barr_type_ahead==1))
        # ___O-_________O-___
        RMinus_Dom_RMinus = np.where((Barr_type==-1) & (Barr_type_ahead==-1))
        # ___O+_________O-___
        RPlus_Dom_RMinus = np.where((Barr_type==1) & (Barr_type_ahead==-1))
        # ___O-_________O+___
        RMinus_Dom_RPlus = np.where((Barr_type==-1) & (Barr_type_ahead==+1))
        # ___O-_________|____
        RMinus_Dom_Barr = np.where((Barr_type==-1) & (Barr_type_ahead==0))
        # ___O+_________|____
        RPlus_Dom_Barr = np.where((Barr_type_ahead==0) & (Barr_type==+1))


        #### And then correct the value of Sigma in each case (before/after)
        corr_sig_Barr_Dom_RPlus = (Dom_size[Barr_Dom_RPlus]-1)/(Dom_size[Barr_Dom_RPlus]) # Sigma decrease x1
        corr_sig_Barr_Dom_RMinus = (Dom_size[Barr_Dom_RMinus]+1)/(Dom_size[Barr_Dom_RMinus]) # Sigma increase x1
        corr_sig_Barr_Dom_Barr = (Dom_size[Barr_Dom_Barr])/(Dom_size[Barr_Dom_Barr]) # Sigma FIX
        corr_sig_RPlus_Dom_RPlus = (Dom_size[RPlus_Dom_RPlus])/(Dom_size[RPlus_Dom_RPlus]) # Sigma FIX
        corr_sig_RMinus_Dom_RMinus = (Dom_size[RMinus_Dom_RMinus])/(Dom_size[RMinus_Dom_RMinus]) # Sigma FIX
        corr_sig_RPlus_Dom_RMinus = (Dom_size[RPlus_Dom_RMinus]+2)/(Dom_size[RPlus_Dom_RMinus]) # Sigma increase x2
        corr_sig_RMinus_Dom_RPlus = (Dom_size[RMinus_Dom_RPlus]-2)/(Dom_size[RMinus_Dom_RPlus]) # Sigma decrease x2
        corr_sig_RMinus_Dom_Barr = (Dom_size[RMinus_Dom_Barr]-1)/(Dom_size[RMinus_Dom_Barr]) # Sigma decrease x1
        corr_sig_RPlus_Dom_Barr = (Dom_size[RPlus_Dom_Barr]+1)/(Dom_size[RPlus_Dom_Barr]) # Sigma increase x1

        ### Multiply Sigma *= Corr (Each sigma value correspond to an specific domaine)
        Barr_sigma[Barr_Dom_RPlus] *= corr_sig_Barr_Dom_RPlus
        Barr_sigma[Barr_Dom_RMinus] *= corr_sig_Barr_Dom_RMinus
        Barr_sigma[Barr_Dom_Barr] *= corr_sig_Barr_Dom_Barr
        Barr_sigma[RPlus_Dom_RPlus] *= corr_sig_RPlus_Dom_RPlus
        Barr_sigma[RMinus_Dom_RMinus] *= corr_sig_RMinus_Dom_RMinus
        Barr_sigma[RPlus_Dom_RMinus] *= corr_sig_RPlus_Dom_RMinus
        Barr_sigma[RMinus_Dom_RPlus] *= corr_sig_RMinus_Dom_RPlus
        Barr_sigma[RMinus_Dom_Barr] *= corr_sig_RMinus_Dom_Barr
        Barr_sigma[RPlus_Dom_Barr] *= corr_sig_RPlus_Dom_Barr

        ### calculate the Supercoiling generated in each domaine
        # RNAPs_genSC_all : contains an array of RNAPs_genSC that should be added or substracted from each domaine
        RNAPs_genSC_all = RNAPs_genSC/Dom_size
        # update the value of sigma
        Barr_sigma[Barr_Dom_RPlus] -= RNAPs_genSC_all[Barr_Dom_RPlus]
        Barr_sigma[Barr_Dom_RMinus] += RNAPs_genSC_all[Barr_Dom_RMinus]
        Barr_sigma[RPlus_Dom_RMinus] += 2*RNAPs_genSC_all[RPlus_Dom_RMinus]
        Barr_sigma[RMinus_Dom_RPlus] -= 2*RNAPs_genSC_all[RMinus_Dom_RPlus]
        Barr_sigma[RMinus_Dom_Barr] -= RNAPs_genSC_all[RMinus_Dom_Barr]
        Barr_sigma[RPlus_Dom_Barr] += RNAPs_genSC_all[RPlus_Dom_Barr]
        # We shall consider the case in which we'll have one RNAPol
        # transcribing from right or the left and without any existing barrier on the other side
        # i : you can relace 'if' with None_Dom_RPlus = np.where((len(Barr_type)==1) & (Barr_type_ahead==1))
        if len(Barr_type)==1:
            # ____________O+_____
            None_Dom_RPlus = np.where(Barr_type_ahead==1)
            # ____________O-_____
            None_Dom_RMinus = np.where(Barr_type_ahead==-1)
            # __O+_______________
            RPlus_Dom_None = np.where(Barr_type_ahead==1)
            # __O-_______________
            RMinus_Dom_None = np.where(Barr_type_ahead==-1)
            # update the value of sigma (one Barr case)
            Barr_sigma[None_Dom_RPlus] -= RNAPs_genSC_all[None_Dom_RPlus]
            Barr_sigma[None_Dom_RMinus] += RNAPs_genSC_all[None_Dom_RMinus]
            Barr_sigma[RPlus_Dom_None] += RNAPs_genSC_all[RPlus_Dom_None]
            Barr_sigma[None_Dom_RMinus] -= RNAPs_genSC_all[None_Dom_RMinus]

        # Now calc_sigma
        Barr_sigma = calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO, TOPO_CTE, DELTA_T)

        try:
            mean_sig_wholeGenome = np.sum(Barr_sigma*Dom_size)/genome
        except ValueError:
            # in case we have one RNAPol transcribing (No Barr_fix)
            # ________O______________
            mean_sig_wholeGenome = (Barr_sigma[0]+Barr_sigma[1])/2

        # Update the initiation rate
        init_rate = f_init_rate(tr_rate, sigma_tr_start, sigma_t, epsilon, m)

        if t%OUTPUT_STEP == 0:
            tt=int(t//OUTPUT_STEP)
            # save all informations to npz file
            # RNAPs_info
            save_RNAPs_info[:, 0, tt] = RNAPs_tr
            save_RNAPs_info[:, 1, tt] = RNAPs_pos
            # save the number of hooked RNAPs
            save_nbr_RNAPs_hooked[tt] = np.size(RNAPs_hooked_id)
            # tr_info
            save_tr_info[:, 0, tt] = tr_nbr
            save_tr_info[:, 1, tt] = init_rate

        save_Dom_sigma.append(Barr_sigma)
        save_Dom_size.append(Dom_size)
        save_Barr_pos.append(Barr_pos)
        save_mean_sig_wholeGenome.append(mean_sig_wholeGenome)

    save_Dom_sigma = np.array(save_Dom_sigma)
    save_Dom_size = np.array(save_Dom_size)
    save_Barr_pos = np.array(save_Barr_pos)
    save_mean_sig_wholeGenome = np.array(save_mean_sig_wholeGenome)

    # before saving we check whether we're resuming or not
    # in order to edit the output_path

    # if the output directory isn't specified and we're not resuming the simulation
    if first_output_path==None and resume==False:
        # create another directory inside 'pth' called
        # 'first_output' and save the results there
        output_dir = "first_output"
        # make sure that the output direcory exists, and create one if it doesn't
        os.makedirs("%s/first_output" %pth, exist_ok=True)
        # set the default output path
        output_path=pth+"first_output"

    # otherwise, if we're resuming the simulation
    elif resume_output_path==None and resume==True:
        # create another directory inside 'pth' called
        # 'resume_output' and save the results there
        output_dir = "resume_output"
        # make sure that the output direcory exists, and create one if it doesn't
        os.makedirs("%s/resume_output" %pth, exist_ok=True)
        # set the default output path
        output_path=pth+"resume_output"

    try:
        # Copy the params to the output folder
        copy(INI_file, output_path)
    except Exception as e:
        print("Input file was not copied")
        sys.exit(1)


    save_files(output_path, Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma, tr_nbr, tr_times, save_RNAPs_info, save_tr_info, save_Dom_sigma, save_Barr_pos, save_mean_sig_wholeGenome, save_Dom_size, DELTA_X, RNAPs_genSC, RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id, RNAPs_hooked_id, RNAPs_strand, ts_beg, ts_remain, save_nbr_RNAPs_hooked, init_rate, Kon, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC)

    print("Simulation completed successfully !! \nNumber of transcripts : \n")
    for i, v in enumerate(tr_nbr):
        print("Transcript{} : {}".format(i, v))

    return (GFF_file, TSS_file, TTS_file,
            SIM_TIME, RNAPS_NB,
            tr_nbr, tr_times, init_rate,
            RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id,
            save_RNAPs_info, save_tr_info, save_Dom_sigma, save_Barr_pos,
            cov_bp, tr_end)


if __name__ == '__main__':
    try:
        INI_file=sys.argv[1] #sys.argv[1]           # e.g "../analysis_scripts/example/params.ini"
        # First simulation
        start_transcribing(INI_file)
        # or you can specify the output path
        #start_transcribing(INI_file, output_path)

        # Resuming
        # uncomment this two lines if you want to resume the simulation
        # 'first_output_path' means the path from which the script will get the npz files
        #first_output_path=sys.argv[2]        # e.g "../analysis_scripts/example/first_output"
        #start_transcribing(INI_file, first_output_path, resume=True)
    except configparser.NoSectionError:
        print("Error ! Please check the path to the paraneters files !")
        sys.exit(1)
    except PermissionError:
        print("Permission denied ! Please check the directory to the output files !")
        sys.exit(1)
    except (FileNotFoundError, NameError):
        print("Error ! Please check the directory to the output files !")
        sys.exit(1)
