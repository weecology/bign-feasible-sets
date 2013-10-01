from __future__ import division
from macroeco_distributions import *
import os
import csv
import numpy as np
from mete import get_mete_rad
import multiprocessing

def get_immediate_subdirectories(in_dir):
    """From http://stackoverflow.com/questions/800197/get-all-of-the-immediate-subdirectories-in-python"""
    return [name for name in os.listdir(in_dir)
            if os.path.isdir(os.path.join(in_dir, name))]

def get_SADs(dataset):
    data_dir = './data/' + dataset + '/' + dataset + '-data-cleaned.txt'
    data = np.genfromtxt(data_dir, dtype = 'S25, f8', delimiter = '\t', 
                         names = ['site', 'obs'])
    return data

def get_pred_geom_logser(dataset):
    out_write_geom = open('./data/' + dataset + '/' + dataset + '-obs-pred-geom.txt', 'w')
    out_write_logser = open('./data/' + dataset + '/' + dataset + '-obs-pred-logser.txt', 'w')
    out_geom = csv.writer(out_write_geom, delimiter = '\t')
    out_logser = csv.writer(out_write_logser, delimiter = '\t')
    
    data = get_SADs(dataset)
    data = data[data['obs'] != 0] # Remove rows with zeros
    site_list = np.sort(list(set(data['site'])))
    for site in site_list:
        data_site = data[data['site'] == site]
        S = len(data_site)
        N = sum(data_site['obs'])
        if S > 4 and round(N) > S: 
            cdf = [(S - i + 0.5) / S for i in range(1, S + 1)]
            pred_geom = trunc_geom.ppf(np.array(cdf), S / N, N)
            pred_logser = get_mete_rad(int(S), int(round(N)))[0]
            results_geom = np.zeros((len(data_site), ), dtype = [('f0', 'S25'), ('f1', float), ('f2', int)])
            results_geom['f0'] = np.array([site] * len(data_site))
            results_geom['f1'] = np.array(sorted(data_site['obs'], reverse = True))
            results_geom['f2'] = np.array(pred_geom)
            out_geom.writerows(results_geom)
            results_logser = np.zeros((len(data_site), ), dtype = [('f0', 'S25'), ('f1', float), ('f2', int)])
            results_logser['f0'] = np.array([site] * len(data_site))
            results_logser['f1'] = np.array(sorted(data_site['obs'], reverse = True))
            results_logser['f2'] = np.array(pred_logser)
            out_logser.writerows(results_logser)
    
    out_write_geom.close()
    out_write_logser.close()
    
if __name__ == '__main__':
    pool = multiprocessing.Pool(4)
    list_of_folders = get_immediate_subdirectories('./data')
    pool.map(get_pred_geom_logser, list_of_folders)