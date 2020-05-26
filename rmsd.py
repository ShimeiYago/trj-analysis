#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from utils import atomicWeightsDecimal as wdict


def main():
    parser = argparse.ArgumentParser(description='This script calculates weighted-RMSD of trr-file inputted.')
    parser.add_argument('-t', '--trr', nargs='+', required=True, help='trr files')
    parser.add_argument('-s', '--top', required=True, help='topology file')
    parser.add_argument('-o', '--out', default='./rmsd.png', help='output path')
    parser.add_argument('-u', '--upper_time', type=int, help='Upper time of trajectory to draw')
    args = parser.parse_args()

    ### load ###
    trj_mdtraj_list = [md.load(trrpath, top=args.top) for trrpath in args.trr]
    trj_list = [trj_mdtraj.xyz[:args.upper_time] for trj_mdtraj in trj_mdtraj_list]
    topo = trj_mdtraj_list[0].topology

    weightlist = make_weightlist(topo)

    ### calucrate weighted-RMSD
    rmsds_list = []
    for trj in trj_list:
        rmsds = [calu_weighted_rmsd(trj[i], trj[0], weightlist) for i in range(trj.shape[0])]
        rmsds_list.append(rmsds)

    ### plot ###
    fig = plt.figure()
    for i, rmsds in enumerate(rmsds_list):
        x = range(len(rmsds))
        plt.plot(x, rmsds, label=os.path.basename(args.trr[i]))
    
    plt.xlabel('time (ps)')
    plt.ylabel('weighted RMSD (nm)')
    plt.legend(loc='upper left')
    fig.savefig(args.out)


def calu_weighted_rmsd(target_struct:np.ndarray, ref_struct:np.ndarray, weightlist:list):
    square_distance_list = np.sum(np.square(target_struct-ref_struct), axis=1)
    return sum([weightlist[i]*square_distance_list[i] for i in range(len(weightlist))]) / sum(weightlist)


def make_weightlist(top, weight_key='standard'):
    atomlist = [atom for atom in top.to_dataframe()[0].name]
    wlist = [float(wdict[atom][weight_key]) for atom in top.to_dataframe()[0].element]
    return wlist


if __name__=='__main__':
    main()