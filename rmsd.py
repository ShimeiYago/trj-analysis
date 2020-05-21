#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from utils import atomicWeightsDecimal as wdict


def main():
    parser = argparse.ArgumentParser(description='This script calculates weighted-RMSD of trr-file inputted.')
    parser.add_argument('-t', '--trr', default='input/fitted.trr', help='trr file')
    parser.add_argument('-s', '--top', default='input/topo.gro', help='topology file')
    parser.add_argument('-o', '--out', default='./rmsd.png', help='output path')
    args = parser.parse_args()

    ### load ###
    trj_mdtraj = md.load(args.trr, top=args.top)
    trj = trj_mdtraj.xyz
    topo = trj_mdtraj.topology

    weightlist = make_weightlist(topo)

    ### calucrate weighted-RMSD
    rmsds = [calu_weighted_rmsd(trj[i], trj[0], weightlist) for i in range(trj.shape[0])]

    ### plot ###
    x = range(len(rmsds))

    fig = plt.figure()
    plt.plot(x, rmsds)
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