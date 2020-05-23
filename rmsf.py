#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='This script calculates RMSF of trr-file inputted.')
    parser.add_argument('-t', '--trr', nargs='+', required=True, help='trr files')
    parser.add_argument('-s', '--top', default='input/topo.gro', help='topology file')
    parser.add_argument('-o', '--out', default='./rmsf.png', help='output path')
    args = parser.parse_args()

    ### load ###
    trj_mdtraj_list = [md.load(trrpath, top=args.top) for trrpath in args.trr]
    trj_list = [trj_mdtraj.xyz for trj_mdtraj in trj_mdtraj_list]
    topo = trj_mdtraj_list[0].topology

    CAindexlist = [atom.index for atom in topo.atoms if atom.name == 'CA']


    ### calucrate RMSF of Ca
    rmsfs_list = []
    for trj in trj_list:
        rmsfs = calu_rmsfs(trj[:, CAindexlist, :])
        rmsfs_list.append(rmsfs)


    ### plot ###
    fig = plt.figure()
    for i, rmsfs in enumerate(rmsfs_list):
        x = range(1, len(CAindexlist)+1)
        plt.plot(x, rmsfs, label=os.path.basename(args.trr[i]))

    plt.xlabel('C-alpha Index')
    plt.ylabel('RMSF (nm)')
    plt.legend()
    fig.savefig(args.out)


def calu_rmsfs(trj):
    mean_structure = trj.mean(axis=0)
    
    return np.sqrt(np.mean(np.sum(np.square(trj - mean_structure), axis=2), axis=0))


if __name__=='__main__':
    main()