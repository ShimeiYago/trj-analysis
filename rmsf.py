#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='This script calculates RMSF of trr-file inputted.')
    parser.add_argument('-t', '--trr', default='input/fitted.trr', help='trr file')
    parser.add_argument('-s', '--top', default='input/topo.gro', help='topology file')
    parser.add_argument('-o', '--out', default='./rmsf.png', help='output path')
    args = parser.parse_args()

    ### load ###
    trj_mdtraj = md.load(args.trr, top=args.top)
    trj = trj_mdtraj.xyz
    topo = trj_mdtraj.topology

    CAindexlist = [atom.index for atom in topo.atoms if atom.name == 'CA']


    ### calucrate RMSF of Ca
    rmsfs = calu_rmsfs(trj[:, CAindexlist, :])


    ### print most largest Ca indexes about RMSF ###
    print('most largest Ca indexes about RMSF')
    print(np.argsort(rmsfs)[::-1][:10])


    ### plot ###
    x = range(1, len(CAindexlist)+1)

    fig = plt.figure()
    plt.plot(x, rmsfs)
    fig.savefig(args.out)


def calu_rmsfs(trj):
    mean_structure = trj.mean(axis=0)
    
    return np.sqrt(np.mean(np.sum(np.square(trj - mean_structure), axis=2), axis=0))


if __name__=='__main__':
    main()