#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:26:33 2015

@author: Changenet Alexandre
"""

if __name__ == "__main__":
    
    import numpy as np
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description="Simulations for validation")
    parser.add_argument("-r", "--Nrepetitions", type=int, default=100,
                        help="number of repetitions")
    parser.add_argument("-n", "--Nsamples", type=int, default=10000,
                        help="number of samples to take; set to 0 to take\
                        all possible samples; default: 0")
    parser.add_argument("-t", "--theta", type=float, default=0.001,
                        help="theta par base") 
    
    args = parser.parse_args()
    
    Nrepetitions = args.Nrepetitions
    Nsamples = args.Nsamples
    theta_base = args.theta
    
    T = [0.1, 1., 2.]
    alpha = [0.1, 2, 10, 50]
    M = [1, 0.1, 50]
    n = [10, 2, 100]
    long_locus = [100, 200]
    dmin = 10000


    for i in range(0,len(long_locus)):
    	for j in range(0,len(T)):
    		for k in range(0,len(alpha)):
			for l in range(0,100):

    				os.system("mkdir ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/")
   
    				os.system("./ms 2 " + str(Nsamples) + " -t " + str(theta_base * long_locus[i]) + " -eN " + str(0.5 * T[j]) + " " + str(alpha[k]) + " >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+"mssim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "-python.ms")
    				os.system("python ./msfile_histo.py -n " + str(Nsamples) + " -t " + str(0.5 * theta_base * long_locus[i]) + " -o ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+ "mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "-python ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+"mssim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "-python.ms")
    				os.system("sed -e '1,5 d' ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/" + "mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "-python.ndff >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/" + "mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "simul_" + str(l) + ".txt")

    				os.system("./scrm 2 " + str(Nsamples) + " -t " + str(theta_base * long_locus[i]) + " -eN " + str(0.5 * T[j]) + " " + str(alpha[k]) + " >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+"scrmsim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "-python.ms")
    				os.system("python ./msfile_histo.py -n " + str(Nsamples) + " -t " + str(0.5 * theta_base * long_locus[i]) + " -o ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+ "scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "-python ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/"+"scrmsim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "-python.ms")
    				os.system("sed -e '1,5 d' ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/" + "scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "-python.ndff >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_T\=" + str(T[j]) +"_alpha\=" + str(alpha[k]) + "/" + "scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_T\=" + str(T[j]) + "_alpha\=" + str(alpha[k]) + "simul_" + str(l) + ".txt")

    	for k in range(0,len(M)):
    		for j in range(0,len(n)):
			for l in range(0,100):


    				os.system("mkdir ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/")

    
    				os.system("./ms 2 " + str(Nsamples) + " -t " + str(theta_base * long_locus[i]) + " -I " + str(n[j]) + " 2" + " 0" * ((n)[j]-1) + " " + str(M[k]) + " >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "mssim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "-python.ms")
    				os.system("python ./msfile_histo.py -n " + str(Nsamples) + " -t " + str(0.5 * theta_base * long_locus[i]) + " -o ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/"+"mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "-python ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "mssim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "-python.ms")
    				os.system("sed -e '1,5 d' ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "-python.ndff >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "mssim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "simul_" + str(l) + ".txt")


    				os.system("./scrm 2 " + str(Nsamples) + " -t " + str(theta_base * long_locus[i]) + " -I " + str(n[j]) + " 2" + " 0" * ((n)[j]-1) + " " + str(M[k]) + " >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "scrmsim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "-python.ms")
    				os.system("python ./msfile_histo.py -n " + str(Nsamples) + " -t " + str(0.5 * theta_base * long_locus[i]) + " -o ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/"+"scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "-python ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "scrmsim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "-python.ms")
    				os.system("sed -e '1,5 d' ./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "-python.ndff >./sim_indep_" + str(Nsamples) + "frag_size" + str(long_locus[i]) + "_n\=" + str(n[j]) +"_M\=" + str(M[k]) + "/" + "scrmsim_indep_histo_" + str(Nsamples) + "frag_size" + str(long_locus[i]) +"_n\=" + str(n[j]) + "_M\=" + str(M[k]) + "simul_" + str(l) + ".txt")

