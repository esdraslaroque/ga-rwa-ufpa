#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
#
# GA: RWA with GOF
# Genetic Algorithm: 
# Routing and Wavelength Assignment with General Objective Function
#
# Copyright 2017 Universidade Federal do Pará (PPGCC UFPA)
#
# Authors: April 2017
# Cassio Trindade Batista - cassio.batista.13@gmail.com

# Last revised on Apr 2017

# REFERENCES:
# [1] 
# Afonso Jorge F. Cardoso et. al., 2010
# A New Proposal of an Efficient Algorithm for Routing and Wavelength 
# Assignment (RWA) in Optical Networks

import sys

import info
import nsf

import numpy as np

from rwa_ga import rwa_ga

def get_wave_availability(k, n):
	return (int(n) & ( 1 << k )) >> k

if __name__ == '__main__':

	for popsize in info.GA_FIX_SIZE_POP:
		for Tc in info.GA_FIX_CROSS_RATE:
			for Tm in info.GA_FIX_MUT_RATE:
				nsf_wave, nsf_adj, nsf_time, nsf_links = nsf.generate()

				blocked_ga  = []
	
				for load in xrange(info.SIM_MIN_LOAD, info.SIM_MAX_LOAD):
					N_ga = nsf_wave.copy()
					T_ga = nsf_time.copy() # holding time
	
					count_block_ga  = 0
	
					for gen in xrange(info.SIM_NUM_GEN):
						until_next   = -np.log(1-np.random.rand())/load
						holding_time = -np.log(1-np.random.rand())
				
						count_block_ga += rwa_ga(N_ga, nsf_adj, T_ga, 
												holding_time, popsize, Tc, Tm)
	
						if info.DEBUG:
							sys.stdout.write('\rLoad: %02d/%02d ' % (load, info.SIM_MAX_LOAD-1))
							sys.stdout.write('Simul: %04d/%04d\t' % (gen+1, info.SIM_NUM_GEN))
							sys.stdout.write('GA:  %04d, ' % count_block_ga)
							sys.stdout.flush()
				
						# Atualiza os todos os canais que ainda estao sendo usados 
						for link in nsf_links:
							i, j = link
							for w in xrange(info.NSF_NUM_CHANNELS):
	
								# GA + GOF
								if  T_ga[i][j][w] > until_next:
									T_ga[i][j][w] -= until_next
									T_ga[j][i][w]  = T_ga[i][j][w]
								else:
									T_ga[i][j][w] = 0
									T_ga[j][i][w] = 0
									if not get_wave_availability(w, N_ga[i][j]):
										N_ga[i][j] += 2**w # free channel
										N_ga[j][i] = N_ga[i][j] 
	
					blocked_ga.append(100.0*count_block_ga/info.SIM_NUM_GEN)
					print 'Done'
	
				if info.DEBUG:
					print '\tGA:  ', blocked_ga
				
				with open('GA_ps%d_tc%.2f_tm%.2f_ch%d.m' % (popsize, Tc, Tm, info.NSF_NUM_CHANNELS), 'a') as f:
					for bp in blocked_ga:
						f.write('%2.2f ' % bp)
					f.write('\n')
			
### EOF ###
