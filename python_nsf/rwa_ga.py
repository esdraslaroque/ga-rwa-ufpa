#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
#
# GA: RWA with GOF
# Genetic Algorithm: 
# Routing and Wavelength Assignment with General Objective Function
#
# Copyright 2016 Empresa Brasileira de Pesquisa Agropecuária (Embrapa)
#
# Authors: April 2016
# Cassio Trindade Batista - cassio.batista.13@gmail.com

# Last revised on April, 2017

# REFERENCES:
# [1] 
# Afonso Jorge F. Cardoso et. al., 2010
# A New Proposal of an Efficient Algorithm for Routing and Wavelength 
# Assignment (RWA) in Optical Networks

import sys
import info
import nsf

#from termcolor import colored
from operator import itemgetter

import math
import random
import numpy as np

# create population of random individuals
def make_chromosome(nsfnet, start_router, end_router, allels):
	count = 0
	# 1. start from source node
	current_router = start_router
	chromosome = [allels.pop(allels.index(current_router))]
	while len(allels):
		# 2. randomly choose, with equal probability, one of the nodes that
		# is SURELY connected to the current node to be the next in path
		next_router = random.choice(allels)

		# SURELY: check whether there is an edge/link/connection or not
		if nsfnet[current_router][next_router] > 0.0: 
			# 3. if the chosen node hasn't been visited before, mark it as the
			# next in the path (gene). Otherwise find another node
			current_router = next_router
			chromosome.append(allels.pop(allels.index(current_router)))

			# 6. do this until the destination node is found
			if current_router == end_router:
				break

			count = 0
		else:
			# max trials to find a valid path: average of 100 chances per gene
			count += 1
			if count > 100:
				chromosome = False
				break

	if chromosome and len(chromosome) > info.NSF_NUM_NODES:
		chromosome = False

	return chromosome

# evaluate via GOF policies
def evaluate(R, N):
	def get_wave_availability(k, n):
		return (int(n) & ( 1 << k )) >> k

	l = len(R)-1
	L = []
	for w in xrange(1, info.NSF_NUM_CHANNELS+1):
		num = 0
		for i in xrange(l):
			rcurr = R[i]
			rnext = R[i+1]
			num += w * get_wave_availability(w-1, N[rcurr][rnext])
		L.append(round(num/float(w*l), 2))

	wl_avail = 0
	for label in L:
		if label == 1.0:
			wl_avail += 1

	first_wl = float('nan')
	if wl_avail:
		first_wl = L.index(1)

	return L, wl_avail, l+1, first_wl # GOF Label, λ available and route's length

# tournament selection
# [[chrom], [L], wl_avail, r_len]
def select_tour(population, Tc, times=3):
	""" Tournament """
	parents = []
	while len(population)-1:
		# choose and pop a random candidate from population
		# it pops because it cannot be selected twice
		candidate = [random.choice(population)]
		if Tc > random.random():
			for tourn_times in xrange(times):
				# choose another candidate and compare two fitnesses
				# the winner becomes the top candidate; loser is eliminated
				candidate.append(random.choice(population))
				if candidate[0][2] >= candidate[1][2]:
					candidate.remove(candidate[1])
				else: # fit[1] > fit[0] .:. candidate[0] eliminated
					candidate.remove(candidate[0])
			parents.append(candidate[0][0])
		population.remove(candidate[0]) # A: no! never the same dad/mum

	return parents

# ranking selection
def select_wheel(population, Tc):
	""" Ranking """
	parents = []
	rank_range = 0
	for ind in population:
		rank_range += ind[5] # cummulative sum of ranks

	max_spins = len(population)
	for spins in xrange(max_spins):
		if Tc > np.random.random():
			cumsum, roulette = 0, np.random.randint(low=0, high=rank_range+1)
			for ind in reversed(population):
				cumsum += ind[5]
				if roulette < cumsum:
					parents.append(ind[0])
					break

	# make it even
	if len(parents) % 2 == 1:
		cumsum, roulette = 0, np.random.randint(low=0, high=rank_range+1)
		for ind in reversed(population):
			cumsum += ind[5]
			if roulette < cumsum:
				parents.append(ind[0])
				break

	return parents

# "one-point" crossover function
def cross(parents):
	""" One Point """
	children = []
	while len(parents)-1 > 0:
		# choose parents and make sure they are differente ones
		dad = random.choice(parents)
		mom = random.choice(parents)
		parents.remove(dad)
		
		# avoid crossing twins: check if parents are the same individual
		for i in xrange(10):
			if dad == mom:
				mom = random.choice(parents)
			else:
				parents.remove(mom)
				break

		# common nodes between father and mother, excluding source and target
		index_routers = []
		for gene in dad[1:len(dad)-1]:
			if gene in mom[1:len(mom)-1]:
				index_routers.append([dad.index(gene), mom.index(gene)])

			if len(index_routers):
				# randomly choose a common node index to be the crossover point
				common_router = random.choice(index_routers)
				children.append(dad[:common_router[0]] + mom[common_router[1]:])
				children.append(mom[:common_router[1]] + dad[common_router[0]:])

		if not len(children):
			children = False

		return children

# "one-point" mutation function
def mutate(nsfnet, normal_chrom):
	# DO NOT perform mutation if:
	# route has only one link which directly connects source to target
	if len(normal_chrom) == 2:
		return normal_chrom

	trans_chrom = list(normal_chrom) # CAN'T NORMALLY COPY THIS

	# choose a random mutation point, excluding the first and the last
	geneid = random.randrange(1, len(normal_chrom)-1)

	# extract or pop() source and target nodes from chromosome
	start_router = trans_chrom.pop(geneid)
	end_router   = trans_chrom.pop()

	# remove all genes after mutation point
	for gene in xrange(geneid, len(trans_chrom)):
		trans_chrom.pop()

	# alphabet: graph vertices that are not in genes before mutation point
	allels = [start_router, end_router]
	allels += [a for a in range(info.NSF_NUM_NODES) if a not in trans_chrom]

	# create a new route R from mutation point to target node
	R = make_chromosome(nsfnet, start_router, end_router, allels)

	# check if new route/path is valid
	if R:
		trans_chrom += R
	else:
		trans_chrom = list(normal_chrom)

	return trans_chrom

# [[chrom], [L], wl_avail, r_len, first_wl, rank]
def cassort(A):
	B = sorted(sorted(sorted(A, 
				key=lambda x:x[3], reverse=False), # route lenght: shortest path first
				key=lambda x:x[4], reverse=False), # least weight: first-fit
				key=lambda x:x[2], reverse=True)   # least congested first 

	# assign rank
	for rank in xrange(len(B)):
		B[rank][5] = len(B)-rank

	return B

# main Genetic Algorithm function
def rwa_ga(N, A, T, holding_time, popsize, Tc, Tm):
	# generates initial population with random but valid chromosomes
	population = [] # [ [[chrom], [L], wl_avail, r_len, first_wl, rank], ... ] 
	trials = 0
	while len(population) < popsize and trials < 300:
		allels = range(info.NSF_NUM_NODES) # router indexes
		chromosome = make_chromosome(A, info.NSF_SOURCE_NODE, info.NSF_DEST_NODE, allels)
		individual = [chromosome, [], 0, 0, float('nan'), float('nan')]
		if chromosome and individual not in population:
			population.append(individual)
			trials = 0
		else:
			trials += 1

	fits = []
	# <GeneticAlgorithm> ------------------------------------------------------
	for generation in range(info.GA_MAX_GEN):
		# perform evaluation (fitness calculation)
		for ind in xrange(len(population)):
			L, wl_avail, r_len, first_wl = evaluate(population[ind][0], N)
			population[ind][1] = L
			population[ind][2] = wl_avail
			population[ind][3] = r_len
			population[ind][4] = first_wl

		# sort population according to wavelength availability
		population = cassort(population)

		# perform selection
		# WARNING: selecting according to ranking through roulette wheel
		#mating_pool = select_tour(list(population), info.GA_MAX_CROSS_RATE)
		mating_pool = select_wheel(list(population), Tc)

		# perform crossover
		offspring = cross(mating_pool)
		if offspring:
			for child in offspring:
				population.pop()
				population.insert(0, [child, [], 0, 0, float('nan'), float('nan')])

		# perform mutation
		for i in xrange(int(math.ceil(Tm*len(population)))):
			normal_ind = random.choice(population)
			trans_ind = mutate(N, normal_ind[0]) # X MEN
			if trans_ind != normal_ind:
				population.remove(normal_ind)
				population.insert(0, [trans_ind, [], 0, 0, float('nan'), float('nan')])

		fit = 0
		for ind in population:
			if ind[2]:
				fit += 1
		fits.append(fit)

	# </GeneticAlgorithm> -----------------------------------------------------

	# perform evaluation (fitness calculation)
	for ind in xrange(len(population)):
		L, wl_avail, r_len, first_wl = evaluate(population[ind][0], N)
		population[ind][1] = L
		population[ind][2] = wl_avail
		population[ind][3] = r_len
		population[ind][4] = first_wl

	# sort population according to wavelength availability
	population = cassort(population)
	
	fit = 0
	for ind in population:
		if ind[2]:
			fit += 1
	fits.append(fit)

	# update NSF graph
	best_route = population[0][0]
	len_route  = population[0][3]
	
	if population[0][2] > 0:
		color = population[0][1].index(1)
		for i in xrange(len_route-1):
			rcurr = best_route[i]
			rnext = best_route[i+1]

			N[rcurr][rnext] -= 2**color
			N[rnext][rcurr] = N[rcurr][rnext] # make it symmetric

			T[rcurr][rnext][color] = holding_time
			T[rnext][rcurr][color] = T[rcurr][rnext][color]

		return 0 # allocated
	else:
		return 1 # blocked

### EOF ###
