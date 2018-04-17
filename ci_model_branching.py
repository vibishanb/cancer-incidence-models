import numpy
import numpy.random as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
from numba import jit

"""V1 [started on 5 March]: Starting the same set of simulations as earlier, with an extended focus of producing comparative predictions for 4 different pathways of evolution-linear, branching, neutral and punctuated. V1 will try to simulate a branching mode of evolution within stem cell populations.

Method used-logistic growth in each sub-clonal line, with the intrinsic growth rate given by the distribution of 'g'. Carrying capacity for each lineage will be n - sum(size of other lineages), which incorporates the effect of competition. Transition rates within lineages are random samples from the same overall distribution. """

"""@jit
def p_calculate(population, p):
	Calculates the probability of at least one mutation for each population of cells.
	Parameters:
	-----------
	population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for extant populations]
	p: intrinsic mutation rate
	
	Returns:
	--------
	parr: 1D-array with probability given by 1-(1-p)**population[:,0]

	return 1-(1-p)**population[:,0]"""
	
"""@jit	
def mutation_simulate(parr):
	Based on parr, simulates a stochastic mutation event in each population of cells. 
	Parameters:
	-----------
	parr: 1D-array containing probability values given by p_calculate()
	
	Returns:
	--------
	mut_index: 1D boolean array, given by parr > rand,  where rand = np.random_sample(len(parr))
	
	rand = np.random_sample(len(parr))
	return parr > rand"""
	
@jit
def generate_mutant_pop(population, index, g_val):
	"""Generates the 2D-array for new mutant populations, to be appended to the existing populations
	
	Parameters:
	-----------
	population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for extant populations]
	index: Boolean array pointing to instances of new mutations, given by mutation_simulate()
	g_val: mut_num number of samples from growth rate distribution, gdist, where mut_num is the number of mutations given by mut_index.sum()
	
	Returns:
	--------
	next_arr: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for nmut populations], where,
	number of new mutants, l = index.sum(),
	growth_rate = random sample of size l from gdist,
	size = 1, if growth_rate > 0,
		0, if growth_rate < 0,
	carrying_capacity = same as corresponding population array, and
	stage = population[index,3] + 1
	"""
	positive_index = numpy.where((g_val>0)*index, True, False)
	
	if positive_index.any():
		l = positive_index.sum()
		next_size = numpy.ones(l)
		k = population[:,2].copy()
		next_stage = population[positive_index,3].copy() + 1
		next_arr = numpy.array([[0,0,0,0]]*l, dtype=numpy.float64)

		for a in range(l):
			next_arr[a] = numpy.array([next_size[a], g_val[positive_index][a], k[positive_index][a], next_stage[a]])
			
		for a in range(l):
			population = numpy.append(population, [next_arr[a]], axis=0)
	
	return population.copy()
	
@jit
def k_calculate(population):
	"""Calculates the carrying capacity for each population of cells, as the current k - sum(sizes of all other populations)
	
	Parameters:
	-----------
	population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for extant populations]
	
	Returns:
	--------
	k_new: 1D-array of length len(population), containing carrying capacity values for corresponding populations
	"""
	population[0,2] = population[0,2].copy() - numpy.sum(population[1:,0].copy())
	other_ind = numpy.arange(1, len(population))
	population[other_ind,2] =  population[other_ind,2].copy() + population[other_ind,0].copy() - numpy.sum(population[~other_ind,0].copy())
	return population

@jit
def grow_logistically(population):
	"""Calculates one step of logistic growth for each population of cells
	
	Parameters:
	-----------
	population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for extant populations]
	
	Returns:
	--------
	m_next: 1D-array with population sizes updated with logistic growth based on parameters corresponding to each population"""
	
	m = population[:,0].copy()
	r = population[:,1].copy()
	k = population[:,2].copy()
	m_next = m + (m*r)*(1-(m/k))
	return m_next


p = 10**-7 #Mutation probability
n = 2*10**5 #Stem cell number
Npop = 100000 #Population size
age = 90 #Maximum lifetime per individual
ndiv = 10 #Number of cell divisions per year of lifespan
time = age*ndiv #Total simulation duration
#wts = [0.013818048555027355, 0.0553159434123515, 0.07253241028642805, 0.07303103455912367, 0.07216711636515384, 0.06647847243710951, 0.06452984736662379, 0.07104508339877749, 0.08076197744193335, 0.08185169462960405, 0.07211714069611326, 0.06271758577923968, 0.048454493422914295, 0.038794489715138394, 0.034263609991378986, 0.03177168658747205, 0.026997999757072733, 0.017842803104216928, 0.01550856249432107] #Weights from the US 2000 standard population

cancer_count = numpy.array([0]*age, dtype=float)
gdist = np.normal(-0.5,0.5,time*Npop) #Distribution of growth_rate, g
flag = 0

for i in range(Npop):
	cell_pop = numpy.array([[n,0,n,0]])
	mut_index = numpy.array([[]], dtype=bool)
	for j in range(time):
		cell_pop[:,0] = grow_logistically(cell_pop)
		parr = 1-numpy.float_power((1-p), cell_pop[:,0]) #Calculate mutation occurence probability for all populations
		mut_index = parr>np.random_sample(len(parr)) #Index of all populations in which mutation has occurred
		if mut_index.any():
			flag += 1
			cell_pop[mut_index,0] -= 1 #Sizes of populations in which mutation has occurred reduces by 1
			g_val = np.choice(gdist, len(mut_index), replace=True)
			cell_pop = generate_mutant_pop(cell_pop, mut_index, g_val) #Newly generated mutant array added to the existing pool of populations
		cell_pop = k_calculate(cell_pop) #Carrying capacity recalculated for all populations, new and old
		
		if (cell_pop[:,3]==5).any():
			cancer_count[int(j/ndiv)] += 1
			break	
		
