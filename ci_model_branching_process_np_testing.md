# Branching evolution model of cancer incidence

This code simulates a branching process of somatic evolution, with stochastic mutation accumulation and cell
populations undergo logistic growth. Competition between populations is captured through a shared carrying capacity.

## V1 [started on 5 March]
- Intrinsic growth rate given by random samples from a normally-distributed 'g'. Carrying capacity calculated as sum(all other existing populations)-carrying capacity of the focal population.
- Since growth/transition rates are samples from the same distribution, this is the context-independent version of the model.
- Working draft finalised on 17 April.

## np testing [re-started on 1 May]
- Testing of the whole simulation for 10 values each of p and n.
- Produces absolute counts for each combination of n and p value, and age-adjusted incidence rates standardised to the US 2000 standard population.


### Module import and function definition
```
import numpy
import numpy.random as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from tqdm import tqdm_notebook as tqdm
from numba import jit

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
	next_arr: 2-D array of the form [[size, growth_rate, n, stage]...for nmut populations], where,
	number of new mutants, l = index.sum(),
	growth_rate = random sample of size l from gdist,
	size = 1, if growth_rate > 0,
		0, if growth_rate < 0, and
	stage = population[index,3] + 1
	"""
	positive_index = numpy.where((g_val>0)*index, True, False)
	
	if positive_index.any():
		l = positive_index.sum()
		next_size = numpy.ones(l)
		k = population[positive_index,2].copy()
		next_stage = population[positive_index,3] + 1
		#next_arr = numpy.array([[0,0,0,0]]*l, dtype=numpy.float64)

		for a in range(l):
			population = numpy.append(population, numpy.array([[next_size[a], g_val[positive_index][a], k[a], next_stage[a]]]), axis=0)
	population[:,2] = k_calculate(population, n)
	return population.copy()
	
@jit
def k_calculate(population, n):
	"""Calculates the carrying capacity for each population of cells, as the current k - sum(sizes of all other populations)
	
	Parameters:
	-----------
	population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for extant populations]
	
	Returns:
	--------
	new_k: 1D-array of length len(population), containing carrying capacity values for corresponding populations
	"""
	current_sizes = population[:,0].copy()
	all_sizes = numpy.array([numpy.sum(current_sizes)]*len(current_sizes))
	current_k = population[:,2].copy()

	new_k = numpy.array([n]*len(current_sizes)) - all_sizes + current_sizes

	return new_k.copy()

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
	return m_next.copy()

@jit
def drop_dead(old_population):
	"""Checks for depleted cell populations and removes from the simulation
    Parameters:
    -----------
    old_population: 2-D array of the form [[size, growth_rate, carrying_capacity, stage]...for all extant populations]
    
    Returns:
    --------
    new_population: 2-D array without those populations are zero or less;
                    new_population = old_population[~numpy.less_equal(size,0)]"""
    
	all_sizes = old_population[:,0].copy()
	dead = numpy.less_equal(all_sizes, 0)
	new_population = old_population[~dead]
	return new_population.copy()

@jit
def cancer_statistics(cancer_count, wts, age, Npop):
	"""Calculates age-normalised incidence rates standardised to the given standard population; resulting wtd_rate is split in as many age classes as in wts.
    Parameters:
    -----------
    cancer_count: crude counts of age-wise cancer incidence,
    wts: weights from the standard population,
    age: maximum allowed age in the simulation, and
    Npop: population size
    
    Returns:
    --------
    wtd_rate: age-normalised incidence rates, standardised to the given population, and
    total_inc: normalised total incidence, calculated as sum(wtd_rate)
    """
	num_surv, cumul_count, cancer_fract, crude_rate = numpy.zeros(age), numpy.zeros(age), numpy.zeros(age), numpy.zeros(age)
	age_rate = numpy.zeros_like(wts)
	num_surv[0]=Npop
	for t in range(1, age):
		num_surv[t]=num_surv[t-1]-cancer_count[t-1] #Number of individuals surviving at given age
		cumul_count[t]=cumul_count[t-1] + cancer_count[t]

	t=0
	while (t <= (age-1) and num_surv[t] != 0 ):
		cancer_fract[t]=cancer_count[t] / (cancer_count[t]+num_surv[t]) #Fraction of surviving population getting cancer
		crude_rate[t]=cancer_fract[t]*100000
		t+=1

	age_rate[0]=crude_rate[0]
	age_rate[1]=sum(crude_rate[1:4])
	age_rate[-1]=sum(crude_rate[85:len(crude_rate)])
	for i in range(2,18):
		age_rate[i]=sum(crude_rate[(5*(i-1)):(5*(i-1)+4)])

	wtd_rate=wts*age_rate #Age adjusted rate calculation-weighted sum of age-specific rates
	total_inc = sum(wtd_rate)
	return wtd_rate, total_inc
```
### Common initialisation
```
Npop = 100000 #Population size
age = 90 #Maximum lifetime per individual
ndiv = 10 #Number of cell divisions per year of lifespan
time = age*ndiv #Total simulation duration
wts = [0.013818048555027355, 0.0553159434123515, 0.07253241028642805, 0.07303103455912367, 0.07216711636515384, 0.06647847243710951, 0.06452984736662379, 0.07104508339877749, 0.08076197744193335, 0.08185169462960405, 0.07211714069611326, 0.06271758577923968, 0.048454493422914295, 0.038794489715138394, 0.034263609991378986, 0.03177168658747205, 0.026997999757072733, 0.017842803104216928, 0.01550856249432107] #Weights from the US 2000 standard population
```
### Main simulation
```
p_arr = numpy.linspace(10**-8, 10**-5, 5) #Mutation probability
n_arr = numpy.linspace(10**5, 10**8, 5) #Stem cell number
nval = len(p_arr)

all_counts = numpy.array([numpy.zeros_like(cancer_count)], dtype=float)
gdist = np.normal(-0.05,0.1,time*Npop) #Distribution of growth_rate, g
# flag = 0

for a in tqdm(range(len(p_arr)), desc='Mutation rate', leave=False):
    p = p_arr[a]
    for b in tqdm(range(len(n_arr)), desc='Cell number', leave=False):
        n = n_arr[b]
        cancer_count = numpy.zeros(age)
        for i in tqdm(range(Npop), desc='Npop', leave=False):
	        cell_pop = numpy.array([[n,0,n,0]])
	        mut_index = numpy.array([[]], dtype=bool)
	        parr = 1-numpy.float_power((1-p), cell_pop[:,0])
	        mut_index = numpy.greater(parr, np.random_sample(len(parr)))
        
	        for j in range(time):        
		        if mut_index.any():
			        flag += 1
			        cell_pop[mut_index,0] -= 1 #Sizes of populations in which mutation has occurred reduces by 1
			        g_val = np.choice(gdist, len(mut_index), replace=True)
			        cell_pop = generate_mutant_pop(cell_pop, mut_index, g_val) #Newly generated mutant array added to the existing pool of populations

		        cell_pop = drop_dead(cell_pop)
		        cell_pop[:,2] = k_calculate(cell_pop, n) #Carrying capacity calculated for all populations, new and old
		        cell_pop[:,0] = grow_logistically(cell_pop) #One step of logistic growth
		        parr = 1-numpy.float_power((1-p), cell_pop[:,0]) #Calculate mutation occurence probability for all populations
		        mut_index = numpy.greater(parr, np.random_sample(len(parr))) #Index of all populations in which mutation has occurred
		        if (cell_pop[:,3]==5).any():
			        cancer_count[int(j/ndiv)] += 1
			        break
        all_counts = numpy.append(all_counts, [cancer_count], axis=0)
all_counts = all_counts[1:]
```
### Statistical calculations
```
std_rate = numpy.zeros_like([wts]*len(all_counts))
tot_incidence = numpy.zeros(len(all_counts))
for i in range(len(all_counts)):
    std_rate[i], tot_incidence[i] = cancer_statistics(all_counts[i], wts, age, Npop)    
```
### Data export
```
df = pd.DataFrame(std_rate, index = [numpy.ravel([[x,x,x,x,x] for x in p_arr]), numpy.array([n_arr]*5).ravel()])
df.index.names = ['Mutation rate', 'Cell number']
df.columns.names = ['Age classes']

a=[[(5*(i-1)),(5*(i-1)+4)] for i in range(2,18)]
b = ['[0]', '[1,4]']
for i in a:
    b.append(str(i))
b.append(str([85, 89]))
df.columns = b

df.to_excel('branching_process_np_testing_std_rate.xlsx')
```
