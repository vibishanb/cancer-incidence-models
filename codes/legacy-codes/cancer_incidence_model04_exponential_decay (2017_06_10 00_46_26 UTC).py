""" Import of functions for drawing random samples and plotting """
import numpy.random as np
import matplotlib.pyplot as plt
from math import log
import xlsxwriter

""" Initial conditions and parameters """
N_pop = 500000 #Size of the whole population of individuals
p = 0.25 * 10 ** -8 #Probability of mutation in a given locus per cell per generation
cancer_count = [0] * 81 #Overall number of cases of cancer observed at every age-initial array of 100 zeros
pop_size = range ( 0, N_pop ) #Iteration variable for overall simulation
growth_adv_dist = np.exponential (0.5, N_pop) #Exponentially decreasing distribution of potential mutation advantage in the population


""" Main simulation of cancer in the population """
for i in pop_size:
    age = 20 #Starting age of the individual
    n = 8 * 10 ** 9 #Number of stem cells in a saturated tissue
    K = n #Carrying capacity for stem cell population
    n_mut = 0 #Initial number of mutations in the cells
    mut_pop = 0 #Number of mutant stem cells arising in the population
    g = growth_adv_dist [i] #One value from the distribution for every individual    
    diff = 0 #index variable for cancer_count
    p_mut = 1 - ( (1-p) ** n ) #Initial probability that at least one mutation arises in the population
    
    
    while age <= 100:
        a = np.random_sample ( )
        if p_mut > a: #Simulation of stochastic mutation acquisition
            n_mut += 1
            mut_pop = mut_pop + 1
            if n_mut < 5:
                mut_pop = mut_pop + ((mut_pop*g) * ( 1 - (mut_pop/K) )) #One step of logistic growth for mutants
                p_mut = 1 - ( (1-p) ** mut_pop ) #Probability changes with the number of mutant cells in the population
            elif n_mut == 5: #Threshold number of mutations for cancer to occur
                cancer_count [diff] += 1 #The case count is added at the corresponding age in parallel with the overall run
                break
        else:
            mut_pop = mut_pop + ((mut_pop*g) * ( 1 - (mut_pop/K) )) #One step of logistic growth for mutants
            p_mut = 1 - ( (1-p) ** mut_pop ) #Probability changes with the number of mutant cells in the population
        age += 1
        diff = age - 20


""" Final calculations and plotting """
#Total number and percentage of cancer cases in the population for each run
count = 0 
for a in cancer_count:
    count += a

cancer_percent = ( count / N_pop ) * 100
print ( cancer_percent )


# Normalization for survivorship
num_surv = [0] * 81
cancer_frac = [0] * 81
num_surv [0] = N_pop

for i in range ( 0, 80 ):
    num_surv [i+1] = num_surv [i] - cancer_count [i] #Number of individuals surviving at given age

for i in range ( 0, 81 ):    
    cancer_frac [i] = cancer_count [i] / ( cancer_count [i] + num_surv [i] ) #Fraction of surviving population getting cancer

""" # Plotting
age = range (20, 101)
plt.plot ( age, cancer_frac )
plt.xlabel ( 'Age' )
plt.ylabel ( 'Normalized fraction of cancer cases' )
plt.show ( ) """

# Data export to Excel
workbook = xlsxwriter.Workbook ( 'Model.xlsx' )
worksheet = workbook.add_worksheet ( )
row = 0
col = 0

for a in cancer_frac:
    worksheet.write ( row, col, a )
    row += 1

workbook.close ( )
