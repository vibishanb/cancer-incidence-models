""" Import of modules """
import numpy.random as np
import xlsxwriter
import matplotlib.pyplot as plt

""" Initialization """
N_pop = 100000 #Population size
prob = [0] * 2 #Range of mutational probabilities
a = len ( prob )
for i in range ( 0, a ):
    prob [i] = 10 ** -(i+7)

time = 100 #Duration of each simulation
growth_rate = np.normal ( -1, 0.5, N_pop )
g_track = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ]

mut_pop = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Transition values of mutant population sizes
product = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Product of p and mut_pop will be stored here
mut_prev = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Threshold values of mutant population sizes
cancer_count = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Age-wise incidence of cancer

n = 50 * 10 ** 7 #Stem cell carrying capacity for the tissue

""" Main simulation """

for i in range ( 0, a ):
    p = prob [i]
    for j in range ( 0, N_pop ):
        g = growth_rate [j]
        t = 0 #Index to track time
        n_mut = [0] * time
        m = 0 #Initial mutant population
        g_total = 0 #Total number of mutant cells with a particular set of mutations
        m_prev = 0 #Number of cells in the previous step
        m_hold = 0 #Holding variable for previous population size
        trans = 0 #Number of transitions in the population
        p_initial  = 1 - ( (1-p) ** n ) #Initial probabiltiy of first mutation arising in the population
        g_track [i][0] = g

        if p_initial > np.random_sample ( ): #At t = 0
            n_mut [0] += 1
            m = 1
            p_mut = 1 - ( (1-p) ** m ) 
        else:
            m = 0
            p_mut  = p_initial

        for t in range ( 1, time ): #From t = 1 to end of time
            n_mut [t] = n_mut [t-1]
            m_hold = m
            m += ( ( m*g ) * ( 1 -  ( m/n ) ) )
            g_track [i][t] += g / N_pop
            
            if n_mut [t] > 0:
                p_mut = 1 - ( (1-p) ** m )
            else:
                p_mut = p_initial

            if p_mut > np.random_sample ( ): #g value is the same for all mutations within an individual => first positive mutation has an all-time advantage
                n_mut [t] += 1
                m_prev += m_hold
                g_total += m
                m = 1
                trans += 1
                
                if n_mut [t] == 5: #Threshold number of mutations for cancer is 5
                    cancer_count [i][t] += 1
                    mut_pop [i][t] += ( g_total / trans ) / N_pop
                    mut_prev [i][t] += ( m_prev / trans ) / N_pop
                    break

""" Calculations """
total_num = [0] * a #Total number of cancer cases
ci_rate = [0] * a #Cancer incidence rate
g1 = [0] * time
g2 = [0] * time
count1 = [0] * time
count2 = [0] * time

for j in range ( 0, time ):
    g1 [j] = g_track [0][j]
    count1 [j] = cancer_count [0][j]
    g2 [j] = g_track [1][j]
    count2 [j] = cancer_count [1][j]

for i in range ( 0, a ):
    for j in range ( 0, time ):
        total_num [i] += cancer_count [i][j]
        ci_rate [i] = ( total_num [i] / N_pop ) * 100000 #CI rate is calculated as fraction of cases per 100000 individuals

print ( ci_rate )

""" Export to Excel """
workbook = xlsxwriter.Workbook ( 'model.xlsx' )
worksheet = workbook.add_worksheet ( )

row = 0
col = 0
for i in range ( 0, a ):
    worksheet.write ( row, col, ci_rate [i] )
    col += 1

row = 0
col = 0
for j in range ( 0, time ):
    worksheet.write ( row+2, col, count1 [j] )
    worksheet.write ( row+2, col+1, count2 [j] )
    row += 1

workbook.close ( )
        

        
        
            
        
