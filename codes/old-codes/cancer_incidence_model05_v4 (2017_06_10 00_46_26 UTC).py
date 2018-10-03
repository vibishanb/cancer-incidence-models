""" Import of modules """
import numpy.random as np
import xlsxwriter
import matplotlib.pyplot as plt

""" Initialization """
N_pop = 5000000 #Population size
prob = [0] * 5 #Range of mutational probabilities
a = len ( prob )
for i in range ( 0, a ):
    prob [i] = ( 2/(i+1) ) * 10 ** -7

growth_rate = np.normal ( 0, 0.5, N_pop )

time = 100 #Duration of each simulation
mut_pop = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Transition values of mutant population sizes
product = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Product of p and mut_pop will be stored here
mut_prev = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Threshold values of mutant population sizes
cancer_count = [ [0 for i in range ( 0, time )] for j in range ( 0, a ) ] #Age-wise incidence of cancer

n = 10 ** 7 #Stem cell carrying capacity for the tissue

""" Main simulation """

for i in range ( 0,a ):
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
        p_mut  = 1 - ( (1-p) ** n ) #Initial probabiltiy of first mutation arising in the population

        if p_mut > np.random_sample ( ): #At t = 0
            n_mut [0] += 1
            m = 1
        else:
            m = 0

        for t in range ( 1, time ): #From t = 1 to end of time
            n_mut [t] = n_mut [t-1]
            m_hold = m
            m += ( ( m*g ) * ( 1 -  ( m/n ) ) )
            p_mut = 1 - ( (1-p) ** m )

            if p_mut > np.random_sample ( ):
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

""" Calculations
total_num = [0] * a
for i in range ( 0, a ):
    for j in range ( 0, time ):
        total_num [i] += cancer_count [i][j]

print (total_num ) """

""" Plotting
age = range ( 0, time )
count1 = [0] * time
count2 = [0] * time
count3 = [0] * time
count4 = [0] * time
count5 = [0] * time
for i in range ( 0, time ):
    count1 [i] = cancer_count [0][i]
    count2 [i] = cancer_count [1][i]
    count3 [i] = cancer_count [2][i]
    count4 [i] = cancer_count [3][i]
    count5 [i] = cancer_count [4][i]
    
plt.plot ( age, count1 )
plt.xlabel ( 'Generation' )
plt.ylabel ( 'Number of cancer cases (un-normalized)' )
plt.show ( ) """

""" Export to Excel """
workbook = xlsxwriter.Workbook ( 'Model4.2.2.xlsx' )
worksheet = workbook.add_worksheet ( )
row = 0
col = 0

for i in range ( 0, a ):
    for j in range ( 0, time ):
        worksheet.write ( row, col, mut_pop [i][j] )
        worksheet.write ( row, col + 7, mut_prev [i][j] )
        worksheet.write ( row, col + 14, cancer_count [i][j] )
        row += 1
    col += 1
    row = 0

workbook.close ( )
        

        
        
            
        
