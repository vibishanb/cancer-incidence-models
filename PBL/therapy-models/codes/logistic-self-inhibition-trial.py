import pandas as pd
import seaborn as sns

r = 0.1
d = 0.005
K = zeros(1000)
K[0] = 1000
x1, x2 = zeros(1000), zeros(1000)
x1[0] = 10
x2[0] = 10

for i in range(1, 1000):
    x1[i] = x1[i-1] + r*x1[i-1]*(1-(x1[i-1]/K[0]))
    x2[i] = x2[i-1] + r*x2[i-1]*(1-(1/(K[i-1]-x1[i-1])))

plot(x1, arange(1000), x2, arange(1000))
plot(arange(1000), x1, arange(1000), x2)
x2

a = arange(100)
plot(a, 1-(1/a))
