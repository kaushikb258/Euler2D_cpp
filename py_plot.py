import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from matplotlib import cm

filename = 'x_output.csv'
X = pd.read_csv(filename,header=None)
X = np.array(X)

filename = 'y_output.csv'
Y = pd.read_csv(filename,header=None)
Y = np.array(Y)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# ii = 1: rho 
# ii = 2: p 
# ii = 3: vel 

ii = raw_input("enter ii (1/2/3): ")

if(ii=='1'):
 filename = 'rho_output.csv'
elif(ii=='2'):
 filename = 'p_output.csv'
elif(ii=='3'):
 filename = 'vel_output.csv'

print("filename=",filename)

Z = pd.read_csv(filename,header=None)
Z = np.array(Z) 

ax.plot_surface(X,Y,Z, rstride=1, cstride=1, linewidth=0, antialiased=False)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('fluid variable')

plt.show()
