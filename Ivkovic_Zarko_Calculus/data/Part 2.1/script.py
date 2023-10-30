from pylab import meshgrid
import numpy as np
import matplotlib.pyplot as plt
x = np.arange(0,5,0.1)
y = np.arange(0,5,0.1)
def z_func(x,y):
 return (np.sin(x+y) + (x-y)**2 -1.5*x + 3.5*y + 3)
X1,X2 = meshgrid(x,y)
xx = np.loadtxt('data',usecols=0)
yy = np.loadtxt('data',usecols=1)
zz = np.loadtxt('data',usecols=2)
Y = z_func(X1,X2)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X1, X2, Y)
ax.scatter(xx,yy,color='red')
ax.plot(xx,yy,color='yellow')
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Gradient descent')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.savefig('Figure_1.png')

