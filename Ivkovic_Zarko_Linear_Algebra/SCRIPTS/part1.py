import numpy as np
from matplotlib import pyplot as plt
if __name__ == '__main__':
    r2 = np.loadtxt('r_2')
    deg = { i:3 for i in [1,2,3,4]}
    deg[5] = 2
    for i in [1,2,3,4,5]:
        points = np.genfromtxt(f"{i}/points_{i}.data") 
        curve = np.genfromtxt(str(i)+"/fit.data")
        fit_np = np.polynomial.polynomial.Polynomial.fit(points[:,0],points[:,1],deg=deg[i])
        xx, yy = fit_np.linspace(n=200,domain=[min(points[:,0])-0.9,max(points[:,0])+0.9])
        fig, ax = plt.subplots()
        ax.annotate('R² = ' + str("{:.4f}".format(r2[i-1][1])),xy=(points[-1,0]-2.2,points[-1,1]+0.9))
        ax.scatter(points[:,0],points[:,1],label='Points',color='blue')
        ax.plot(curve[:,0],curve[:,1],label='Fit',color='red')
        ax.plot(xx,yy,color='green',label='Fit Numpy',linestyle='dotted')
        ax.legend()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.plot()
        plt.savefig(str(i)+f"/plot_{i}.png")
