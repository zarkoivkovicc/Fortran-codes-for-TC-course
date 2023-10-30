import matplotlib.pyplot as plt
import numpy as np
from energydiagram import ED
for k in ['benzene','cyclopentyl', 'butadiene', 'naphthalene', 'pyrene'] :
    a = ED()
    a.top_text_fontsize='xx-small'
    energies = np.loadtxt(f"{k}/{k}.eval",skiprows=1)
    for i,j in enumerate(energies) :
        if (i%2 == 0) :
            a.add_level(j,top_text='',right_text='$\pi$' + str(i+1),position='last')
    for i,j in enumerate(energies) :
        if (i%2 != 0 and (i==1)) :
            a.add_level(j,top_text='',right_text='$\pi$' + str(i+1))
        elif (i%2 != 0 and (i!=1)) :
            a.add_level(j,top_text='',right_text='$\pi$' + str(i+1),position='last')
    a.dimension = 0.6
    a.space = 0.5
    a.offset_ratio = 5
    a.plot(ylabel="Energy / eV")
    plt.savefig(f"{k}/energy_{k}.png")
