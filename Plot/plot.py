import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#plt.rcParams["font.family"] = "Times New Roman"

mpl.rc('xtick', labelsize=22) 
mpl.rc('ytick', labelsize=22)

data = np.genfromtxt("state_space.txt", skip_header = 9,unpack = True)

fig = plt.figure()

plt.plot(data[0],data[1],'bo', ms=2)


fig.savefig("state_space.jpg")
