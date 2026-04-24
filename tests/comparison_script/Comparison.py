import numpy as np


if __name__ == '__main__':
    
    uRISTRA = np.fromfile('bar_linear_displacements_SimStudio.bin', dtype='<f8').reshape(-1, 3)
    uFlorence = np.loadtxt('bar_linear_displacements_florence.txt', delimiter=',', skiprows=1)
    

    uDiff = uFlorence - uRISTRA
    print(uDiff)
    print("Max absolute difference in displacements:", np.max(np.abs(uDiff)))
    print("Mean absolute difference in displacements:", np.mean(np.abs(uDiff)))