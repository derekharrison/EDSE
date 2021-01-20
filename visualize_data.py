import numpy as np
import matplotlib.pyplot as plt

num_data_file = 'data.txt'
x_p, vec_ana, eig_vectors, vec_ana_1, eig_vectors_1, vec_ana_2, eig_vectors_2 = np.genfromtxt(num_data_file, unpack=True)


fig = plt.figure()
plt.plot(x_p, vec_ana, 'r--', x_p, eig_vectors, 'bs')
plt.show()
plt.plot(x_p, vec_ana_1, 'r--', x_p, eig_vectors_1, 'bs')
plt.show()
plt.plot(x_p, vec_ana_2, 'r--', x_p, eig_vectors_2, 'bs')
plt.show()
