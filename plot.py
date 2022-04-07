import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt('results/p.txt')
rho = np.loadtxt('results/rho.txt')
u = np.loadtxt('results/u.txt')
e = np.loadtxt('results/e.txt')

t = p[:, 0]
p = p[:, 1:]
rho = rho[:, 1:]
u = u[:, 1:]
e = e[:, 1:]


x = np.linspace(0.2, 20 - 0.02, 100)

i = 1500
# print(t[i])
print(t.shape)


fig, ax = plt.subplots(2, 2)

ax[0, 0].plot(x, rho[i, :], linestyle='none', marker='.')
ax[0, 0].grid()
ax[0, 0].set_title('rho')

ax[0, 1].plot(x, u[i, :], linestyle='none', marker='.')
ax[0, 1].grid()
ax[0, 1].set_title('u')

ax[1, 0].plot(x, e[i, :], linestyle='none', marker='.')
ax[1, 0].grid()
ax[1, 0].set_title('e')

ax[1, 1].plot(x, p[i, :], linestyle='none', marker='.')
ax[1, 1].grid()
ax[1, 1].set_title('p')

plt.show()
