#!/usr/bin/env python
# coding: utf-8

# # Partícula en el anillo

# Es el sistema de una partícula moviéndose en una trayectoria de radio constante tal que $x^2 + y^2 = r^2$.
# 
# 
# ```{admonition} Inserto matemático: Hamiltoniano del sistema
# :class: dropdown
# El Hamiltoniano para una partícula en dos dimensiones, tanto en coordenadas cartesianas como en coordenadas polares es
# 
# $$
# H = -\frac{\hbar^2}{2m} \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right) = -\frac{\hbar^2}{2mr^2} \left( \frac{\partial^2}{\partial r^2} + \frac{1}{r} \frac{\partial}{\partial r} + \frac{1}{r^2} \frac{\partial^2}{\partial \phi^2} \right)
# $$
# 
# Como r es constante en el anillo, entonces se eliminan las derivadas respecto a r, y el Hamiltoniano se vuelve más simple
# 
# $$
# H = -\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2}
# $$
# 
# Al sustituir el Hamiltoniano en la ecuación de Schrödinger, se obtiene
# 
# $$
# -\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2} \psi = E \psi
# $$
# ```
# 
# La solución a la ecuación diferencial tiene la forma
# 
# $$
# \psi = Ae^{im_l\phi} + Be^{-im_l\phi}
# $$
# 
# con $m_l = ( 2mr^2 E/\hbar^2 )^{1/2}$. Despejando se obtiene que 
# 
# $$
# E = \frac{\hbar^2 m_l^2}{2mr^2}
# $$
# 
# ```{admonition} Inserto matemático: Condiciones a la frontera
# :class: dropdown
# 
# Debido a que la función de onda debe ser contínua, y a que la partícula se mueve en un anillo, debe de cumplirse la condición cíclica $\psi(\phi) = \psi(\phi+2\pi)$, es decir que al dar una vuelta, la función de onda debe terminar en el mismo punto donde comenzó. Sustituyendo la función de onda en la condición cíclica se obtiene
# 
# $$
# Ae^{im_l\phi} + Be^{-im_l\phi} = Ae^{im_l\phi}e^{im_l2\pi} + Be^{-im_l\phi}e^{-im_l2\pi}
# $$
# 
# Para que la igualdad anterior pueda cumplirse, $m_l$ debe ser un número entero, tal que se cumpla $e^{im_l2\pi}=e^{-im_l2\pi}=1$. Esto origina la cuantización $m_l = {0, \pm 1, \pm 2, \cdots}$.
# ```
# 
# Tras aplicar las condiciones a la frontera, normalizar y con B=0,
# 
# $$
# \psi = \left( \frac{1}{2\pi} \right)^{1/2} e^{i m_l \phi} = \left( \frac{1}{2\pi} \right)^{1/2} cos(m_l \phi) + i \left( \frac{1}{2\pi} \right)^{1/2} sin(m_l \phi)
# $$

# **Grafique la función de onda y su cuadrado para $m_l=1$**

# In[1]:


# Gráfica


# In[2]:


# Grafica

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

#Numero cuantico
ml=1

#Coordenadas polares
phi = np.linspace(0, 2 * np.pi, 100)
r=1.0
psi_r = np.sqrt(1/(2*np.pi))*np.cos(ml*phi)
psi_i = np.sqrt(1/(2*np.pi))*np.sin(ml*phi)

#Grafica de la funcion de onda
fig = plt.figure()
ax = fig.gca(projection='3d')
x = r * np.sin(phi)
y = r * np.cos(phi)
ax.plot(x, y, 0, color='k') #Eje de la grafica
ax.plot(x, y, psi_r, label='psi-R')
ax.plot(x, y, psi_i, label='psi-I')
ax.legend()
plt.show()

#Grafica del cuadrado de la funcion de onda
fig = plt.figure()
ax = fig.gca(projection='3d')
x = r * np.sin(phi)
y = r * np.cos(phi)
ax.plot(x, y, 0, color='k') #Eje de la grafica
ax.plot(x, y, ((psi_r+1J*psi_i)*(psi_r-1J*psi_i)).real, label='$psi^2$')
ax.legend()
plt.show()


# **Referencias**

# - Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005.
# - Zettili, N. Quantum Mechanics: Concepts and Applications, 2nd ed.; Wiley: Chichester, U.K, 2009.
# - Levine, I. N. Quantum Chemistry, 5th ed.; Prentice Hall: Upper Saddle River, N.J, 2000.
# - McQuarrie, D. A.; Simon, J. D. Physical Chemistry: A Molecular Approach; University Science Books: Sausalito, Calif, 1997.
