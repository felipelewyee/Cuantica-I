#!/usr/bin/env python
# coding: utf-8

# # Dimerización de $NO_2$

# El $N_2O_4$ es un compuesto que se encuentra en la atmósfera. Con los cambios de temperatura se descompone en el radical $NO_2$ mediante la reacción:
# 
# $N_2O_4 <=> 2NO_2$
# 
# Las geometrías de $N_2O_4$ y $NO_2$ se dan a continuación:
# 
# Molécula: $N_2O_4$ Carga: 0 Multiplicidad: 1
# 
# |Átomo |x (Å)  |y (Å)  |z (Å) |
# |------|-------|-------|------|
# |N     |0.0000 |0.0000 |0.0000|
# |N     |-1.7820|0.0000 |0.0000|
# |O     |0.4516 |1.1010 |0.0000|
# |O     |0.4516 |-1.1010|0.0000|
# |O     |-2.2336|1.1010 |0.0000|
# |O     |-2.2336|-1.1010|0.0000|
# 
# Molécula: $NO_2$ Carga: 0 Multiplicidad: 2
# 
# |Átomo |x (Å)  |y (Å)  |z (Å) |
# |------|-------|-------|------|
# |N     | 0.0000| 0.0000|0.0000|
# |O     | 0.0000| 1.0989|0.4653|
# |O     | 0.0000|-1.0989|0.4653|
# 

# **Pregunta 1.** Calcule la energía de la molécula de $N_2O_4$ con HF y la base aug-cc-pvdz. 

# In[1]:


import psi4
psi4.set_memory("2 gb")
psi4.geometry("""
0 1
   N       -4.84638        1.76109        0.00000
   N       -3.46888        1.78415        0.00000
   O       -2.82385        2.93169       -0.00000
   O       -2.85055        0.76276        0.00000
   O       -5.46471        2.78248        0.00000
   O       -5.49141        0.61355        0.00000
""")
n2o4,wfn=psi4.opt("HF/aug-cc-pvdz", return_wfn=True)
print(n2o4)


# **Pregunta 2.** Calcule la energía de la molécula de $NO_2$ con HF y la base aug-cc-pvdz. [Complete donde haga falta - Reemplace las X]

# In[2]:


import psi4
NO2 = psi4.geometry("""
0 2
   N       -4.39539        1.87380        0.00000
   O       -3.90978        3.09520       -0.00000
   O       -3.65594        0.93810        0.00000
units angstrom
""")
psi4.set_options({'reference': 'uhf'})
no2=psi4.optimize("HF/aug-cc-pvdz")
print(no2)


# **Pregunta a.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según HF.

# In[3]:


(2*no2-n2o4)*2625.5


# **Pregunta 3.** Calcule la energía de la molécula de $N_2O_4$ con DFT B3LYP y la base aug-cc-pvdz.

# In[4]:


import psi4
psi4.set_memory("2 gb")
psi4.geometry("""
0 1
N     0.0000  0.0000  0.0000
N     -1.7820 0.0000  0.0000
O     0.4516  1.1010  0.0000
O     0.4516  -1.1010 0.0000
O     -2.2336 1.1010  0.0000
O     -2.2336 -1.1010 0.0000
""")
n2o4,wfn=psi4.opt("B3LYP/aug-cc-pvdz", return_wfn=True)
print(n2o4)


# **Pregunta 4.** Calcule la energía de la molécula de $NO_2$ con DFT B3LYP y la base aug-cc-pvdz.

# In[5]:


import psi4
NO2 = psi4.geometry("""
0 2
N      0.0000  0.0000 0.0000 
O      0.0000  1.0989 0.4653 
O      0.0000 -1.0989 0.4653 
units angstrom
""")
psi4.set_options({'reference': 'uhf'})
no2=psi4.optimize("B3LYP/aug-cc-pvdz")
print(no2)


# **Pregunta b.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según DFT B3LYP.

# In[6]:


(2*no2-n2o4)*2625.5


# ## Referencias

# - Parrish, R. M.; Burns, L. A.; Smith, D. G. A.; Simmonett, A. C.; DePrince, A. E.; Hohenstein, E. G.; Bozkaya, U.; Sokolov, A. Y.; Di Remigio, R.; Richard, R. M.; et al. **Psi4 1.1: An Open-Source Electronic Structure Program Emphasizing Automation, Advanced Libraries, and Interoperability.** Journal of Chemical Theory and Computation 2017, 13 (7), 3185–3197.
# 
