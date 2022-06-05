#!/usr/bin/env python
# coding: utf-8

# # Dimerización de $NO_2$

# ```{warning}
# Si está utilizando Google Colab o la ejecución en línea, debe de ejecutar al inicio el siguiente código
# ~~~python
# !pip install pyscf
# !pip install geomeTRIC
# ~~~
# ```

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

# In[1]:


# Importe PySCF y su Optimizador


# In[2]:


# Descomentar estas líneas si está en modo online

#!pip install pyscf
#!pip install geomeTRIC

import pyscf
from pyscf import scf
from pyscf import dft
from pyscf.geomopt.geometric_solver import optimize


# **Pregunta 1.** Calcule la energía de la molécula de $N_2O_4$ con HF y la base aug-cc-pvdz. 

# In[3]:


N2O4 = pyscf.gto.M(atom="""
   N       -4.84638        1.76109        0.00000
   N       -3.46888        1.78415        0.00000
   O       -2.82385        2.93169       -0.00000
   O       -2.85055        0.76276        0.00000
   O       -5.46471        2.78248        0.00000
   O       -5.49141        0.61355        0.00000
""",basis="aug-cc-pVDZ")
N2O4.max_memory = 2000
N2O4.build()

rhf = scf.RHF(N2O4)

N2O4_eq = optimize(rhf)

rhf = scf.RHF(N2O4_eq)
n2o4 = rhf.kernel()


# **Pregunta 2.** Calcule la energía de la molécula de $NO_2$ con HF y la base aug-cc-pvdz. [Complete donde haga falta - Reemplace las X]

# In[7]:


NO2 = pyscf.gto.M(atom="""
   N       -4.39539        1.87380        0.00000
   O       -3.90978        3.09520       -0.00000
   O       -3.65594        0.93810        0.00000
""",spin=1,basis="aug-cc-pVDZ")
NO2.max_memory = 2000
NO2.build()

uhf = scf.UHF(NO2)

NO2_eq = optimize(uhf)

uhf = scf.UHF(NO2_eq)
no2 = uhf.kernel()


# **Pregunta a.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según HF.

# In[8]:


(2*no2-n2o4)*2625.5


# **Pregunta 3.** Calcule la energía de la molécula de $N_2O_4$ con DFT B3LYP y la base aug-cc-pvdz.

# In[9]:


N2O4 = pyscf.gto.M(atom="""
   N       -4.84638        1.76109        0.00000
   N       -3.46888        1.78415        0.00000
   O       -2.82385        2.93169       -0.00000
   O       -2.85055        0.76276        0.00000
   O       -5.46471        2.78248        0.00000
   O       -5.49141        0.61355        0.00000
""",basis="aug-cc-pVDZ")
N2O4.max_memory = 2000
N2O4.build()

rks = dft.RKS(N2O4)
rks.xc = "B3LYP"
N2O4_eq = optimize(rks)

rks = dft.RKS(N2O4_eq)
rks.xc = "B3LYP"
n2o4 = rks.kernel()


# **Pregunta 4.** Calcule la energía de la molécula de $NO_2$ con DFT B3LYP y la base aug-cc-pvdz.

# In[10]:


NO2 = pyscf.gto.M(atom="""
   N       -4.39539        1.87380        0.00000
   O       -3.90978        3.09520       -0.00000
   O       -3.65594        0.93810        0.00000
""",spin=1,basis="aug-cc-pVDZ")
NO2.max_memory = 2000
NO2.build()

uks = dft.UKS(NO2)
uks.xc = "B3LYP"
NO2_eq = optimize(uks)

uks = dft.UKS(NO2_eq)
uks.xc = "B3LYP"
no2 = uks.kernel()


# **Pregunta b.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según DFT B3LYP.

# In[11]:


(2*no2-n2o4)*2625.5


# ## Referencias

# - Parrish, R. M.; Burns, L. A.; Smith, D. G. A.; Simmonett, A. C.; DePrince, A. E.; Hohenstein, E. G.; Bozkaya, U.; Sokolov, A. Y.; Di Remigio, R.; Richard, R. M.; et al. **Psi4 1.1: An Open-Source Electronic Structure Program Emphasizing Automation, Advanced Libraries, and Interoperability.** Journal of Chemical Theory and Computation 2017, 13 (7), 3185–3197.
# 
