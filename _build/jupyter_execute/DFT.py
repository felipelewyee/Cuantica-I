#!/usr/bin/env python
# coding: utf-8

# # Teoría de funcionales de la densidad (DFT)

# La teoría de los funcionales de la densidad se basa en que la densidad electrónica contiene toda la información del sistema, por lo que se puede extraer de ella la energía mediante un funcional
# 
# $$
# E[\rho] = T_s[\rho] + U[\rho] + V_{nuc}[\rho] + E_{xc}[\rho]
# $$
# 
# Donde
# - $E[\rho]$ es la energía del sistema 
# - $T_s[\rho]$ es la energía cinética de Kohn-Sham
# - $U[\rho]$ es la energía de Hartree
# - $V[\rho]$ es la interacción núcleo electrón
# - $E_{xc} [\rho]$ es la energía de intercambio correlación
# 
# En la formulación de Kohn y Sham, la densidad electrónica ($\rho(r)$) se calcula a partir de los orbitales de Kohn-Sham ($\psi^{KS}_i(r)$)
# 
# $$
# \rho(r) = \sum_i^{N} |\psi^{KS}_i(r)|^2
# $$
# 
# 
# ```{warning}
# Los teoremas de Hohenberg y Kohn, y Kohn y Sham prueban que existe un funcional universal que conecta a la densidad electrónica con la energía del sistema, pero este funcional no se conoce, particularmente por la parte de $E_{xc}[\rho]$. Se han realizado diversas propuestas de como construir este funcional, por lo que en la práctica existen cientos de funcionales de DFT.
# ```
# 
# ```{note}
# No existe el mejor funcional, sino que depende del sistema químico que se esté estudiando.
# ```
# 
# En general los funcionales aproximan de diferente manera el intercambio y la correlación, dependiendo de como se aproxime el funcional, estos se pueden clasificar en diferentes categorías. La clasificación fue propuesta por Perdew y se conoce como la `escalera de Jacob`.
# 
# ```{warning}
# Evaluar la contribución del funcional de intercambio-correlación requiere de realizar integración numérica. Esto se hace con un mallado entorno a la molécula. Existen diversos esquemas para colocar los puntos en el mallado, y para seleccionar cuantos puntos poner, estos puede cambiar de un software a otro, e incluso entre diferentes versiones de un mismo software.
# ```

# **Importe PySCF**

# ```{warning}
# Si está utilizando Google Colab o la ejecución en línea, debe de ejecutar al inicio el siguiente código
# ~~~
# !pip install pyscf
# ~~~
# ```

# In[1]:


#Importe PySCF


# In[2]:


# Descomentar estas líneas si está en modo online

#!pip install pyscf

import pyscf
from pyscf import dft


# Para ejemplificar el uso de estos funcionales, declare la molécula de agua.
# 
# ```
# h2o = pyscf.gto.Mole(atom="""
#     O    0.0000    0.0000    0.1173
#     H    0.0000    0.7572   -0.4692
#     H    0.0000   -0.7572   -0.4692 
# """,basis="6-311G")
# h2o = h2o.build()
# ```

# In[3]:


# h2o


# In[4]:


h2o = pyscf.gto.Mole(atom="""
    O    0.0000    0.0000    0.1173
    H    0.0000    0.7572   -0.4692
    H    0.0000   -0.7572   -0.4692 
""",basis="6-311G")
h2o = h2o.build()


# ## Aproximación Local de la Densidad (LDA)
# 
# Fue una de las primeras aproximaciones y ya ha sido superada. En este caso, el intercambio se calcula mediante
# 
# $$
# E_x^{LDA} = -\frac{3}{4} \left( \frac{3}{\pi} \right)^{1/3} \int \rho^{4/3} (r) dr 
# $$
# 
# y la correlación se calcula mediante el funcional de Vosko, Wilk y Nusair
# 
# $$
# E_c^{LDA} = \int \varepsilon_c^{VWN} dr
# $$
# 
# Haga un cálculo de energía con LDA y la base 6-311G con la siguiente instrucción
# ```
# rks = dft.RKS(h2o)
# rks.xc = "LDA,VWN"
# rks.kernel()
# ```

# In[5]:


# LDA


# In[6]:


rks = dft.RKS(h2o)
rks.xc = "LDA"
rks.kernel()


# ## Aproximación de Gradientes Generalizados (GGA)
# 
# Esto calcula la energía con base en la densidad electrónica y su gradiente
# 
# $$
# E_{xc}^{GGA} = -\int \varepsilon_{xc}^{GGA} (\rho,\nabla \rho) dr
# $$
# 
# Para ello separa los funcionales en intercambio y correlación.
# 
# $$
# \varepsilon_{xc}^{GGA} = \varepsilon_{x}^{GGA} + \varepsilon_{c}^{GGA}
# $$
# 
# Algunos funcionales de intercambio GGA son
# - PWx86: Perdew-Wang 1986
# - B88: Becke 1988
# - PWx91: Perdew-Wang 1991
# - PBE: Perdew-Burke-Ernzerhof
#     
# Algunos funcionales de correlación GGA son
# - LYP: Lee-Yang-Parr
# - Pc86: Perdew 1986
# - PWc91: Perdew-Wang 1991
# - PBE: Perdew-Burke-Ernzerhof
#     
# La combinación de estos funcionales genera los funcionales GGA. **Haga un cálculo de energía con PBE y la base 6-311G con la siguiente instrucción**
# ```
# rks = pyscf.dft.RKS(h2o)
# rks.xc = "LDA"
# rks.kernel()
# ```

# In[7]:


# PBE


# In[8]:


rks = dft.RKS(h2o)
rks.xc = "PBE"
rks.kernel()


# ## Aproximación meta-GGA
# 
# Los funcionales meta-GGA usan la densidad electrónica, el gradiente de la densidad electrónica, y el laplaciano de la densidad electrónica.
# 
# $$
# E_{xc}^{meta-GGA} = -\int \varepsilon_{xc}^{meta-GGA} (\rho,\nabla \rho,\nabla^2 \rho) dr
# $$
# 
# Algunos ejemplos son:
# - B95: Becke 1995
# - TPSS: Tau-Perdew-Staroverov-Scuseria
# 
# **Haga un cálculo de energía con TPSS y la base 6-311G con la siguiente instrucción**
# ```
# rks = pyscf.dft.RKS(h2o)
# rks.xc = "TPSS"
# rks.kernel()
# ```

# In[9]:


# TPSS


# In[10]:


rks = dft.RKS(h2o)
rks.xc = "TPSS"
rks.kernel()


# ## Funcionales Híbridos
# 
# Mezclan un funcional de intercambio con el `intercambio de Hartree-Fock` en alguna proporción.
# 
# Algunos ejemplos de estos funcionales son:
# 
# - B3LYP
# - PBE0
# - M05-2X y M06-2X
# - TPSSh
# 
# **Haga un cálculo de energía con M06-2X y la base 6-311G con la siguiente instrucción**
# ```
# rks = pyscf.dft.RKS(h2o)
# rks.xc = "M062X"
# rks.kernel()
# ```

# In[11]:


#M062X


# In[12]:


rks = dft.RKS(h2o)
rks.xc = "M062X"
rks.kernel()


# ## Referencias
# 
# - P. Hohenberg y W. Kohn, Inhomogeneous Electron Gas, Physical Review 136, B864 (1964).
# P. W. Atkins, y R. Friedman, Molecular Quantum Mechanics (Oxford University Press, 2005).
# - D. Rappoport, N. R. M. Crawford, F. Furche, y K. Burke, Which functional should I choose?, (2008).
# - W. Koch y M.C. Holthausen, A Chemist’s Guide to Density Functional Theory, (2001).
# - K. Burke y L.O. Wagner, DFT in a nutshell, Int. J. Quantum Chem. 113, 96 (2013).
# - K. Burke, Perspective on density functional theory, J. Chem. Phys. 136, (2012).
