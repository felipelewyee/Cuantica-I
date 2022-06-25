#!/usr/bin/env python
# coding: utf-8

# # Rotor Rígido

# Este sistema consiste en dos partículas de masa $m_1$ y $m_2$ moviéndose con una separación constante $r=r_2-r_1$.
# 
# El Hamiltoniano consta de los términos de energía cinética de ambas partículas
# 
# $$
# H = -\frac{\hbar^2}{2m_1} \nabla^2_1 -\frac{\hbar^2}{2m_2} \nabla^2_2
# $$
# 
# ```{margin}
# La `masa reducida` tiene masa $\mu$ y coordenadas $r$, el `centro de masa` tiene masa $m_T$ y coordenadas $R_{\rm cm}$.
# ```
# 
# ```{admonition} Inserto matemático: Sistema de masa reducida
# :class: dropdown
# 
# Este sistema es equivalente al de una partícula de masa reducida ($\mu$) girando en torno al centro de masa de ambas partículas (${\rm cm}$). La masa total del sistema está dada por la suma de las masas de las partículas
# 
# $$
# m_T = m_1 + m_2
# $$
# 
# La masa reducida de la nueva partícula está dada por
# 
# $$
# \frac{1}{\mu} = \frac{1}{m_1} + \frac{1}{m_2}
# $$
# 
# El centro de masa del sistema se calcula mediante
# 
# $$
# R_{\rm cm} = \left( \frac{m_1}{m_T} \right) r_1 + \left( \frac{m_2}{m_T} \right) r_2
# $$
# 
# y el Hamiltoniano se calcula con la energía cinética de la masa reducida y del centro de masa, es decir
# 
# $$
# H = -\frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} - \frac{\hbar^2}{2\mu} \nabla^2_{\mu}
# $$
# 
# La ecuación de Schrödinger a resolver es
# 
# $$
# \left(-\frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} - \frac{\hbar^2}{2\mu} \nabla^2_{\mu}\right) \psi = E \psi
# $$
# 
# Se propone que la eigenfunción se puede separar en el producto de una eigenfunción del centro de masa y una eigenfunción de la partícula de masa reducida
# 
# $$
# \psi=\psi_{cm}\psi_{\mu}
# $$
# 
# Al sustituir en la ecuación de Schrödinger se obtiene
# 
# $$
# \left(-\frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} - \frac{\hbar^2}{2\mu} \nabla^2_{\mu}\right) \psi_{\rm cm}\psi_{\mu} = E \psi_{cm}\psi_{\mu}
# $$
# 
# Si consideramos que la energía está dada por $E_T = E_{\rm cm} + E_{r}$ y distribuimos, resulta
# 
# $$
# -\psi_{\mu} \frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} \psi_{\rm cm} - \psi_{\rm cm}\frac{\hbar^2}{2\mu} \nabla^2_{\mu} \psi_{\mu} = \psi_{\mu} E_{\rm cm} \psi_{\rm cm} +  \psi_{\rm cm} E_{\mu} \psi_{\mu}
# $$
# 
# Si multiplicamos ambos lados de la ecuación anterior por $\frac{1}{\psi_{\mu}\psi_{cm}}$, resulta
# 
# $$
# -\frac{1}{\psi_{\rm cm}} \left( \frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} \psi_{\rm cm} - E_{\rm cm} \psi_{\rm cm} \right) = \frac{1}{\psi_{\mu}} \left( \frac{\hbar^2}{2\mu} \nabla^2_{\mu} \psi_{\mu} + E_{\mu} \psi_{\mu} \right)
# $$
# 
# ya que el lado izquierdo solo depende de las coordenadas del centro de masa, y el lado derecho solo depende de las coordenadas de la masa reducida, ambos lados deben ser igual a una constante. Si elegimos esta constante como cero, y despejamos lo que está dentro de cada paréntesis se obtienen dos ecuaciones independientes.
# ```
# 
# Después de cambiar a un sistema de masa reducida se obtienen dos ecuaciones.
# 
# La primera ecuación corresponde al movimiento del `centro de masa del sistema` y la hemos estudiado previamente en el movimiento de la `partícula libre`
# 
# $$
# -\frac{\hbar^2}{2m_T} \nabla^2_{\rm cm} \psi_{\rm cm} = E_{\rm cm} \psi_{\rm cm}
# $$
# 
# esta ecuación tiene como soluciones
# 
# $$
# \psi_{\rm cm} = A e^{ikx} + B e^{-ikx}
# $$
# 
# con $k^2=2m_TE/\hbar^2$, y simplemente nos dice que el sistema en conjunto se mueve libremente por el espacio.
# 
# La segunda ecuación corresponde a la `masa reducida`, y la hemos estudiado en la `partícula en la esfera`
# 
# $$
# -\frac{\hbar^2}{2\mu} \nabla^2_{\mu} \psi_{\mu} = E_{\mu} \psi_{\mu}
# $$
# 
# sabemos por tanto que su solución son los armónico esféricos $Y_l^{m_l}(\theta,\phi)$ con $E = \frac{l(l+1)\hbar^2}{2\mu |r|^2}$.

# In[ ]:





# ## Referencias

# - C.C. Pye, On the Solution of the Quantum Rigid Rotor, J. Chem. Educ. 83, 460 (2006).
# - P. W. Atkins, y R. Friedman, Molecular Quantum Mechanics (Oxford University Press, 2005).
# - I.N. Levine, D.H. Busch, y H. Shull, Quantum chemistry (Pearson Prentice Hall Upper Saddle River, NJ, 2009).
# - D.A. McQuarrie y J.D. Simon, Physical Chemistry: A Molecular Approach (University Science Books, 1997).
