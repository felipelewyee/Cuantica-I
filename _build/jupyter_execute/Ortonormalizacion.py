#!/usr/bin/env python
# coding: utf-8

# # Normalización y Ortogonalidad

# ## Normalización

# Frecuentemente se busca que la función de onda esté normalizada, es decir, que se cumpla la siguiente integral
# 
# $$
# \int \psi_{normalizada} \psi_{normalizada}* = \int |\psi_{normalizada}|^2 = 1
# $$
# 
# Las funciones de onda no siempre están normalizadas, pero podemos realizar los siguientes pasos para construir una `función de onda normalizada`.
# 
# ```{tip}
# 1 Integrar el cuadrado de la función de onda original en el dominio.
# 
# $$
# N^2 = \int_{x_1}^{x_2} \psi_{original}^2 dr
# $$
# 
# 2 Obtener la norma de la función de onda.
# 
# $$
# N = \sqrt{N^2} = \sqrt{\int_{x_1}^{x_2} |\psi_{original}|^2 dx}
# $$
# 
# 3 Multiplicar la función de onda original por el inverso de su norma.
# 
# $$
# \psi_{normalizada} = \frac{1}{N} \psi_{original}
# $$
# ```
# 
# ```{margin}
# El inverso de la norma también recibe el nombre de constante de normalización, $C$, es decir, $C = \frac{1}{N}$. En este caso, la función de onda normalizada se obtiene mediante $\psi_{normalizada} = C \psi_{original}$.
# ```
# 
# ```{margin}
# Si $N=1$, la función de onda original ya está normalizada.
# ```

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# También podemos hacer integrales con Python. Por ejemplo, vamos a integrar $y=x^2$ en el dominio $-3 \leq x \leq3$. Para ello importaremos el subpaquete integrate de la librería scipy.
# ~~~python
# from scipy import integrate
# ~~~
# y luego realizaremos la integral con la función quad
# ~~~python
# integrate.quad(lambda x: x**2,-3,3)
# ~~~
# Aquí "lambda" indica las variables de la ecuación, seguido de la ecuación y los límites de la integral. **Pruebe a realizar la integral.**

# In[1]:


# Integre y = x**2 de -3 a 3


# In[2]:


import scipy.integrate as integrate
integrate.quad(lambda x: x**2,-3,3)


# Considere la función de onda $\psi = x$, con $x \epsilon [-3,3]$, **proponga una función de onda normalizada y evalúe la integral con integrate.quad para comprobar que la norma es 1.**
# 
# ```{tip}
# Recuerde que para normalizar:
# 1 Integre el cuadrado de la función de onda en el dominio.
# 
# $$
# N^2 = \int_{x_1}^{x_2} \psi_original^2 dr
# $$
# 
# 2 Obtenga la norma.
# 
# $$
# N = \sqrt{N^2} = \sqrt{\int_{x_1}^{x_2} |\psi_{original}|^2 dx}
# $$
# 
# 3 Multiplique la función de onda original por el inverso de su norma.
# 
# $$
# \psi_{normalizada} = \frac{1}{N} \psi_{original}
# $$
# ```

# In[3]:


#Normalice función de onda


# In[4]:


import numpy as np

norm2 = integrate.quad(lambda x: x*x,-3,3)[0]
print("N**2 =",norm2)

norm = np.sqrt(norm2)
print("N =",norm)

integrate.quad(lambda x: (x/norm)*(x/norm),-3,3)


# **OPCIONAL**
# 
# También es posible realizar integrales triples, el siguiente ejemplo solo es demostrativo, **simplemente copie y pegue para realizar la integral.**
# 
# Sea la función de onda
# 
# $$
# \psi = e^{-(x^2+y^2+z^2)}=e^{-r^2}
# $$
# 
# La integral de su cuadrado será
# 
# $$
# N^2 = \int \psi^* \psi d\textbf{r} = \int\limits_{0}^\infty \int\limits_{0}^{2\pi} \int\limits_{0}^\pi\psi^* \psi r^2 sin\theta dr d\phi d\theta = \int\limits_{0}^{\pi} sin \theta d\theta \int\limits_{0}^{2\pi} d\phi \int\limits_{0}^\infty e^{-2r^2} r^2 dr = \left(2\pi\right)\left(2\right)\left(\frac{1}{8}\sqrt{\frac{\pi}{2}}\right) = \left(\frac{\pi}{2}\right)^{3/2}
# $$
# 
# El siguiente código evalúa la integral y la guarda en la variable "norm2".
# ~~~python
# norm2 = integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
# ~~~

# In[5]:


#Encuentre el cuadrado de la norma de la funión de onda


# In[6]:


norm2 = integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]

print("N**2 =",norm2)


# **Compruebe que numéricamente esto es $\left(\frac{\pi}{2}\right)^{3/2}$.**

# In[7]:


#Evalúe (pi/2)^(3/2)


# In[8]:


print((np.pi/2)**(3/2))


# Evalúe la norma

# In[9]:


# Constante de normalización


# In[10]:


norm = np.sqrt(norm2)
print("N =",norm)


# **OPCIONAL**
# 
# Proponga una $\psi$ normalizada, y compruebe que la integral de su cuadrado da 1. Recuerde la siguiente relación
# 
# $$
# \psi_{normalizada} = \frac{1}{\sqrt{N}} \psi_{original}
# $$
# 
# y que para este ejercicio
# 
# $$
# \psi_{original} = e^{-r}
# $$

# In[11]:


#Comprobación


# In[12]:


norm = integrate.tplquad(lambda theta,phi,r: 1/norm**2*r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(norm)


# **OPCIONAL**
# 
# **Repita el proceso de normalización para $\psi=re^{-r^2}$.** La integral sin normalizar será
# 
# $$
# \int \psi^* \psi d\textbf{r} = \frac{3}{4}\left(\frac{\pi}{2}\right)^{3/2}
# $$

# In[13]:


#Normalice


# In[14]:


norm2=integrate.tplquad(lambda theta,phi,r: r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
print("N**2 =",norm2)

norm = np.sqrt(norm2)
print("N =",norm)

norm=integrate.tplquad(lambda theta,phi,r: 1/norm**2*r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(norm)

