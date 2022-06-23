#!/usr/bin/env python
# coding: utf-8

# # Método Variacional

# ````{admonition} Ver Antes
# :class: seealso
# 
# Antes de resolver estos ejercicios, puede resultarle de utilidad revisar el Notebook de [](Algebra-Simbolica.ipynb).
# ````

# ```{margin}
# La expresión del valor esperado de la energía mostrada considera el sistema de una partícula en tres dimensiones.
# ```
# 
# ````{admonition} Principio Variacional
# :class: Note
# 
# Si se usa una función de prueba normalizada que satisfaga las condiciones de frontera del sistema, el valor esperado del Hamiltoniano es un límite superior a la energía, es decir,
# 
# \begin{eqnarray*}
#     E_{\rm prueba} &=& \int_{0}^{\phi=2\pi} \int_{0}^{\theta=\pi} \int_{0}^{r=\infty} r^2 \sin \theta \left( \psi_{\rm prueba}^* \hat{H} \psi_{\rm prueba} \right) dr d\theta d\phi \geq E_{\rm exacta} \>\>.
# \end{eqnarray*}
# 
# Esto significa que se puede proponer una función de prueba que satisfaga las condiciones del sistema, y el valor esperado de la energía siempre estará por arriba del valor exacto, y sólo será igual si la función de prueba corresponde a la función de onda exacta del sistema.
# ````

# ## Átomo de hidrógeno

# La parte radial del Hamiltoniano del átomo de hidrógeno es
# \begin{equation*}
# H = -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r} \>\>.
# \end{equation*}
# 
# ```{margin}
# El Hamiltoniano del átomo de hidrógeno contiene una parte radial y una parte angular, sin embargo, estas son separables y la energía sólo depende de la parte radial.
# ```
# 
# Considere la función de prueba con $\alpha>0$
# \begin{equation*}
# \psi_{\rm prueba} = \left( \frac{2\alpha}{\pi} \right)^{3/4} e^{-\alpha r^2}  \>\>,
# \end{equation*}
# ```{margin}
# Note que la función de prueba propuesta no es la eigenfunción del estado base.
# ```
# 
# calcule la energía correspondiente al átomo de hidrógeno usando la función de prueba,
# \begin{eqnarray*}
# E_{\rm prueba} &=& \int_{0}^{\phi=2\pi} \int_{0}^{\theta=\pi} \int_{0}^{r=\infty} r^2 \sin \theta \left( \psi_{\rm prueba}^* \hat{H} \psi_{\rm prueba} \right) dr d\theta d\phi \>\>.
# \end{eqnarray*}

# **Paso 0.** Importe sympy.

# In[1]:


# Importe sympy


# In[2]:


import sympy as sym


# **Paso 1.** Identifique las variables y declárelas como símbolos. En este caso: $r$ y $alpha$.

# In[3]:


# Declare las variables r y alpha como símbolos


# In[4]:


r=sym.symbols("r")
alpha=sym.symbols("alpha",positive="True") #positive=True indica que alpha solo puede ser positivo.


# **Paso 2.** Declare la función de prueba, $\psi_{\rm prueba}$, con las variables anteriores.

# In[5]:


# Declare la función de prueba


# In[6]:


psi=(2*alpha/sym.pi)**(sym.S(3)/4)*sym.exp(-alpha*r**2)
sym.pprint(psi)


# **Paso 3.** Opcional. Identifique si hay partes de la ecuación que pueda agrupar o distribuir con operaciones simples, como multiplicación y expréselas por separado.
# 
# En este caso
# \begin{eqnarray*}
# E_{\rm prueba} &=& \int_{0}^{\phi=2\pi} \int_{0}^{\theta=\pi} \int_{0}^{r=\infty} r^2 \sin \theta \bigg( \psi_{\rm prueba}^* \hat{H} \psi_{\rm prueba} \bigg) dr d\theta d\phi\\
# &=& \int_{0}^{\phi=2\pi} \int_{0}^{\theta=\pi} \int_{0}^{r=\infty} r^2 \sin \theta \left\{ \psi_{\rm prueba}^* \left[ -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r} \right] \psi_{\rm prueba} \right\} dr d\theta d\phi\\
# &=& \int_{0}^{\theta=\pi} \sin \theta d\theta \int_{0}^{\phi=2\pi} d\phi \int_{0}^{r=\infty} r^2 \left\{ \psi_{\rm prueba}^* \left[ -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r} \right] \psi_{\rm prueba} \right\} dr\\
# &=& 4\pi \int_{0}^{r=\infty} r^2 \left[ \psi_{\rm prueba}^* \left( \color{red}{ -\frac{1}{2}  \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \psi_{\rm prueba}} \color{blue}{-\frac{1 }{r}\psi_{\rm prueba}} \right)  \right] dr \\
# \end{eqnarray*}
# 
# Sugerimos que a la parte en $\color{red}{rojo}$ le llamemos A, y a la parte en $\color{blue}{azul}$ le llamemos B. Esto es opcional pero hará más fácil escribir la integral.
# 
# ```{margin}
# Recuerde que puede hacer derivadas con la función `sym.diff()`
# ```

# In[7]:


# Declare A (rojo) y B (azul)


# In[8]:


# Primera parte (en rojo)
d_psi=sym.diff(psi,r)
dd_psi=sym.diff(r**2*d_psi,r)
A=(-sym.S(1)/2)*(1/r**2)*dd_psi

# Primera parte (en azul)
B=-(1/r)*psi


# **Paso 4.** Exprese la función a integrar y resuelva la integral con `sp.integrate()`

# In[9]:


# Declare la función e integre para obtener la energía


# In[10]:


f=r**2*psi*(A+B)
E=sym.integrate(4*sym.pi*f,(r,0,sym.oo))
sym.pprint(E)


# In[11]:


from OptMultiple import MultipleChoice


# In[12]:


question = "¿Cómo difiere la energía obtenida con la función de prueba respecto a la energía exacta del átomo de hidrógeno?"
answers = [
    "La energía de la función de prueba es mayor que la exacta.",
    "La energía de la función de prueba es igual a la exacta.",
    "La energía de la función de prueba es menor que la exacta.",
]
explanation = (
    "Dado que la función de prueba no es la exacta, por el principio variacional se cumple que \( E_{\mathrm{prueba}} > E_{\mathrm{exacta}} \)." 
    "Note que para cualquier \(\\alpha>0\), la energía de prueba obtenida es más positiva que la energía exacta del átomo de hidrógeno de \(-0.5\) Hartree."
)
MultipleChoice(
    question, answers, correct_answer=0, explanation=explanation
)


# ```{admonition} Pregunta
# :class: note
# 
# ¿Cómo podría obtener el valor de $\alpha$ que de la energía más cercana a la exacta para la función de prueba propuesta?
# 
# **Idea.** Recuerde máximos y mínimos de Cálculo I.
# 
# ```

# ## Referencias

# - F.L. Pilar, Elementary Quantum Chemistry (Dover ed., 2001).
# -  A. Szabo y N.S. Ostlund, Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory (Dover ed., 1996).
