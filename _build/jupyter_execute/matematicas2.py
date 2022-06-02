#!/usr/bin/env python
# coding: utf-8

# # Álgebra simbólica

# Sympy es una librería de álgebra simbólica. Se importa así:

# In[1]:


import sympy as sym


# En vez de realizar operaciones de manera numérica, el álgebra simbólica guarda las operaciones como símbolos. Obtenga la siguiente operación:
# \begin{equation}
# \sqrt{11}
# \end{equation}
# Usando python estándar
# ~~~python
# import numpy as np
# a=np.sqrt(11)
# print(a)
# ~~~
# y usando sympy
# ~~~python
# import sympy as sym
# a=sympy.sqrt(11)
# print(a)
# ~~~

# In[2]:


import numpy as np
a=np.sqrt(11)
a


# In[3]:


import sympy as sym
a=sym.sqrt(11)
a


# Podemos mejorar la apariencia de una expresión si en vez de usar print utilizamos:
# ~~~python
# sym.pprint(a)
# ~~~
# o si simplemente poniendo la expresión.

# In[4]:


sym.pprint(a)


# Podemos obtener el valor numérico con evalf(). Obtenga el valor numérico de $\sqrt{11}$

# In[5]:


a.evalf()


# Las variables comunes de álgebra son símbolos. Si quiero definir x:

# In[6]:


x=sym.symbols("x")


# Imprimimos x para verificar.

# In[7]:


x


# Una variable puede almacenar cosas distintas a su nombre. Pruebe con
# ~~~python
# y=x+1
# ~~~

# In[8]:


y=x+1
y


# Es posible multiplicar expresiones y hacerles cualquier tipo de operación. Por ejemplo, haga el producto $xy$, guárdelo en una variable $z$ e imprima el resultado.

# In[9]:


z=y*x
z


# Se puede obtener el resultado de desarrollar las operaciones con 
# ~~~python
# sym.expand(expresión)
# ~~~
# Desarrolle la operación contenida en z, guárdelo en una variable c e imprímala.

# In[10]:


c=sym.expand(z)
sym.pprint(c)


# Podemos obtener expresiones en su forma factorizada con sp.simplify(). Pruebe con la variable c.
# 

# In[11]:


sym.simplify(c)


# Se pueden utilizar álgebra simbólica para derivar, esto se hace con 
# ~~~python
# sym.diff(funcion, variable a derivar1, variable a derivar2, variable a derivar3,...).
# ~~~
# 
# 
# Exprese de manera simbólica la función
# \begin{equation}
# f=x^2+sin(x)
# \end{equation}
# y obtenga $\frac{df}{dx}$ y $\frac{d^2f}{dx^2}$.

# In[12]:


x=sym.symbols("x")
f=x**2+sym.sin(x)
sym.pprint(f)
sym.pprint(sym.diff(f,x))
sym.pprint(sym.diff(f,x,x))


# Podemos integrar con 
# ~~~python
# sym.integrate(funcion, variable a integrar1, variable a integrar2, variable a integrar3,...)
# ~~~
# Defina la función
# \begin{equation}
# f=e^{-r}
# \end{equation}
# Y realice la integral indefinida
# \begin{equation}
# \int e^{-r} dr
# \end{equation}
# 

# In[13]:


r=sym.symbols("r")
f=sym.exp(-r)
sym.pprint(f)
sym.pprint(sym.integrate(f,r))


# Podemos integrar de manera definida con 
# ~~~python
# sym.integrate(funcion, (variable_a_integrar1,lim_inf1,lim_sup1), (variable_a_integrar2,lim_inf2,lim_sup2),(variable_a_integrar3,lim_inf3,lim_sup3),...)
# ~~~
# Defina la función
# \begin{equation}
# f=e^{-\alpha r^2}
# \end{equation}
# Y realice la integral
# \begin{equation}
# \int_0^\infty e^{-\alpha r^2} dr
# \end{equation}
# 

# In[14]:


r=sym.symbols("r")
alpha=sym.symbols("alpha",positive=True)
f=sym.exp(-alpha*r**2)
sym.pprint(sym.integrate(f,(r,0,sym.oo)))


# También se pueden resolver ecuaciones diferenciales. dsolve indica resolver una ecuación, y Eq indica la ecuación.
# Resuelva
# \begin{equation}
# -\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi(x) = E \psi(x)
# \end{equation}
# Con $k^2=\frac{2mE}{\hbar^2}$
# \begin{equation}
# \frac{d^2}{dx^2} \psi(x) + k^2 \psi(x) = 0
# \end{equation}

# In[15]:


x=sym.symbols("x")
k=sym.symbols("k")

psi=sym.Function("psi")
eq=sym.Eq(psi(x).diff(x,x)+k**2*psi(x),0)
sym.pprint(eq)

sym.pprint(sym.dsolve(eq,psi(x)))


# La parte radial del Hamiltoniano del átomo de hidrógeno es
# \begin{equation}
# H = -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r}
# \end{equation}
# 
# Se usa la función de prueba
# \begin{equation}
# \psi^{prueba} = \left( \frac{2\alpha}{\pi} \right)^{3/4} e^{-\alpha r^2}
# \end{equation}
# 
# La energía se obtiene como:
# \begin{equation}
# E^{prueba} = \int_{0}^{r=\infty} \int_{0}^{\theta=\pi} \int_{0}^{\phi=2\pi} r^2 sin \theta \left( {\psi^{prueba}}^* \hat{H} \psi^{prueba} \right) dr d\theta d\phi
# \end{equation}
# 

# In[16]:


r=sym.symbols("r")
alpha=sym.symbols("alpha",positive=True)
psi=(2*alpha/sym.pi)**(sym.S(3)/4)*sym.exp(-alpha*r**2)
#sym.pprint(psi)

sym.pprint(sym.integrate(4*sym.pi*r**2*psi*(-1/2*1/r**2*sym.diff(r**2*sym.diff(psi,r),r)-psi/r),(r,0,sym.oo)))


# ## Átomo de hidrógeno

# La parte radial del Hamiltoniano del átomo de hidrógeno es
# \begin{equation}
# H = -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r}
# \end{equation}
# 
# Se usa la función de prueba
# \begin{equation}
# \psi^{prueba} = \left( \frac{2\alpha}{\pi} \right)^{3/4} e^{-\alpha r^2}
# \end{equation}
# 
# La energía se obtiene como:
# \begin{eqnarray}
# E^{prueba} &=& \int_{0}^{r=\infty} \int_{0}^{\theta=\pi} \int_{0}^{\phi=2\pi} r^2 sin \theta \left( {\psi^{prueba}}^* \hat{H} \psi^{prueba} \right) dr d\theta d\phi
# \end{eqnarray}
# 

# **Paso 1.** Identifique las variables y declárelas como símbolos. En este caso: $r$ y $alpha$.

# In[17]:


r=sym.symbols("r")
alpha=sym.symbols("alpha",positive="True") #positive=True indica que alpha solo puede ser positivo.


# **Paso 2.** Identifique si existe alguna función que pueda ser expresada con las variables anteriores. En este caso sí la hay y es $\psi^{prueba}$.

# In[18]:


psi=(2*alpha/sym.pi)**(sym.S(3)/4)*sym.exp(-alpha*r**2)
sym.pprint(psi)


# **Paso 3.** Opcional. Identifique si hay partes de la ecuación que pueda agrupar o distribuir con operaciones simples, como multiplicación y expréselas por separado.
# 
# En este caso
# \begin{eqnarray}
# E^{prueba} &=& \int_{0}^{r=\infty} \int_{0}^{\theta=\pi} \int_{0}^{\phi=2\pi} r^2 sin \theta \left( {\psi^{prueba}}^* \hat{H} \psi^{prueba} \right) dr d\theta d\phi\\
# &=& \int_{0}^{r=\infty} \int_{0}^{\theta=\pi} \int_{0}^{\phi=2\pi} r^2 sin \theta \left( {\psi^{prueba}}^* \left( -\frac{1}{2} \left( \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \right)-\frac{1}{r} \right) \psi^{prueba} \right) dr d\theta d\phi\\
# &=& 4\pi \int_{0}^{r=\infty} r^2 \left( {\psi^{prueba}}^* \left( \color{red}{ -\frac{1}{2}  \frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} \psi^{prueba}} \color{blue}{-\frac{1 }{r}\psi^{prueba}} \right)  \right) dr \\
# \end{eqnarray}
# 
# A la parte en $\color{red}{rojo}$ le llamaremos A, y a la parte en $\color{blue}{azul}$ le llamaremos B.

# In[19]:


# Primera parte (en rojo)
d_psi=sym.diff(psi,r)
dd_psi=sym.diff(r**2*d_psi,r)
A=(-sym.S(1)/2)*(1/r**2)*dd_psi

# Primera parte (en azul)
B=-(1/r)*psi


# **Paso 4.** Exprese la función a integrar y haga la integral con $sp.integrate()$

# In[20]:


f=r**2*psi*(A+B)
E=sym.integrate(4*sym.pi*f,(r,0,sym.oo))
sym.pprint(E)


# ## Referencias

# Página oficial de sympy
# www.sympy.org
