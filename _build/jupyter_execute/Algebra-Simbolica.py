#!/usr/bin/env python
# coding: utf-8

# # Álgebra simbólica

# `Sympy` es una librería de álgebra simbólica. Se importa así:

# In[1]:


import sympy as sym


# En vez de realizar operaciones de manera numérica, el álgebra simbólica guarda las operaciones como símbolos. Obtenga la siguiente operación:
# \begin{equation*}
# \sqrt{11}
# \end{equation*}
# 
# ````{admonition} Aprendizaje de código
# :class: importantUsando python estándar
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
# ````

# In[2]:


# raiz de 11 en Python estándar


# In[3]:


import numpy as np
a=np.sqrt(11)
print(a)


# In[4]:


# raíz de 11 en Python simbólico


# In[5]:


import sympy as sym
a=sym.sqrt(11)
print(a)


# ```{admonition} Para pensar
# ¿Qué diferencias observa entre los resultados de `print(a)` de python estándar y de python simbólico?
# ```

# ````{admonition} Aprendizaje de código
# :class: important
# Podemos mejorar la apariencia de una expresión si en vez de usar print utilizamos:
# ~~~python
# sym.pprint(a)
# ~~~
# 
# o simplemente escribimos el nombre de la variable.
# ````

# In[6]:


# Pruebe la alternativa de impresión.


# In[7]:


sym.pprint(a)


# ````{admonition} Aprendizaje de código
# :class: important
# Podemos obtener el valor numérico con evalf().
# ````
# **Obtenga el valor numérico de $\sqrt{11}$**

# In[8]:


# Valor numérico de raíz de 11


# In[9]:


a.evalf()


# ````{admonition} Aprendizaje de código
# :class: important
# Las variables comunes de álgebra son símbolos. Si queremos definir x, escribimos
# 
# ~~~python
# x=sym.symbols("x")
# ~~~
# ````

# In[10]:


# Defina la variable simbólica x


# In[11]:


x=sym.symbols("x")


# **Imprimima x para verificar.**

# In[12]:


# Imprima x


# In[13]:


x


# Una variable puede almacenar cosas distintas a su nombre. Pruebe con
# ````{admonition} Aprendizaje de código
# :class: important
# ~~~python
# y=x+1
# ~~~
# ````
# e imprima $y$. Observe que el resultado no es un número, sino una expresión.

# In[14]:


# Defina "y" simbólica e imprímala.


# In[15]:


y=x+1
y


# Es posible multiplicar expresiones y hacerles cualquier tipo de operación. Por ejemplo, haga el producto $xy$, guárdelo en una variable $z$ e imprima el resultado.

# In[16]:


# z = xy


# In[17]:


z = x*y
z


# ````{admonition} Aprendiendo código
# :class: important
# Se puede obtener el resultado de desarrollar las operaciones con 
# ~~~python
# sym.expand(expresión)
# ~~~
# ````
# 
# Desarrolle la operación contenida en z, guárdelo en una variable c e imprímala.

# In[18]:


# Expanda el contenido de la variable z


# In[19]:


c=sym.expand(z)
sym.pprint(c)


# ````{admonition} Aprendizaje de codigo
# :class: important
# Podemos obtener expresiones en su forma factorizada con sp.simplify().
# ````
# 
# Pruebe con la variable c.
# 

# In[20]:


# Simplifique c


# In[21]:


sym.simplify(c)


# ````{admonition} Aprendizaje de código
# :class: important
# Se pueden utilizar álgebra simbólica para derivar, esto se hace con 
# ~~~python
# sym.diff(funcion, variable a derivar, orden de la derivada).
# ~~~
# o si vamos a derivar respecto a diferentes variables
# ~~~python
# sym.diff(funcion, variable a derivar1, variable a derivar2, variable a derivar3,...).
# ~~~
# ````
# 
# Exprese de manera simbólica la función
# \begin{equation}
# f=x^2+sin(x) \nonumber
# \end{equation}
# y obtenga $\frac{df}{dx}$ y $\frac{d^2f}{dx^2}$.

# In[22]:


# Exprese f(x) y derive a primer y segundo orden.


# In[23]:


x=sym.symbols("x")
f=x**2+sym.sin(x)
sym.pprint(f)
sym.pprint(sym.diff(f,x))
sym.pprint(sym.diff(f,x,2))


# ````{admonition} Aprendizaje de código
# :class: important
# Podemos integrar simbólicamente con 
# ~~~python
# sym.integrate(funcion, variable a integrar1, variable a integrar2, variable a integrar3,...)
# ~~~
# Defina la función
# \begin{equation*}
# f=e^{-r}
# \end{equation*}
# Y realice la integral indefinida
# \begin{equation*}
# \int e^{-r} dr
# \end{equation*}
# ````

# In[24]:


# Defina f(r) = e^{-r} y realice la integral indicada.


# In[25]:


r=sym.symbols("r")
f=sym.exp(-r)
sym.pprint(f)
sym.pprint(sym.integrate(f,r))


# ````{admonition} Aprendizaje de código
# :class: important
# Podemos integrar de manera definida con
# ~~~python
# sym.integrate(funcion, (variable_a_integrar1,lim_inf1,lim_sup1), (variable_a_integrar2,lim_inf2,lim_sup2),(variable_a_integrar3,lim_inf3,lim_sup3),...)
# ~~~
# ````
# Defina la función
# \begin{equation*}
# f=e^{-\alpha r^2}
# \end{equation*}
# Y realice la integral
# \begin{equation*}
# \int_0^\infty e^{-\alpha r^2} dr
# \end{equation*}

# In[26]:


# Realice integral definida


# In[27]:


r=sym.symbols("r")
alpha=sym.symbols("alpha",positive=True)
f=sym.exp(-alpha*r**2)
sym.pprint(sym.integrate(f,(r,0,sym.oo)))


# También se pueden resolver ecuaciones diferenciales. dsolve indica resolver una ecuación, y Eq indica la ecuación.
# Resuelva
# \begin{equation*}
# -\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi(x) = E \psi(x)\,,
# \end{equation*}
# definiendo $k^2=\frac{2mE}{\hbar^2}$
# \begin{equation*}
# \frac{d^2}{dx^2} \psi(x) + k^2 \psi(x) = 0\,.
# \end{equation*}

# In[28]:


x=sym.symbols("x")
k=sym.symbols("k")

psi=sym.Function("psi")
eq=sym.Eq(psi(x).diff(x,x)+k**2*psi(x),0)
sym.pprint(eq)

sym.pprint(sym.dsolve(eq,psi(x)))


# ## Referencias

# Página oficial de sympy
# www.sympy.org
