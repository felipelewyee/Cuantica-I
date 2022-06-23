#!/usr/bin/env python
# coding: utf-8

# # Tutorial Básico

# A continuación aprenderá las funciones básicas de Python como declarar una variable, hacer operaciones básicas y utilizar librerías, así como lo necesario para seguir el manual. **Para usarlo es necesario escribir en los cuadros vacíos y presionar el botón Run o las teclas "Shift" + "Enter".**
# 
# ## Declarando variables
# 
# ````{admonition} Aprendizaje de código
# :class: important
# Al igual que en álgebra, en Python se pueden usar variables para **almacenar valores numéricos**. Por ejemplo:
# ```python
# x = 5.0
# ```
# ````
# 
# **En el siguiente recuadro declare una variable llamada "x" y asígnele el valor que desee, recuerde dar "shift"+"enter" o Run.**

# In[1]:


# Declare una variable x


# ````{admonition} Aprendizaje de código
# :class: important
# Es posible **ver (o imprimir) el valor de cualquier variable escribiendo**
# ~~~python
# print(variable)
# ~~~
# 
# En nuestro ejemplo
# ~~~python
# print(x)
# ~~~
# imprimiría el valor 5.0. 
# ````
# 
# **Imprima el valor de la variable que declaró anteriormente**

# In[2]:


#Imprima el valor de la variable x


# ## Operaciones Básicas

# ````{admonition} Aprendizaje de código
# :class: important
# 
# Las operaciones comunes del álgebra también pueden usarse en Python con una sintaxis muy similar.
# 
# | Operación     |    Álgebra    |    Python    |
# |:-------------:|:-------------:|:------------:|
# |Suma           |     $z=x+y$   |       $z=x+y$|
# |Resta          |     $z=x-y$   |       $z=x-y$|
# |Multiplicación |   $z=(x)(y)$  |       $x=x*y$|
# |División       |$z=\frac{x}{y}$|       $z=x/y$|
# |Potencia       |     $z=x^{y}$ |      $z=x**y$|
# 
# En nuestro ejemplo, podemos declarar una nueva variable, hacer las operaciones anteriores, guardar el resultado e imprimirlo.
# ~~~python
# #Declaramos la variable x
# x = 5.0
# #Declaramos la variable y
# y = 2.0
# 
# #Hacemos algunas operaciones con nuestras variables x, y; y guardamos el resultado
# z_suma  = x+y
# z_resta = x-y
# z_mult  = x*y
# z_div   = x/y
# z_pot   = x**y
# 
# #Imprimimos los resultados
# print("x",x,"y",y)
# print("suma",z_suma)
# print("resta",z_resta)
# print("multiplicación",z_mult)
# print("división",z_div)
# print("potencia",z_pot)
# ~~~
# 
# Los resultados se almacenan en las variables z_suma, z_resta, z_mult y z_div. Es importante notar que una variable puede tener un nombre más largo que una sola letra; y que también podemos imprimir texto.
# 
# ````
# 
# **Declare una variable llamada "y", asígnele el valor 10.0, calcule las 5 operaciones utilizando su variable ¨ "x" previamente definida, e imprima el resultado.**

# In[3]:


#Declare variable


# ## Librerías

# La funcionalidad de Python se puede extender utilizando **librerías**. Existen 4 importantes para el curso, **numpy**, **matplotlib**, **scipy** y **sympy**.
# 
# ````{admonition} Aprendizaje de código
# :class: important
# 
# Para activar una librería utilizamos la palabra **import**, seguida del nombre de la librería, por ejemplo
# ~~~python
# import scipy
# ~~~
# 
# Frecuentemente utilizamos un alias para asignarle un nombre más corto a la librería. Esto se realiza con la palabra **as** seguida del **alias**.
# ~~~python
# import numpy as np
# import sympy as sym
# ~~~
# 
# También podemos tomar partes específicas de una librería, para ello utilizamos la palabra **from**. Por ejemplo, importaremos el subpaquete pyplot de la librería matplotlib, la cual nos sirve para graficar, y lo renombraremos como plt.
# ~~~python
# from matplotlib import pyplot as plt
# ~~~
# 
# ````
# 
# **Copie y pegue las instrucciones anteriores en el recuadro para importar las librerías** numpy, matplotlib, scipy y sympy.

# In[4]:


#Importe librerías


# ````{admonition} Aprendizaje de código
# :class: important
# Utilizaremos numpy para obtener la raíz cuadrada de la variable x. Para ello llamaremos a la función sqrt de la librería numpy, esto se realiza escribiendo el nombre de la librería, seguido por un punto, seguido del nombre de la función a usar.
# ~~~python
# x = 5.0
# raiz = np.sqrt(x)
# print(raiz)
# ~~~
# ````
# 
# **Defina una variable x y obtenga e imprima su raíz cuadrada.**

# In[5]:


#Raíz cuadrada de x

