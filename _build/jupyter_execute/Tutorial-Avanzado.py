#!/usr/bin/env python
# coding: utf-8

# # Tutorial Avanzado

# En esta parte del tutorial aprenderá a utilizar Python para manejar vectores y matrices, crear gráficas y resolver integrales. Estas funciones serán trascendentales para el curso de `Química Cuántica`. **Recuerde escribir en los cuadros vacíos y presionar el botón Run o las teclas "Shift" + "Enter".**
# 
# Utilizaremos las librerías numpy, matplotlib y scipy que aprendió a utilizar en el Tutorial Básico, **impórtelas en el siguiente recuadro**.

# In[1]:


# Importe librerías


# In[2]:


import numpy as np
from matplotlib import pyplot as plt
import scipy


# ## Vectores y Matrices

# Python permite manejar fácilmente objetos multidimensionales, como vectores y matrices.
# ````{admonition} Aprendizaje de código
# :class: important
# Empezaremos creando una lista, la cual es un conjunto de cosas. Las listas se crean poniendo elementos entre corchetes cuadrados, por ejemplo
# ~~~python
# lista = [5, 9, 4]
# ~~~
# 
# Podemos imprimir todos los elementos de la lista con
# ~~~python
# print(lista)
# ~~~
# También podemos imprimir un elemento específico indicando el número del elemento entre paréntesis cuadrados. En Python los elementos se cuentan desde "cero". Para imprimir el número 9 en la lista anterior ejecutaríamos:
# ~~~python
# print(lista[1])
# ~~~
# ````
# 
# **Cree la siguiente lista [8,25,32,41] y guárdela en la variable "lista". Imprima la lista completa, así como el primer elemento del vector.**

# In[3]:


#Cree una lista


# In[4]:


lista = [8,25,32,41]
print(lista)
print(lista[0])


# La lista puede convertirse en un arreglo de numpy, esto facilita el manejo de cantidades vectoriales.
# ````{admonition} Aprendizaje de código
# :class: important
# Para declarar la lista anterior como un arreglo de numpy escribiremos
# ~~~python
# vector = np.array([5, 9, 4])
# ~~~
# 
# Una ventaja es que podemos declarar dos arreglos, por ejemplo
# ~~~python
# v1 = np.array([2, 4, 6])
# v2 = np.array([3, 5, 7])
# ~~~
# 
# y sumarlos, de una forma muy parecida a como se haría en física.
# 
# ~~~python
# print(v1+v2)
# ~~~
# ````
# 
# **Cree dos vectores como arreglos de numpy, uno que contenga los elementos "1,8,6", y otro que contenga "5,2,7", e imprima su suma**

# In[5]:


#Suma de dos vectores


# In[6]:


v1 = np.array([1,8,6])
v2 = np.array([5,2,7])

print(v1+v2)


# ````{admonition} Aprendizaje de código
# :class: important
# También podemos crear matrices como una "lista de listas", donde cada lista representa un renglón de la matriz. A continuación generemos una matriz de 3x3 que contenga ceros en todos sus elementos
# ~~~python
# matriz = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
# ~~~
# ````
# **Genere la matriz de ceros anterior, e imprímala.**

# In[7]:


#Cree matriz


# ````{admonition} Aprendizaje de código
# :class: important
# Podemos usar numpy para hacer nuestra matriz de ceros de forma más fácil
# ~~~python
# matriz = np.zeros((3,3))
# print(matriz)
# ~~~
# ````
# **Prueba usar numpy para definir una matriz de ceros de 5x5, e imprímala.**

# In[8]:


#Cree matriz de ceros con numpy


# In[9]:


matriz = np.zeros((5,5))
print(matriz)


# **A continuación genere una matriz de 3x3, guárdela en la variable "matriz1", y asígnele los valores del 1 al 9, tal que**
# 
# $$
# \text{matriz}1 = \begin{pmatrix}
# 1 & 2 & 3\\
# 4 & 5 & 6\\
# 7 & 8 & 9\\
# \end{pmatrix}
# $$
# 
# Puede ayudarse del siguiente código.
# ~~~python
# matriz1 = np.zeros((3,3))
# matriz1[0][0]=1.0
# matriz1[0][1]=2.0
# matriz1[0][2]=3.0
# matriz1[1][0]=4.0
# matriz1[1][2]=5.0
# ...
# ~~~

# In[10]:


#Llene una matriz de 3x3


# In[11]:


matriz1 = np.zeros((3,3))
matriz1[0][0]=1.0
matriz1[0][1]=2.0
matriz1[0][2]=3.0
matriz1[1][0]=4.0
matriz1[1][1]=5.0
matriz1[1][2]=6.0
matriz1[2][0]=7.0
matriz1[2][1]=8.0
matriz1[2][2]=9.0


# ````{admonition} Aprendizaje de código
# :class: important
# Hacer esto puede tomar mucho tiempo con matrices grandes. Existe una instrucción llamada **for** que nos permite hacer cosas repetitivas. Por ejemplo:
# ~~~python
# matriz2 = np.zeros((3,3))
# val=0.0
# for i in range(3):
#     for j in range(3):
#         val=val+1.0
#         matriz2[i][j]=val
# ~~~
# ````
# **Prueba la instrucción anterior e imprime las dos matrices (matriz1 y matriz2) en el siguiente recuadro.**

# In[12]:


#Matriz de 3x3 con for


# ````{admonition} Aprendizaje de código
# :class: important
# Se pueden multiplicar dos matrices con
# ~~~python
# matriz3=np.matmul(matriz1,matriz2)
# ~~~
# ````
# **Multiplique las dos matrices anteriores e imprima el resultado.**

# In[13]:


#Matriz3 = Matriz1 x Matriz2


# Los valores propios ($\lambda$) y vectores propios ($v$) de una matriz ($M$) son muy importantes en química cuántica. Cumplen la propiedad de que al multiplicar la matriz por el vector propio resulta el mismo vector multiplicado por una constante, es decir:
# 
# $$
# M v = \lambda v
# $$
# 
# ````{admonition} Aprendizaje de código
# :class: important
# Los valores y vectores propios se pueden calcular con
# ~~~python
# val,vec=np.linalg.eig(matriz)
# ~~~
# ````
# 
# **Calcule los valores propios y vectores propios de la matriz 1, e imprímalos.**
# 
# ```{margin}
# Si la matriz es simétrica o hermitiana es más conveniente usar `np.linalg.eigh(matriz)`.
# ```
# 
# ```{margin}
# `np.linalg.eig(matriz)` regresa valores propios desordenados mientras que `np.linalg.eigh(matriz)` los regresa ordenados.
# ```

# In[14]:


# Valores y vectores propios


# In[15]:


val,vec=np.linalg.eig(matriz1)


# ## Gráficas

# Vamos a graficar la función $y=\sin(x)$ de $-3$ a $3$.
# 
# ````{admonition} Aprendizaje de código
# :class: important
# Primero crearemos el dominio de la función (los valores de x), en nuestro ejemplo le daremos 50 puntos. Utilizaremos la función linspace(a,b,n) de numpy, esta crea un conjunto de n números distribuidos en el intervalo de a hasta b.
# ~~~python
# x=np.linspace(-3,3,50)
# ~~~
# ````
# 
# **Genere el dominio tal que $x\in[-3,3]$, utilice 50 puntos, guárdelo en la variable x e imprímalo.**

# In[16]:


# Cree dominio de x


# ````{admonition} Aprendizaje de código
# :class: important
# Obtendremos el valor de "y" usando numpy. Al escribir la siguiente instrucción Python toma cada valor de x, le aplica la función y lo guarda en la variable "y".
# ~~~python
# y=np.sin(x)
# ~~~
# ````
# 
# **Evalúe la función**
# 
# $$
# y = \sin(x)
# $$

# In[17]:


# Evalúe función


# ````{admonition} Aprendizaje de código
# :class: important
# Finalmente, **realice la gráfica con**
# 
# ~~~python
# plt.scatter(x,y)
# plt.show()
# ~~~
# ````
# 
# ```{margin}
# Prueba cambiar `plt.scatter(x,y)` por `plt.plot(x,y)` para tener una gráfica continua.
# ```

# In[18]:


# Genere gráfica


# In[19]:


from matplotlib import pyplot as plt
import numpy as np

x=np.linspace(-3,3,50)
y=np.sin(x)

plt.scatter(x,y)
#plt.plot(x,y)
plt.show()


# Las instrucciones anteriores son los pasos básicos para generar una gráfica. **En el siguiente recuadro genere la gráfica de las funciones $y_1=e^{-|x|}$ y $y_2=e^{-|x|^2}$ con 100 puntos desde -3 hasta 3.**

# In[20]:


# Gráfica


# In[21]:


x  = np.linspace(-3,3,100)
y1 = np.exp(-np.abs(x))
y2 = np.exp(-np.abs(x)**2)

plt.scatter(x,y1)
plt.scatter(x,y2)

plt.show()


# ## Integrales

# ````{admonition} Aprendizaje de código
# :class: important
# También podemos hacer integrales con Python. Por ejemplo, vamos a integrar $y=x^2$ en el dominio $-3 \leq x \leq 3$. Para ello importaremos el subpaquete integrate de la librería scipy.
# ~~~python
# from scipy import integrate
# ~~~
# y luego realizaremos la integral con la función quad
# 
# ~~~python
# integrate.quad(lambda x: x**2,-3,3)[0]
# ~~~
# Aquí "lambda" indica las variables de la ecuación, seguido de la ecuación y los límites de la integral.
# ````
# 
# ```{margin}
# El resultado de `integrate.quad` es un vector, por eso al final de la línea escribimos `[0]` para acceder al primer elemento, que contiene el resultado de la integral. Pruebe escribiendo `integrate.quad(lambda x: x**2,-3,3)` y observe la diferencia.
# ```
# 
# **Pruebe a realizar la integral.**

# In[22]:


# Integre y = x**2 de -3 a 3


# In[23]:


import scipy.integrate as integrate
integrate.quad(lambda x: x**2,-3,3)[0]


# Considere la función de onda $\psi = x$, con $x \epsilon [-3,3]$, **proponga una función de onda normalizada y evalúe la integral con integrate.quad para comprobar que la norma es 1.**
# 
# ```{tip}
# Recuerde que para normalizar:
# 
# 1 Integre el cuadrado de la función de onda en el dominio.
# 
# $$
# N^2 = \int_{x_1}^{x_2} |\psi_ {\rm original}|^2 dr
# $$
# 
# 2 Obtenga la norma.
# 
# $$
# N = \sqrt{N^2} = \sqrt{\int_{x_1}^{x_2} |\psi_{\rm original}|^2 dx}
# $$
# 
# 3 Multiplique la función de onda original por el inverso de su norma.
# 
# $$
# \psi_{\rm normalizada} = \frac{1}{N} \psi_{\rm original}
# $$
# ```

# In[24]:


#Normalice función de onda


# In[25]:


norm2 = integrate.quad(lambda x: x*x,-3,3)[0]
print("N**2 =",norm2)

norm = np.sqrt(norm2)
print("N =",norm)

integrate.quad(lambda x: (x/norm)*(x/norm),-3,3)[0]


# **OPCIONAL**
# 
# ````{admonition} Aprendizaje de código
# :class: important
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
# \begin{eqnarray*}
# N^2 &=& \int \bigg( \psi^* \psi \bigg)\,d\textbf{r}  = \int\limits_{0}^\pi \int\limits_{0}^{2\pi} \int\limits_{0}^\infty \bigg(\psi^* \psi\, r^2 \sin\theta \bigg) dr d\phi d\theta \\
# &=& \int\limits_{0}^{\pi} \sin \theta d\theta \int\limits_{0}^{2\pi} d\phi \int\limits_{0}^\infty e^{-2r^2} r^2 dr \\
# &=& \left(2\pi\right)\left(2\right)\left(\frac{1}{8}\sqrt{\frac{\pi}{2}}\right) = \left(\frac{\pi}{2}\right)^{3/2}
# \end{eqnarray*}
# 
# El siguiente código evalúa la integral y la guarda en la variable "norm2".
# ~~~python
# norm2 = integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
# ~~~
# ````
# 
# ```{margin}
# Nuestro objetivo en este punto es que conozca lo mucho que se puede hacer en Python. Exploraremos estas ideas más adelante.
# ```

# In[26]:


#Determine el cuadrado de la norma de la función de onda


# In[27]:


norm2 = integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]

print("N**2 =",norm2)


# **Compruebe que numéricamente esto es $\left(\frac{\pi}{2}\right)^{3/2}$.**

# In[28]:


#Evalúe (pi/2)^(3/2)


# In[29]:


print((np.pi/2)**(3/2))


# Evalúe la norma

# In[30]:


# Constante de normalización


# In[31]:


norm = np.sqrt(norm2)
print("N =",norm)


# **OPCIONAL**
# 
# Proponga una $\psi$ normalizada, y compruebe que la integral de su cuadrado da 1. Recuerde la siguiente relación
# 
# $$
# \psi_{\rm normalizada} = \frac{1}{\sqrt{N}} \psi_{\rm original}
# $$
# 
# y que para este ejercicio
# 
# $$
# \psi_{\rm original} = e^{-r}
# $$

# In[32]:


#Comprobación


# In[33]:


norm = integrate.tplquad(lambda theta,phi,r: 1/norm**2*r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
print(norm)


# **OPCIONAL**
# 
# **Repita el proceso de normalización para $\psi=re^{-r^2}$.** La integral sin normalizar será
# 
# $$
# \int \psi^* \psi \,d\textbf{r} = \frac{3}{4}\left(\frac{\pi}{2}\right)^{3/2}
# $$

# In[34]:


#Normalice


# In[35]:


norm2=integrate.tplquad(lambda theta,phi,r: r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
print("N**2 =",norm2)

norm = np.sqrt(norm2)
print("N =",norm)

norm=integrate.tplquad(lambda theta,phi,r: 1/norm**2*r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)[0]
print(norm)

