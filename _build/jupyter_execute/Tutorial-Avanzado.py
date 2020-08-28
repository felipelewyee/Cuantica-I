# Tutorial Avanzado.

En esta parte aprenderá a utilizar Python para manejar vectores y matrices, crear gráficas y resolver integrales. Estas funciones nos serán trascendentales para el curso de **Química Cuántica I**. **Recuerde escribir en los cuadros vacíos y presionar el botón Run o las teclas "Shift" + "Enter".**

Utilizaremos las librerías numpy, matplotlib y scipy que aprendió a utilizar en el Tutorial Básico, impórtelas utilizando las líneas
~~~python
import numpy as np
from matplotlib import pyplot as plt
import scipy
~~~

# Importe librerías


## Vectores y Matrices

Python permite manejar fácilmente objetos multidimensionales, como vectores y matrices. Empezaremos creando una lista, la cual es un conjunto de cosas. Las listas se crean poniendo elementos entre corchetes cuadrados, por ejemplo
~~~python
lista = [5, 9, 4]
~~~

Podemos imprimir todos los elementos de la lista con
~~~python
print(lista)
~~~
También podemos imprimir un elemento específico indicando el número del elemento entre paréntesis cuadrados. En Python los elementos se cuentan desde "cero". Para imprimir el número 9 en la lista anterior ejecutaríamos:
~~~python
print(vector[1])
~~~

**Cree la siguiente lista [8,25,32,41] y guárdela en la variable "lista". Imprima la lista completa, así como el primer elemento del vector.**

#Cree una lista


La lista puede convertirse en un arreglo de numpy, esto facilita el manejo de cantidades vectoriales. Para declarar la lista anterior como un arreglo de numpy escribiremos
~~~python
vector = np.array([5, 9, 4])
~~~

Una ventaja es que podemos declarar dos arreglos, por ejemplo
~~~python
v1 = np.array([2, 4, 6])
v2 = np.array([3, 5, 7])
~~~

y sumarlos, de una forma muy parecida a como se haría en física.

~~~python
print(v1+v2)
~~~

**Cree dos vectores como arreglos de numpy, uno que contenga los elementos "1,8,6", y otro que contenga "5,2,7", e imprima su suma**

#Suma de dos vectores


También podemos crear matrices como una "lista de listas", donde cada lista representa un renglón de la matriz. A continuación generemos una matriz de 3x3 que contenga ceros en todos su elementos
~~~python
matriz = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
~~~
**Genere la matriz de ceros anterior, e imprímala.**

#Cree matriz


Podemos usar numpy para hacer nuestra matriz de ceros de forma más fácil
~~~python
matriz = np.zeros((3,3))
print(matriz)
~~~
**Prueba usar numpy para definir una matriz de ceros de 5x5, e imprímala.**

#Cree matriz de ceros con numpy

**A continuación genere una matriz de 3x3, guárdela en la variable "matriz1", y asígnele los valores del 1 al 9, tal que **

$$
matriz1 = \begin{pmatrix}
1 & 2 & 3\\
4 & 5 & 6\\
7 & 8 & 9\\
\end{pmatrix}
$$

Puede ayudarse del siguiente código.
~~~python
matriz1 = np.zeros((3,3))
matriz1[0][0]=1.0
matriz1[0][1]=2.0
matriz1[0][2]=3.0
matriz1[1][0]=4.0
matriz1[1][2]=5.0
...
~~~

#Llene una matriz de 3x3

Hacer esto puede tomar mucho tiempo con matrices grandes. Existe una instrucción llamada **for** que nos permite hacer cosas repetitivas. Por ejemplo:
~~~python
matriz2 = np.zeros((3,3))
val=0.0
for i in range(3):
    for j in range(3):
        val=val+1.0
        matriz2[i][j]=val
~~~
**Prueba la instrucción anterior e imprime las dos matrices (matriz1 y matriz2) en el siguiente recuadro.**

#Matriz de 3x3 con for

Se pueden multiplicar dos matrices con
~~~python
matriz3=np.matmul(matriz1,matriz2)
~~~
**Multiplique las dos matrices anteriores e imprima el resultado.**

#Matriz3 = Matriz1 x Matriz2

Los valores propios ($\lambda$) y vectores propios ($v$) de una matriz ($M$) son muy importantes en química cuántica. Cumplen la propiedad de que al multipicar la matriz por el vector propio resulta el mismo vector multiplciado por una constante, es decir:

$$
M v = \lambda v
$$

se pueden encontrar con
~~~python
val,vec=np.linalg.eig(matriz)
~~~

**Encuentre los valores propios y vectores propios de la matriz 1, e imprimalos.**



## Gráficas

Vamos a graficar la función y=sin(x) de -3 a 3.

Primero crearemos el dominio de la función (los valores de x), en nuestro ejemplo le daremos 50 puntos. Utilizaremos la función linspace(a,b,n) de numpy, esta crea un conjunto de n números distribuidos en el intervalo de a hasta b.
~~~python
x=np.linspace(-3,3,50)
~~~

**Genere el dominio tal que $x\in[-3,3]$, utilice 50 puntos, guárdelo en la variable x e imprímalo.**

# Cree dominio de x

Obtendremos el valor de "y" usando numpy, al escribir la siguiente instrucción Python toma cada valor de x, le aplica la función y lo guarda en la ariable "y".
~~~python
y=np.sin(x)
~~~

**Evalúe la función**

$$
y = sin(x)
$$

# Evalúe función

Finalmente, **realice la gráfica con**

~~~python
plt.scatter(x,y)
plt.show()
~~~

# Genere gráfica

Las instrucciones anteriores son los pasos básicos para generar una gráfica. **En el siguiente recuadro genere la gráfica de las funciones $y_1=e^{-|x|}$ y $y_2=e^{-|x|^2}$ con 100 puntos desde -3 hasta 3.**

# Gráfica

## Integrales

También podemos hacer integrales con Python. Por ejemplo, vamos a integrar $y=x^2$ en el dominio $-3 \leq x \leq3$. Para ello importaremos el subpaquete integrate de la librería scipy.
~~~python
from scipy import integrate
~~~
y luego realizaremos la ingeral con la función quad
~~~python
integrate.quad(lambda x: x**2,-3,3)
~~~
Aquí "lambda" indica las variables de la ecuación, seguido de la ecuación y los límites de la integral. **Pruebe a realizar la integral.**

# Integre y = x**2 de -3 a 3

Considere la funión de onda $\psi = x$, con $x \epsilon [-3,3]$, **proponga una función de onda normalizada y evalúe la integral con integrate.quad para comprobar que la norma es 1.**

Recuerde que para normalizar:
1. Integre el cuadrado de la función de onda en el dominio.
2. Calcule la constante de normalización

$$
N = \frac{1}{\sqrt{norm}} = \frac{1}{\sqrt{\int_{x_1}^{x_2} |\psi_{original}|^2 dx}}
$$

3. Evalúe la integral del cuadrado de la función de onda normalizada para comprobar que es igual a 1. Recuerde que la función de onda normalizada es:

$$
\psi_{normalizada} = N \psi_{original}
$$

#Normalice función de onda

**OPCIONAL**

También es posible realizar integrales triples, el siguiente ejemplo solo es demostrativo, **simplemente copie y pegue para realizar la integral.**

Sea la función de onda
\begin{equation}
\psi = e^{-(x^2+y^2+z^2)}=e^{-r^2}
\end{equation}

La integral de su cuadrado será
\begin{equation}
\int \psi^* \psi d\textbf{r} = \int\limits_{0}^\infty \int\limits_{0}^{2\pi} \int\limits_{0}^\pi\psi^* \psi r^2 sin\theta dr d\phi d\theta = \int\limits_{0}^{\pi} sin \theta d\theta \int\limits_{0}^{2\pi} d\phi \int\limits_{0}^\infty e^{-2r^2} r^2 dr = \left(2\pi\right)\left(2\right)\left(\frac{1}{8}\sqrt{\frac{\pi}{2}}\right) = \left(\frac{\pi}{2}\right)^{3/2}
\end{equation}

El siguiente código evalúa la integral y la guarda en la variable "norm".
~~~python
norm=integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
~~~

#Encuentre norma de la funión de onda

**Compruebe que numéricamente esto es $\left(\frac{\pi}{2}\right)^{3/2}$.**

#Evalúe (pi/2)^(3/2)

**OPCIONAL**

Proponga una $\psi$ normalizada, y compruebe que la integral de su cuadrado da 1. Recuerde la siguiente relación

$$
\psi_{normalizada} = \frac{1}{\sqrt{norm}} \psi_{original}
$$

y que para este ejercicio

$$
\psi_{original} = e^{-r}
$$

#Comprobación

**OPCIONAL**

**Repita el proceso de normalización para $\psi=re^{-r^2}$.** La integral sin normalizar será

$$
\int \psi^* \psi d\textbf{r} = \frac{3}{4}\left(\frac{\pi}{2}\right)^{3/2}
$$

#Normalice

## Autores

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.