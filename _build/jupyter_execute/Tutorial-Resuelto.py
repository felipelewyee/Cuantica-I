# Tutorial.

Hola!, si estás viendo esto es que instalaste **Anaconda** y **Python** correctamente, **felicidades!**. Este es un tutorial de como usar Python3 y de algunas de las funciones que nos servirán para el curso de **química cuántica I**. **Para usarlo solo tienes que escribir en los cuadros en blanco y presionar el botón Run o "Shift" + "Enter".**

Al igual que en álgebra, en python se pueden usar variables para **almacenar valores numéricos**, por ejemplo
~~~python
x=5.0
~~~
**En el siguiente recuadro prueba definir una variable y asignarle un valor, puedes usar el nombre que tu quieras, recuerda dar "shift"+"enter" o Run.**

x=5.0

Puedes **ver (o imprimir) el valor de cualquier variable escribiendo**
~~~python
print(variable)
~~~

En nuestro ejemplo
~~~python
print(x)
~~~
Imprimiría 5.0. **Prueba imprimir el valor de la variable que hiciste arriba**

print(x)

## Operacions Básicas

También aplican las operaciones conocidas del álgebra, por ejemplo


| Operación     |    Álgebra    |    Python    |
|---------------|---------------|--------------|
|Suma|     $z=x+y$   |       $z=x+y$|
|Resta|     $z=x-y$   |       $z=x-y$|
|Multiplicación|     $z=(x)(y)$   |       $x=x*y$|
|División|     $z=\frac{x}{y}$   |       $z=x/y$|
|Potencia|     $z=x^{y}$   |       $z=x**y$|
|---------|-----------------------|------------------|

En nuestro ejemplo, podemos escribir una nueva variable, hacer las operaciones, guardar el resultado e imprimirlo
~~~python
#Definimos la variable y
y=2.0
#Hacemos algunas operaciones con nuestras variables x, y; y guardamos el resultado
z_suma=x+y
z_resta=x-y
z_mult=x*y
z_div=x/y
z_pot=x**y
#Imprimimos los resultados
print("x",x,"y",y)
print("suma",z_suma)
print("resta",z_resta)
print("multiplicación",z_mult)
print("división",z_div)
print("potencia",z_pot)
~~~

Los resultados se almacenan en las variables z_suma, z_resta, z_mult y z_div. Nota que una variable puede tener un nombre más largo que una letra; y que también podemos imprimir texto.

**Define una variable extra, calcula las 5 operaciones y obtén el resultado.**


y=2.0
z_suma=x+y
z_resta=x-y
z_mult=x*y
z_div=x/y
z_pot=x**y
print("x=",x,"y=",y)
print("suma",z_suma)
print("resta",z_resta)
print("multiplicación",z_mult)
print("división",z_div)
print("potencia",z_pot)

## Librerías

La funcionalidad de python se puede extender utilizando **librerías**. Existen 2 importantes para el curso, **numpy** y **matplotlib**.
**Para activar una librería utilizamos la palabra import**, por ejemplo
~~~python
import numpy as np
~~~

También podemos tomar partes específicas de un módulo con la palabra **from**. Por ejemplo, importaremos la función pyplot de la librería matplotlib, la cual nos sirve para graficar, y la renombraremos como plt.
~~~python
from matplotlib import pyplot as plt
~~~

**Copia y pega las instrucciones anteriores en el recuadro para obtener las funciones de las librerias.**

import numpy as np
from matplotlib import pyplot as plt

Utilizaremos numpy para obtener la raíz cuadrada de nuestra variable, en el ejemplo nuestra variable x
~~~python
x=5.0
raiz = np.sqrt(x)
print(raiz)
~~~

**Define una variable y obtén e imprime su raiz cuadrada.**

x=5.0
raiz=np.sqrt(x)
print(raiz)

## Vectores y Matrices

Podemos hacer vectores poniéndolos entre corchetes cuadrados, por ejemplo
~~~python
vector = [5, 9, 4]
~~~
Podemos imprimir todo el vector con
~~~python
print(vector)
~~~
O podemos imprimir un elemento específico. En python los elementos se cuentan desde "cero", para imprimir el número 9 haríamos:
~~~python
print(vector[1])
~~~

**Crea un vector de dimensión 3 (osea con 3 números), imprime el vector completo e imprime también el primer elemento del vector.**

vector = [5,9,4]
print(vector)
print(vector[0])

También podemos crear matrices, por ejemplo, vamos a crear una matriz de 3x3 que contenga ceros en todos su elementos
~~~python
matriz = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
~~~
**Prueba haciendo una matriz con valores distintos de cero, e imprime todos los elementos de la matriz.**

matriz = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0]]
print(matriz)

Podemos usar numpy para hacer nuestra matriz de ceros de forma más fácil
~~~python
matriz = np.zeros((3,3))
print(matriz)
~~~
**Prueba usar numpy para definir una matriz de ceros de 5x5, e imprimela.**

matriz = np.zeros((5,5))
print(matriz)

A continuación define una matriz de 3x3 y asígnale los valores del 1 al 9. Por ejemplo.
~~~python
matriz1 = np.zeros((3,3))
matriz1[0][0]=1.0
matriz1[0][1]=2.0
matriz1[0][2]=3.0
matriz1[1][0]=4.0
matriz1[1][2]=5.0
...
~~~
**Escribe y completa el código anterior en el siguiente recuadro, e imprime la matriz.**

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
print(matriz1)

Hacer esto puede tomar mucho tiempo con matrices grandes. Existe una instrucción llamada "for" que nos permite hacer cosas repetitivas. Por ejemplo:
~~~python
matriz2 = np.zeros((3,3))
val=0.0
for i in range(3):
    for j in range(3):
        val=val+1.0
        matriz2[i][j]=val
~~~
**Prueba la instrucción anterior e imprime las dos matrices (matriz1 y matriz2) en el siguiente recuadro.**

matriz2 = np.zeros((3,3))
val=0.0
for i in range(3):
    for j in range(3):
        val=val+1.0
        matriz2[i][j]=val
print(matriz1)
print(matriz2)

Puedes multiplicar dos matrices con
~~~python
matriz3=np.matmul(matriz1,matriz2)
~~~
**Multiplica tus dos matrices anteriores e imprime el resultado.**

matriz3=np.matmul(matriz1,matriz2)
print(matriz3)

Y puedes encontrar sus eigenvalores y eigenvectores con la instrucción
~~~python
val,vec=np.linalg.eig(matriz)
~~~
**Encuentra los eigenvectores y los eigenvalores de tu matriz anterior, e imprimelos.**

val,vec=np.linalg.eig(matriz1)
print("Eigenvalores")
print(val)
print("Eigenvectores")
print(vec)

## Gráficas

Vamos a graficar la función y=sin(x) de -3 a 3.
Primero crearemos el dominio de la función (los valores de x), en nuestro ejemplo le daremos 50 puntos. Utilizaremos linspace(a,b,n) que crea un conjunto de n números distribuidos desde a hasta b.
~~~python
x=np.linspace(-3,3,50)
~~~

Luego obtendremos el valor de y usando numpy
~~~python
y=np.sin(x)
~~~
y graficamos con la instrucción
~~~python
plt.scatter(x,y)
~~~
Para abrir la gráfica
~~~python
plt.show()
~~~

x=np.linspace(-3,3,50)
y=np.sin(x)
plt.scatter(x,y)
plt.show()

**En el siguiente recuadro haga la gráfica de las funciones $y=e^{-|x|}$ y de $y=e^{-|x|^2}$ con 100 puntos desde -3 hasta 3.**

x=np.linspace(-3,3,200)
y1=np.exp(-np.abs(x))
y2=np.exp(-np.abs(x)**2.0)
plt.scatter(x,y1)
plt.scatter(x,y2)
plt.show()

## Integrales

También podemos hacer integrales con python. Por ejemplo, vamos a integrar y=x^2 en el dominio $-3 \leq x \leq3$. Para ello importaremos scipy.integrate.quad con
~~~python
import scipy.integrate as integrate
~~~
y luego realizaremos la ingeral
~~~python
scipy.integrate.quad(lambda x: x^2,-3,3)
~~~
Aquí "lambda" indica las variables de la ecuación, seguido de la ecuación y los límites de la integral. **Pruebe a realizar la integral.**

Si tenemos la funión de onda $\psi = x$, definida con $x \epsilon [-3,3]$, **proponga una función de onda normalizada y evalúe la integral con scipy.integrate.quad para comprobar que la norma es 1.**

import scipy.integrate as integrate
integrate.quad(lambda x: x**2,-3,3)



norm2 = integrate.quad(lambda x: x*x,-3,3)
print(norm2)
#Constante de normalizacion
c=np.sqrt(1.0/18.0)
integrate.quad(lambda x: (c*x)*(c*x),-3,3)

Sea la función de onda
\begin{equation}
\psi = e^{-(x^2+y^2+z^2)}=e^{-r^2}
\end{equation}
Su integral de normalización será
\begin{equation}
\int\limits_{-\infty}^\infty \int\limits_{-\infty}^\infty \int\limits_{-\infty}^\infty \psi^* \psi dx dy dz = \int \psi^* \psi d\textbf{r} = \int\limits_{0}^\infty \int\limits_{0}^{2\pi} \int\limits_{0}^\pi\psi^* \psi r^2 sin\theta dr d\phi d\theta = \int\limits_{0}^{\pi} sin \theta d\theta \int\limits_{0}^{2\pi} d\phi \int\limits_{0}^\infty e^{-2r^2} r^2 dr = \left(2\pi\right)\left(2\right)\left(\frac{1}{8}\sqrt{\frac{\pi}{2}}\right) = \left(\frac{\pi}{2}\right)^{3/2}
\end{equation}

Esta es una triple integral que se resuelve como sigue:
~~~python
gaussian_norm=integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
~~~
**En el cuadro de abajo resuelve la integral copiando la línea mostrada e imprime el resultado**

gaussian_norm=integrate.tplquad(lambda theta,phi,r: r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(gaussian_norm)

**Comprueba que numéricamente esto es $\left(\frac{\pi}{2}\right)^{3/2}$ y propón una $\psi'$ normalizada, compruebe que su norma es uno.**

print((np.pi/2)**(3/2))
gaussian_norm=integrate.tplquad(lambda theta,phi,r: (2/np.pi)**(3/2)*r**2.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(gaussian_norm)

**Repita el proceso de normalización para $\psi=re^{-r^2}$. La integral sin normalizar será
\begin{equation}
\int \psi^* \psi d\textbf{r} = \frac{3}{4}\left(\frac{\pi}{2}\right)^{3/2}
\end{equation}**

gaussian_norm=integrate.tplquad(lambda theta,phi,r: r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(gaussian_norm)
gaussian_norm=integrate.tplquad(lambda theta,phi,r: (3.0/4.0*(np.pi/2.0)**(3/2))**(-1.0)*r**4.0*np.sin(theta)*np.exp(-2.0*r**2.0), 0, np.inf, lambda r: 0, lambda theta: 2*np.pi,lambda r, theta: 0, lambda r, theta: np.pi)
print(gaussian_norm)

## Autores

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.