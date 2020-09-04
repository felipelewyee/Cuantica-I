# Partícula Libre y Partícula en la caja

A continuación estudiaremos sistemas sencillos que nos permiten entender como surge la cuantización. Además, nos permitirán familiarizarnos con los pasos para resolver los problemas de química cuántica. Podemos resumir estos como:
1. Identificar las interacciones y restricciones del sistema.
2. Escrbir el Hamiltoniano ($\mathcal{H}$) y la ecuación de Schrodinger ($\mathcal{H}\psi = \varepsilon \psi$).
3. Encontrar la función de onda ($\psi$).
4. Estudiar las condiciones de cuantización.

## Partícula Libre

La partícula libre es el sistema más sencillo posible, consiste en una partícula que se mueve en un potencial **uniforme**, **conservativo** y **sin restricciones**. 

El primer paso consiste en escribir el Hamiltoniano, el cuál contiene las interacciones del sistema. En 1D el Hamiltoniano es
\begin{equation}
\mathcal{H} = - \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V(x)
\end{equation}
donde el primer término es la energía cinética y el segundo el potencial.

La ecuación de Schrodinger a resolver se obtiene de susituir el Hamiltoniano en $\mathcal{H}\psi=\varepsilon\psi$:
\begin{equation}
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} +V(x)\right)\psi(x)=\varepsilon\psi(x)
\end{eqaution}

Esta se puede reescribir como una ecuación diferencial homogénea de segundo orden
\begin{equation}
-\frac{d^2 \psi(x)}{dx^2} = \frac{2m(\varepsilon-V)}{\hbar^2}\psi(x) = \frac{2mE}{\hbar^2}\psi(x)
\end{equation}

La solución a la ecuación anterior es una función cuya segunda derivada es la misma función multiplicada por alguna constante. Esta característica la cumplen las exponenciales, por tanto se propone que la solución tenga la forma:
\begin{equation}
\psi = A'e^{ikx} + B'e^{-ikx}
\end{equation}

Al sustituir en la ecuación diferencial
\begin{equation}
-\frac{d^2 \psi(x)}{dx^2} = k^2(A'e^{ikx} + B'e^{-ikx}) = k^2\psi
\end{equation}

Podemos encontrar la forma de $k$ por simple comparación:
\begin{equation}
k^2 = \frac{2mE}{\hbar^2}
\end{equation}


Por la identidad de Euler, la función de onda también puede escribirse como:
\begin{equation}
\psi = Acos(kx) + Bsin(kx)
\end{equation}

Los pasos anteriores nos han permitido encontrar una forma para la función de onda, sin embargo, no ha surgido cuantización.

## Partícula en una caja

La versión 1D de este sistema consiste en una partícula que se mueve en el espacio con un potencial definido en tres regiones, tal que

$$
V(x) = \left\{
  \begin{array}{lll}
  \infty      & \mathrm{si\ } x < 0 & I\\
  V & \mathrm{si\ } 0 \le x \le L & II \\
  \infty     & \mathrm{si\ } x > L & III
  \end{array}
  \right.
$$

Esto significa que la partícula está confinada a un intervalo en $x \epsilon [0,L]$. La función de onda se puede dividir por regiones. Es imposible que la partícula se encuentre en la región I y en la región III, ya que el potencial es infinito, por lo tanto:

$$
\psi(x) = 0 \left\{
  \begin{array}{lll}
  \mathrm{si\ } x < 0 & I\\
  \mathrm{si\ } x > L & III
  \end{array}
  \right.
$$

Para encontrar la función de onda en la región II hay que escribir la de Schrodinger
\begin{equation}
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} +V(x)\right)\psi(x)=\varepsilon\psi(x)
\end{equation}


cuya solución, como se vio anteriormente, es
\begin{equation}
  \begin{array}{lll}
  \psi(x) = Acos(kx) + Bsin(kx)      & \mathrm{si\ } 0 \leq x \leq L; & k^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
\end{equation}


La función de onda debe ser continua, esto significa que la región I y la región II deben unirse en el mismo punto, es decir, $\psi_{I}(0) = \psi_{II}(0) = 0$. Esto implica que $A=0$, ya que
\begin{equation}
\psi_{II}(0) = 0 = Acos(0) + Bsin(0) = A
\end{equation}


Por la continuidad con la región III también se cumple $\psi_{II}(L) = \psi_{III}(L) = 0$, es decir
\begin{equation}
\psi_{II}(L) = 0 = Bsin(kL)
\end{equation}


Ya obtuvimos que $A$ vale cero, sin embargo, $B$ no puede ser cero porque $\psi_{II}$ se anularía. La única forma de que se cumpla la ecuación anterior es que $kL$ sea un múltiplo de $\pi$, es decir $kL = n \pi$, o lo que es lo mismo
\begin{equation}
  \begin{array}{ll}
  k=\frac{n\pi}{L};      & n=1,2,3,...\\
  \end{array}
\end{equation}


Hemos obtenido que
\begin{equation}
\psi = B sin \left( \frac{n \pi x}{L}\right)
\end{equation}


Para encontrar el valor de B hay que normalizar la función de onda, resultando que:
\begin{equation}
\psi_n(x)=\left(\frac{2}{L}\right)^{1/2}sin\left(\frac{n\pi x}{L}\right)
\end{equation}


Al sustituir la función de onda se obtiene la expresión de la energía, que es
\begin{equation}
  \begin{array}{ll}
  E = \frac{\hbar^2 k^2}{2m} = \frac{h^2 n^2}{8mL^2};      & n=1,2,3,...\\
  \end{array}
\end{equation}

**Grafique la función de onda y el cuadrado de la función de onda para n=1 para L=4.0 A**

Sugerencias
1. Importe numpy y matplotlib.pyplot
2. Declare la variable L y asígnele un valor, por ejemplo $L=1$
3. Cree el dominio de x de 0 a L con numpy.linspace, utilice una cantidad de puntos, por ejemplo 50.
4. Evalúe la función de onda en el dominio
5. Calcule el cuadrado de la función de onda en el dominio
6. Grafique la función de onda y su cuadrado usando matplotlib y pyplot.

# Gráfica de psi_1 y su cuadrado

import numpy as np
from matplotlib import pyplot as plt

L=4.0
n=1.0

x=np.linspace(0,L,100)

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi

plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

**Grafique la función de onda y el cuadrado de la función de onda para n=1,2,3,4 para L=4.0 A**
Sugerencias:
1. Importe numpy y matplotlib.pyplot
2. Declare la variable L y asígnele un valor, por ejemplo $L=1$
3. Cree el dominio de x de 0 a L con numpy.linspace, utilice una cantidad de puntos, por ejemplo 50.
4. Evalúe las 4 funciones de onda en el dominio
5. Calcule el cuadrado de las 4 funciones de onda en el dominio
6. Grafique las funciones y su cuadrado usando matplotlib y pyplot.

# Gráfica de psi_1, psi_2, psi_3, psi_4 y su cuadrado

import numpy as np
from matplotlib import pyplot as plt

L=4.0
x=np.linspace(0,L,100)

n=1.0

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=2.0

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=3.0

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=4.0

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

Como estamos haciendo una secuencia de gráficas donde aumentamos n de uno en uno, podemos hacerlo con un ciclo for. **Repita la gráfica de la función de onda con $n=1,2,3,4$ utilizando un ciclo for**.

# Gráfica de psi_1, psi_2, psi_3, psi_4 y su cuadrado con for

import numpy as np
from matplotlib import pyplot as plt

L=4.0
x=np.linspace(0,L,100)

for n in range(1,5):
    psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
    psi2=psi*psi
    plt.plot(x,psi,label="psi")
    plt.plot(x,psi2,label="$psi^2$")
    plt.legend()
    plt.axhline(y=0, color='k')
    plt.show()

**Haga la gráfica de E en función de n para los primeros 10 niveles energéticos de un electrón en una caja.**

Recuerde que

\begin{equation}
  E = \frac{h^2 n^2}{8mL^2}
\end{equation}

Con $n=1,2,3,...$

Considere h=1 y m=1. Esto se llama unidades atómicas. Tome L=1 A

from matplotlib import pyplot as plt

for n in range(1,11):
    plt.hlines(n**2.0/(8.0*(1.0**2.0)),0,1)

plt.xlim(0,4)
plt.show()

**Muestre que $\psi_1$ y $\psi_3$ son ortonormales (Tome $L=2.0$).**

Ayuda. Haga la integral

import numpy as np
from scipy import integrate

L=2.0

psi_1psi_1 = integrate.quad(lambda x: np.sqrt(2.0/L)*np.sin(np.pi*1.0*x/L)*np.sqrt(2.0/L)*np.sin(np.pi*1.0*x/L),0,L)
psi_3psi_3 = integrate.quad(lambda x: np.sqrt(2.0/L)*np.sin(np.pi*3.0*x/L)*np.sqrt(2.0/L)*np.sin(np.pi*3.0*x/L),0,L)

psi_1psi_3 = integrate.quad(lambda x: np.sqrt(2.0/L)*np.sin(np.pi*1.0*x/L)*np.sqrt(2.0/L)*np.sin(np.pi*3.0*x/L),0,L)

print("Integrales")
print("psi_1psi_1",psi_1psi_1)
print("psi_3psi_3",psi_3psi_3)
print("psi_1psi_3",psi_1psi_3)

También existe la partícula en una caja para 2-Dimensiones. Se confina la partícula en $x\varepsilon[0,L_x]$ y $y\varepsilon[0,L_y]$.

La ecuación de Schrodinger a resolver es
\begin{equation}
-\frac{\hbar^2}{2m} \left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)\psi(x,y)=E\psi(x,y)
\end{equation}

Que se resuelve por el método de separación de variables y se obtiene:
\begin{equation}
  \begin{array}{ll}
  k_x=\frac{n_x\pi}{L_x};      & n_x=1,2,3,...\\
  k_y=\frac{n_y\pi}{L_y};      & n_y=1,2,3,...\\
  \end{array}
\end{equation}

\begin{equation}
\psi_{n_x,n_y}(x)=\left(\frac{2}{L_x}\right)^{1/2} \left(\frac{2}{L_y}\right)^{1/2}sin\left(\frac{n_x\pi x}{L_x}\right)sin\left(\frac{n_y\pi y}{L_y}\right)
\end{equation}

\begin{equation}
E = \frac{h^2 }{8m} \left(\frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} \right)
\end{equation}

**Obtenga la gráfica de $\psi_{1,1}$, es decir $n_x=1$ y $n_y=1$, y de $|\psi_{1,1}|^2$ con $L_x = L_y = 4.0$**

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

Lx=4.0
Ly=4.0

nx=1.0
ny=1.0

x = np.linspace(0, Lx, 30)
y = np.linspace(0, Ly, 30)
X, Y = np.meshgrid(x, y)

psi = np.sqrt(2.0/Lx)*np.sqrt(2.0/Ly)*np.sin(nx*np.pi*X/Lx)*np.sin(ny*np.pi*Y/Ly)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, psi, rstride=1, cstride=1,
                cmap='YlGnBu', edgecolor='none')
ax.set_title("$\Psi$")
plt.show()

psi = np.sqrt(2.0/Lx)*np.sqrt(2.0/Ly)*np.sin(nx*np.pi*X/Lx)*np.sin(ny*np.pi*Y/Ly)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, psi**2.0, rstride=1, cstride=1,
                cmap='YlGnBu', edgecolor='none')
ax.set_title("$|\Psi|^2$")
plt.show()

**Obtenga las gráfica de $\psi_{3,3}$ y $|\psi_{3,3}|^2$ con $L_x = L_y = 4.0$**

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

Lx=4.0
Ly=4.0

nx=3.0
ny=3.0

x = np.linspace(0, Lx, 70)
y = np.linspace(0, Ly, 70)
X, Y = np.meshgrid(x, y)

psi = np.sqrt(2.0/Lx)*np.sqrt(2.0/Ly)*np.sin(nx*np.pi*X/Lx)*np.sin(ny*np.pi*Y/Ly)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, psi, rstride=1, cstride=1,
                cmap='YlGnBu', edgecolor='none')
ax.set_title("$\Psi$")
plt.show()

psi = np.sqrt(2.0/Lx)*np.sqrt(2.0/Ly)*np.sin(nx*np.pi*X/Lx)*np.sin(ny*np.pi*Y/Ly)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, psi**2.0, rstride=1, cstride=1,
                cmap='YlGnBu', edgecolor='none')
ax.set_title("$|\Psi|^2$")
plt.show()

## Referencias

(1) Pilar, F. L. Elementary Quantum Chemistry, Dover ed.; Dover Publications: Mineola, N.Y, 2001.

(2) Breneman, G. L. The Two-Dimensional Particle in a Box. Journal of Chemical Education 1990, 67 (10), 866.