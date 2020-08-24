# Partícula Libre

Es la descripción de una partícula moviéndose en un potencial **uniforme**, **conservativo** y **sin restricciones**.

En 1D el Hamiltoniano es
\begin{equation}
\mathcal{H} = - \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V(x)
\end{equation}

La ecuación de Schrodinger a resolver se obtiene de susituir el Hamiltoniano en $\mathcal{H}\psi=\varepsilon\psi$:
\begin{equation}
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} +V(x)\right)\psi(x)=\varepsilon\psi(x)
\end{equation}

Esta se puede reescribir como una ecuación diferencial homogénea de segundo orden
\begin{equation}
-\frac{d^2 \psi(x)}{dx^2} = \frac{2m(\varepsilon-V)}{\hbar^2}\psi(x) = \frac{2mE}{\hbar^2}\psi(x)
\end{equation}

Se necesita una función cuya segunda derivada sea la misma función multiplicada por alguna constante. Esta característica la cumplen las exponenciales, por tanto:
\begin{equation}
\psi = A'e^{ikx} + B'e^{-ikx}
\end{equation}

Al sustituir
\begin{equation}
-\frac{d^2 \psi(x)}{dx^2} = k^2(A'e^{ikx} + B'e^{-ikx}) = k^2\psi
\end{equation}

Por comparación:
**\begin{equation}
k^2 = \frac{2mE}{\hbar^2}
\end{equation}**

Por la identidad de Euler, la función de onda también puede escribirse como:
**\begin{equation}
\psi = Acos(kx) + Bsin(kx)
\end{equation}**

# Partícula en una caja

En este caso el potencial está definido por:
\begin{equation}
V(x) = \left\{
  \begin{array}{lll}
  \infty      & \mathrm{si\ } x < 0 & I\\
  V & \mathrm{si\ } 0 \le x \le L & II \\
  \infty     & \mathrm{si\ } x > L & III
  \end{array}
  \right.
\end{equation}

Esto significa que la partícula está confinada a un intervalo en $x \epsilon [0,L]$. La función de onda se puede dividir por regiones, tal que
\begin{equation} 
  \begin{array}{ll}
  \psi(x) = 0      & \mathrm{si\ } x < 0 \mathrm{\ o\ } x>L\\
  \end{array}
\end{equation}

La región II tiene una función de onda dada por el Hamiltoniano
\begin{equation}
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} +V(x)\right)\psi(x)=\varepsilon\psi(x)
\end{equation}
cuya solución conocida es
\begin{equation} 
  \begin{array}{lll}
  \psi(x) = Acos(kx) + Bsin(kx)      & \mathrm{si\ } 0 \leq x \leq L; & k^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
\end{equation}

Por la continuidad con la región I, se cumple $A=0$ ya que
\begin{equation}
\psi(0) = 0 = Acos(0) + Bsin(0) = A
\end{equation}
Por la continuidad con la región III se tiene
\begin{equation}
\psi(L) = 0 = Bsin(kL)
\end{equation}

$B \neq 0$ porque $\psi$ se anularía, entonces $kL$ debe ser un múltiplo de $n\pi$, es decir
\begin{equation} 
  \begin{array}{ll}
  k=\frac{n\pi}{L};      & n=1,2,3,...\\
  \end{array}
\end{equation}
Al normalizar la función de onda:
\begin{equation}
\psi_n(x)=\left(\frac{2}{L}\right)^{1/2}sin\left(\frac{n\pi x}{L}\right)
\end{equation}

La energía será
\begin{equation} 
  \begin{array}{ll}
  E = \frac{\hbar^2 k^2}{2m} = \frac{h^2 n^2}{8mL^2};      & n=1,2,3,...\\
  \end{array}
\end{equation}


**Grafique la función de onda y el cuadrado de la función de onda para n=1 para L=4.0 A**

Sugerencias
1. Importe numpy y matplotlib.pyplot
2. Defina la variable L
3. Cree el dominio de 0 a L con numpy.linspace
4. Calcule la función de onda en el dominio
5. Calcule el cuadrado de la función de onda en el dominio
6. Grafique la función de onda y su cuadrado usando matplotlib y pyplot.




**Grafique la función de onda y el cuadrado de la función de onda para n=1,2,3,4 para L=4.0 A**
Sugerencias:
1. Importe numpy y matplotlib.pyplot
2. Defina la variable L
3. Cree el dominio de 0 a L con numpy.linspace
4. Calcule las 4 funciones de onda en el dominio
5. Calcule el cuadrado de las 4 funciones de onda en el dominio
6. Grafique las funciones y su cuadrado usando matplotlib y pyplot.



Como estamos haciendo una secuencia de gráficas donde aumentamos n de uno en uno, podemos hacerlo con un ciclo for



**Haga la gráfica de E en función de n para los primeros 10 niveles energéticos de un electrón en una caja.**

Recuerde que
\begin{equation} 
  E = \frac{h^2 n^2}{8mL^2}
\end{equation}
Con $n=1,2,3,...$

Considere h=1 y m=1. Esto se llama unidades atómicas. Tome L=1 A



**Muestre que $\psi_1$ y $\psi_3$ son ortonormales (Tome $L=2.0$).**

Ayuda. Haga la integral



También existe la partícula en una caja para 2-Dimensiones. Se confina la partícula en $x\varepsilon[0,L_x]$ y $y\varepsilon[0,L_y]$.

La ecuación de Schrodinger a resolver es
\begin{equation}
- \frac{\hbar^2}{2m} \left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)\psi(x,y)=E\psi(x,y)
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



# Referencias

(1) Pilar, F. L. Elementary Quantum Chemistry, Dover ed.; Dover Publications: Mineola, N.Y, 2001.

(2) Breneman, G. L. The Two-Dimensional Particle in a Box. Journal of Chemical Education 1990, 67 (10), 866.

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.