# Partícula en la caja

A continuación estudiaremos sistemas sencillos que nos permiten entender como surge la cuantización. Además, nos permitirán familiarizarnos con los pasos para resolver los problemas de química cuántica. Podemos resumir estos como:
1. Identificar las interacciones y restricciones del sistema.
2. Escrbir el Hamiltoniano ($\mathcal{H}$) y la ecuación de Schrodinger ($\mathcal{H}\psi = \varepsilon \psi$).
3. Encontrar la función de onda ($\psi$).
4. Estudiar las condiciones de cuantización.

## Caja 1D

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

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} +V(x)\right)\psi(x)=\varepsilon\psi(x)
$$

cuya solución es

$$
\begin{array}{lll}
  \psi(x) = Acos(kx) + Bsin(kx)      & \mathrm{si\ } 0 \leq x \leq L; & k^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
$$

posteriormente hay que recurrir a las condiciones a la frontera.

```{admonition} Inserto matemático: Condiciones a la frontera
:class: dropdown

La función de onda debe ser continua, esto significa que la región I y la región II deben unirse en el mismo punto, es decir, $\psi_{I}(0) = \psi_{II}(0) = 0$. Esto implica que $A=0$, ya que

$$
\psi_{II}(0) = 0 = Acos(0) + Bsin(0) = A
$$

Por la continuidad con la región III también se cumple $\psi_{II}(L) = \psi_{III}(L) = 0$, es decir

$$
\psi_{II}(L) = 0 = Bsin(kL)
$$

Ya obtuvimos que $A$ vale cero, sin embargo, $B$ no puede ser cero porque $\psi_{II}$ se anularía. La única forma de que se cumpla la ecuación anterior es que $kL$ sea un múltiplo de $\pi$, es decir $kL = n \pi$, o lo que es lo mismo

$$
\begin{array}{ll}
  k=\frac{n\pi}{L};      & n=1,2,3,...\\
  \end{array}
$$

Hemos obtenido que

$$
\psi = B sin \left( \frac{n \pi x}{L}\right)
$$
```

Para encontrar el valor de B hay que normalizar la función de onda, resultando que:

$$
\psi_n(x)=\left(\frac{2}{L}\right)^{1/2}sin\left(\frac{n\pi x}{L}\right)
$$

Al sustituir la función de onda se obtiene la expresión de la energía, que es

$$
  \begin{array}{ll}
  E = \frac{\hbar^2 k^2}{2m} = \frac{h^2 n^2}{8mL^2};      & n=1,2,3,...\\
  \end{array}
$$

**Importe las siguientes librerías**
- numpy
- pyplot de matplotlib

#librerias

import numpy as np
from matplotlib import pyplot as plt

Considere un electrón dentro de una caa de longitud 4 angstroms. Defina las siguientes constantes
```
hbar = 1
m = 1
L = 4
```

# Constantes

hbar = 1
m = 1
L = 4

**Grafique la función de onda ($\psi$) y su cuadrado ($\psi^2$) para n=1 y L=4.0 A**

```{tip}
1. Declare la variable n asígnele su valor.
2. Cree el dominio de x de 0 a L con numpy.linspace, utilice una cantidad de puntos, por ejemplo 100.
3. Evalúe la función de onda en el dominio
4. Calcule el cuadrado de la función de onda en el dominio
5. Grafique la función de onda y su cuadrado usando matplotlib y pyplot.
```

# Inserte código para gráfica

# Gráfica de psi_1 y su cuadrado

n=1

x=np.linspace(0,L,100)

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi

plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

**Grafique la función de onda ($\psi$) y su cuadrado ($\psi^2$) para n=1,2,3,4 para L=4.0 A**
```{tip}
1. Cree el dominio de x de 0 a L con numpy.linspace, utilice una cantidad de puntos, por ejemplo 100.
2. Evalúe las 4 funciones de onda en el dominio
3. Calcule el cuadrado de las 4 funciones de onda en el dominio
4. Grafique las funciones y su cuadrado usando matplotlib y pyplot.
```

# Inserte código para gráfica

# Gráfica de psi_1, psi_2, psi_3, psi_4 y su cuadrado

x=np.linspace(0,L,100)

n=1

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=2

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=3

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

n=4

psi=np.sqrt(2.0/L)*np.sin(n*np.pi*x/L)
psi2=psi*psi
plt.plot(x,psi,label="psi")
plt.plot(x,psi2,label="psi^2")

plt.legend()
plt.axhline(y=0, color='k')
plt.show()

Como estamos haciendo una secuencia de gráficas donde aumentamos n de uno en uno, podemos hacerlo con un ciclo for. **Repita la gráfica de la función de onda con $n=1,2,3,4$ utilizando un ciclo for**.

# Inserte código para 4 gráficas en las que solo cambia el valor de n, use un for

# Gráfica de psi_1, psi_2, psi_3, psi_4 y su cuadrado con for

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

```{tip}
$$
E = \frac{\hbar^2 \pi^2 n^2}{2mL^2}
$$

Con $n=1,2,3,...$
```

Considere h=1 y m=1. Esto se llama unidades atómicas. Tome L=1 A

# Inserte código para gráfica

L = 1

for n in range(1,11):
    plt.hlines(hbar**2*np.pi**2*n**2/(2*m*L**2),0,1)

plt.xlim(0,4)
plt.show()

**Muestre que $\psi_1$ y $\psi_3$ son ortonormales (Tome $L=2.0$).**

```{tip}
Haga las integrales

$$
\int_0^L \psi_1 \psi_1 dx = 1
\int_0^L \psi_3 \psi_3 dx = 1
\int_0^L \psi_1 \psi_3 dx = 0
$$

```

# Integral

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

## Caja 2D

También existe la partícula en una caja para 2-Dimensiones. Se confina la partícula en $x\varepsilon[0,L_x]$ y $y\varepsilon[0,L_y]$.

La ecuación de Schrodinger a resolver es

$$
-\frac{\hbar^2}{2m} \left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right)\psi(x,y)=E\psi(x,y)
$$

Que se resuelve por el método de separación de variables y se obtiene:

$$
\begin{array}{ll}
  k_x=\frac{n_x\pi}{L_x};      & n_x=1,2,3,...\\
  k_y=\frac{n_y\pi}{L_y};      & n_y=1,2,3,...\\
  \end{array}
$$

$$
\psi_{n_x,n_y}(x)=\left(\frac{2}{L_x}\right)^{1/2} \left(\frac{2}{L_y}\right)^{1/2}sin\left(\frac{n_x\pi x}{L_x}\right)sin\left(\frac{n_y\pi y}{L_y}\right)
$$

$$
E = \frac{h^2 }{8m} \left(\frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} \right)
$$

Para hacer gráficas 3D, importe la siguiente librería
```
from mpl_toolkits import mplot3d
```

# librería

from mpl_toolkits import mplot3d

**Obtenga la gráfica de $\psi_{1,1}$, es decir $n_x=1$ y $n_y=1$, y de $|\psi_{1,1}|^2$ con $L_x = L_y = 4.0$**

# Inserte código para gráfica

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

# Inserte código para gráfica

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