# Efecto Túnel

## Barrera de potencial $V$ de longitud infinita

Considere una partícula moviéndose hacia una barrera de potencial V en $x=0$. Antes de $x=0$ el potencial vale cero, y después vale $V$, es decir:

$$
V(x) = \left\{
  \begin{array}{lll}
  0      & \mathrm{si\ } x < 0 & I\\
  V & \mathrm{si\ } 0 \le x < \infty & II \\
  \end{array}
  \right.
$$

En este sistema consideraremos el caso en el que la partícula tiene menor energía, $E$, que el potencial, $V$, es decir, $E < V$.

```{admonition} Para pensar
:class: tip
De manera clásica, la partícula no podría pasar del lado izquierdo (región I) al lado derecho (región II) de la caja, porque no tiene suficiente energía. Por esta razón, no podríamos encontrar a la partícula en la región II. ¿Qué pasará cuánticamente?
```

En la región I, la ecuación de Schrödinger es

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} \right)\psi_I(x) = E \psi_I(x)
$$

cuya solución es

$$
\begin{array}{lll}
  \psi_I(x) = Acos(k_1x) + Bsin(k_1x)      & \mathrm{si\ } 0 \le x; & k_1^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
$$

En la región II, la ecuación de Schrödinger es

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V \right)\psi_{II}(x) = E \psi_{II}(x)
$$

cuya solución, si $E < V$ es

$$
\begin{array}{lll}
  \psi_{II}(x) = C e^{-k_2x} + De^{k_2x}      & \mathrm{si\ } 0 \le x; & k_2^2 = \frac{2m(V-E)}{\hbar^2}\\
  \end{array}
$$

Los coeficientes pueden obtenerse a partir de la condicióón de continuidad de la función de onda en $x=0$

```{admonition} Inserto matemático: Condiciones de continuidad
:class: dropdown

Para que la función de onda converga al infinito, se requiere que $D$ sea cero, ya que $e^{k_2x}$ es una función siempre creciente, es decir

$$
D = 0
$$

Se requiere que la función de onda de la zona I y de la zona II sean iguales en $x=0$, es decir

$$
\psi_I(0) = \psi_{II}(0)
$$

$$
A = C
$$

También se requiere que las derivadas de la función de onda de la zona I y de la zona II sean iguales en $x=0$, es decir

$$
\frac{d\psi_I}{dx}(0) = \frac{d\psi_{II}}{dx}(0)
k_1 B = - k_2 C
$$
```
Se obtiene

$$
B = -\frac{k_2}{k_1} A
$$

$$
C = A
$$

$$
D = 0
$$

**Importe numpy y pyplot de matplotlib.**

# Importe librerías

import numpy as np
from matplotlib import pyplot as plt

**De valores a las constantes del sistema**. Considere $m=1$, $\hbar=1$. Asigne algún valor a la energía y al potencial, respetando que $V > E$, observe que en este caso la energía no está cuantizada, por lo que puede tomar cualquier valor. A manera de ejemplo, considere $E=1$ y $V=10$. 

# Valores de m,hbar,E,V

m = 1
hbar = 1
E = 1
V = 10

Defina $k_1$ y $k_2$, recuerde que

$$
k_1 = \frac{\sqrt{2mE}}{\hbar}
$$

$$
k_2 = \frac{\sqrt{2m(V-E)}}{\hbar}
$$

# k1 y k2

k1 = np.sqrt(2*m*E)/hbar
k2 = np.sqrt(2*m*(V-E))/hbar

Defina las constantes

$$
B = -\frac{k_2}{k_1} A
$$

$$
C = A
$$

Por conveniencia, defina 

$$
A=1
$$

# A, B, C

A = 1
B = -k2/k1*A
C = A

Defina el dominio de $x$ para la región I y para la región II, recuerde que ambos se separan en $x=0$.

# x1 y x2

x1 = np.linspace(-10,0,100)
x2 = np.linspace(0,10,100)

Genere la función de onda para la región I y para la región II. Recuerde

$$
\psi_I = A cos(k_1 x) + B sin(k_1 x)
$$

$$
\psi_{II} = C e^{-k_2 x}
$$

# psi_I y psi_II

psi_I = A*np.cos(k1*x1) + B *np.sin(k1*x1)
psi_II = C*np.exp(-k2*x2)

Grafique $\psi_I$ y $\psi_{II}$

# Grafica

plt.plot(x1,psi_I)
plt.plot(x2,psi_II)

Note que la partícula puede ser encontrada dentro de la región clásicamente prohibida, esto se conoce como penetración. 

```{admonition} Concepto: Longitud de decaimiento
:class: note

La longitud de decaimiento, $\frac{1}{k_2}$ es la distancia dentro de la barrera a la cual la función de onda a deacído a  $\frac{1}{e}$.
```

**Calcule la longitud de decaimiento de este sistema.**

## Barrera de potencial $V$ de longitud finita

Considere una partícula moviéndose hacia una barrera de potencial V en $0 \leq x \leq  L$, es decir:

$$
V(x) = \left\{
  \begin{array}{lll}
  0      & \mathrm{si\ } x < 0 & I\\
  V & \mathrm{si\ } 0 \leq x \leq L & II \\
  0 & \mathrm{si\ } L \le x & II \\
  \end{array}
  \right.
$$

En este sistema consideraremos el caso en el que la partícula tiene menor energía, $E$, que el potencial, $V$, es decir, $E < V$.

En la región I, la ecuación de Schrödinger es

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} \right)\psi_I(x) = E \psi_I(x)
$$

cuya solución es

$$
\begin{array}{lll}
  \psi_I(x) = Acos(k_1x) + Bsin(k_1x)      & \mathrm{si\ } 0 \le x; & k_1^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
$$

En la región II, la ecuación de Schrödinger es

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V \right)\psi_{II}(x) = E \psi_{II}(x)
$$

cuya solución, si $E < V$ es

$$
\begin{array}{lll}
  \psi_{II}(x) = C e^{-k_2x} + De^{k_2x}      & \mathrm{si\ } 0 \le x; & k_2^2 = \frac{2m(V-E)}{\hbar^2}\\
  \end{array}
$$

En la región III, la ecuación de Schrödinger es

$$
\left(- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} \right)\psi_{III}(x) = E \psi_{III}(x)
$$

cuya solución es

$$
\begin{array}{lll}
  \psi_I(x) = Fcos(k_1x) + Gsin(k_1x)      & \mathrm{si\ } 0 \le x; & k_1^2 = \frac{2mE}{\hbar^2}\\
  \end{array}
$$

Los coeficientes pueden obtenerse a partir de la condicióón de continuidad de la función de onda en $x=0$

```{admonition} Inserto matemático: Condiciones de continuidad
:class: dropdown

Se requiere que la función de onda de la zona I y de la zona II sean iguales en $x=0$, es decir

$$
\psi_I(0) = \psi_{II}(0)
$$

$$
A = C + D
$$

También se requiere que las derivadas de la función de onda de la zona I y de la zona II sean iguales en $x=0$, es decir

$$
\frac{d\psi_I}{dx}(0) = \frac{d\psi_{II}}{dx}(0)
$$

$$
k_1 B = - k_2 (C - D)
$$

Se requiere que la función de onda de la zona II y de la zona III sean iguales en $x=L$, es decir

$$
\psi_{II}(L) = \psi_{III}(L)
$$

$$
Ce^{k_2L} + De^{-k_2L} = F cos(k_1L) + G sin(k_1L)
$$


También se requiere que las derivadas de la función de onda de la zona II y de la zona III sean iguales en $x=L$, es decir

$$
\frac{d\psi_{II}}{dx}(L) = \frac{d\psi_{III}}{dx}(L)
$$

$$
-k_2 C e^{-k_2L} + k_2 D e^{k_2L}  = - k_1 F sin(k_1L) + k_1 G cos(k_1L) 
$$

Al dividir la última ecuación entre $k_2$, y sumarla con la penúltima se obtiene

$$
D = \frac{F(cos(k_1L) - \frac{k_1}{k_2}sin(k_1L) ) + G(sin(k_1L) + \frac{k_1}{k_2}cos(k_1L) )}{2 e^{k_2L}}
$$

y al restarla

$$
C = \frac{F(cos(k_1L) + \frac{k_2}{k_1}sin(k_1L) ) + G(sin(k_1L) - \frac{k_2}{k_1}cos(k_1L) )}{2 e^{-k_2L}}
$$
```

Se obtiene

$$
A = C + D
$$

$$
B = - \frac{k_2}{k_1} (C - D)
$$

$$
C = \frac{F(cos(k_1L) + \frac{k_1}{k_2}sin(k_2L) ) + G(sin(k_1L) - \frac{k_1}{k_2}cos(k_1L) )}{2 e^{-k_2L}}
$$

$$
D = \frac{F(cos(k_1L) - \frac{k_1}{k_2}sin(k_1L) ) + G(sin(k_1L) + \frac{k_1}{k_2}cos(k_1L) )}{2 e^{k_2L}}
$$


**De valores a las constantes del sistema**. Considere $m=1$, $\hbar=1$. Asigne algún valor a la energía y al potencial, respetando que $V > E$, observe que en este caso la energía no está cuantizada, por lo que puede tomar cualquier valor. A manera de ejemplo, considere $E=1$ y $V=10$. Asigne un valor a la longitud de la berrera de potencial, por ejemplo, $L=1$.

# m, hbar, L, E, v

m = 1
hbar = 1
L = 1
E = 1
V = 10

Asigne valores a las constantes $A$, $B$, $C$, $D$, $F$, $G$. Recuerde.

$$
A = C + D
$$

$$
B = - \frac{k_2}{k_1} (C - D)
$$

$$
C = \frac{F(cos(k_1L) + \frac{k_1}{k_2}sin(k_1L) ) + G(sin(k_1L) - \frac{k_1}{k_2}cos(k_1L) )}{2 e^{-k_2L}}
$$

$$
D = \frac{F(cos(k_1L) - \frac{k_1}{k_2}sin(k_1L) ) + G(sin(k_1L) + \frac{k_1}{k_2}cos(k_1L) )}{2 e^{k_2L}}
$$

Por conveniencia, considere $F=1$ y $G=1$.

# A, B, C, D, F, G

F = 1
G = 1

C = (F*(np.cos(k1*L) + k1/k2*np.sin(k1*L)) + G*(np.sin(k1*L) - k1/k2*np.cos(k1*L)))/(2*np.exp(-k2*L))
D = (F*(np.cos(k1*L) - k1/k2*np.sin(k1*L)) + G*(np.sin(k1*L) + k1/k2*np.cos(k1*L)))/(2*np.exp(k2*L))

A = C + D

B = -k1/k2*(C-D)

Defina un dominio para x.

# Dominio de x para región I, II y III

x1 = np.linspace(-10,0,100)
x2 = np.linspace(0,L,100)
x3 = np.linspace(L,10,100)

Defina la función de onda para las tres regiones. Recuerde

$$
\psi_I = A cos(k_1x) + B sin(k_1x)
$$

$$
\psi_{II} = C e^{-k_2 x} + D e^{k_2 x}
$$

$$
\psi_{III} = F cos(k_1x) + G sin(k_1x)
$$

# Función de onda para región I, II y III

psi_I = A*np.cos(k1*x1) + B*np.sin(k1*x1)
psi_II = C*np.exp(-k2*x2) + D*np.exp(k2*x2)
psi_III = F*np.cos(k1*x3) + G*np.sin(k1*x3)

Grafique la función de onda.

# Grafique función de onda

plt.plot(x1,psi_I)
plt.plot(x2,psi_II)
plt.plot(x3,psi_III)