# Penetración

Considere una partícula moviéndose hacia una barrera de potencial de longitud infinita de valor V en $x=0$. Antes de $x=0$ el potencial vale cero, y después vale $V$, es decir:

$$
V(x) = \left\{
  \begin{array}{lll}
  0      & \mathrm{si\ } x < 0 & I\\
  V & \mathrm{si\ } 0 \le x < \infty & II \\
  \end{array}
  \right.
$$

<img src="images/tunel-barrera-infinita.png" alt="Figura de tunel de barrera infinita" width="300"/>

En este sistema consideraremos el caso en el que la partícula tiene menor energía, $E$, que el potencial, $V$, es decir, $E < V$.

```{admonition} Para pensar
:class: tip
De manera clásica, la partícula no podría pasar del lado izquierdo (región I) al lado derecho (región II) de la caja, porque no tiene suficiente energía. Por esta razón, no podríamos encontrar a la partícula en la región II. ¿Qué pasará cuánticamente?
```

Para resolver el sistema hay que planear el Hamiltoniano por regiones y resolver una función de onda para cada región.

```{admonition} Inserto matemático: Hamiltoniano por regiones
:class: dropdown

| Región      | Hamiltoniano | Función de onda | Constantes |
|:----------------:|:---------:|:--------:|:--------:|
| I | $- \frac{\hbar^2}{2m} \frac{d^2}{dx^2}$ | $\psi_I(x) = Ae^{ik_1x} + Be^{-ik_1x}$ | $k_1^2 = \frac{2mE}{\hbar^2}$ |
| II| $- \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V$ | $\psi_{II}(x) = C e^{-k_2x} + De^{k_2x}$ | $k_2^2 = \frac{2m(V-E)}{\hbar^2}$ |
```

Se obtienen las funciones de onda

$$
\psi_I(x) = Ae^{ik_1x} + Be^{-ik_1x}
$$

$$
\psi_{II}(x) = C e^{-k_2x} + De^{k_2x}
$$

Los coeficientes pueden obtenerse a partir de la condición de continuidad de la función de onda en $x=0$

```{admonition} Inserto matemático: Condiciones de continuidad
:class: dropdown

| Regiones | Condición | Ecuación |
|:---: |:---: | :---:|
| II | $\psi_{II}(\infty) = 0$ | $D = 0$|
| I y II | $\psi_{I}(0) = \psi_{II}(0)$ | $A + B = C$ |
| I y II | $\frac{\psi_{I}}{dx}(0) = \frac{\psi_{II}}{dx}(0)$ | $ik_1 (A - B) = - k_2 C$|
```

Se obtiene

$$
B = -\left(\frac{k_2 + ik_1}{k_2 - ik_1}\right) A
$$

$$
C = \frac{2ik_1}{ik_1 - k_2} A
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
B = -\left(\frac{k_2 + ik_1}{k_2 - ik_1}\right) A
$$

$$
C = \frac{2ik_1}{ik_1 - k_2} A
$$

Por conveniencia, defina 

$$
A=1
$$

# A, B, C

A = 1
B = -((k2 + 1j*k1)/(k2 - 1j*k1))*A
C = 2*1j*k1/(1j*k1-k2)*A

Defina el dominio de $x$ para la región I y para la región II, recuerde que ambos se separan en $x=0$.

# x1 y x2

x1 = np.linspace(-2,0,100)
x2 = np.linspace(0,2,100)

Genere la función de onda para la región I y para la región II. Recuerde

$$
\psi_I = A e^{ik_1 x} + B e^{-ik_1 x}
$$

$$
\psi_{II} = C e^{-k_2 x}
$$

# psi_I y psi_II

psi_I = A*np.exp(1j*k1*x1) + B*np.exp(-1j*k1*x1)
psi_II = C*np.exp(-k2*x2)

Grafique $|\psi_I|^2$ y $|\psi_{II}|^2$

# Grafica

plt.plot(x1,abs(psi_I)**2)
plt.plot(x2,abs(psi_II)**2)

Note que la partícula puede ser encontrada dentro de la región clásicamente prohibida, esto se conoce como penetración. 

```{admonition} Concepto: Longitud de decaimiento
:class: note

La longitud de decaimiento, $\frac{1}{k_2}$ es la distancia dentro de la barrera a la cual la función de onda a deacído a  $\frac{1}{e}$.
```

**Calcule la longitud de decaimiento de este sistema.**