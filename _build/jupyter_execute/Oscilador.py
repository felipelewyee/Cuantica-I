# Oscilador Armónico Cuántico

El sistema es una partícula que se mueve en un potencial dado por

$$
V = \frac{1}{2} kx^2
$$

donde $k$ es una constante. 

```{note}
En mecánica clásica podemos encontrar sistemas con este potencial como los resortes, mientras que en química cuántica esto sirve para modelar la vibración de los enlaces.
```

En este caso, el Hamiltoniano contiene la energía cinética y el potencial

$$
H = -\frac{\hbar^2}{2m} \frac{d^2}{dx^2} + \frac{1}{2} kx^2
$$

y la ecuación de Schrödinger a resolver es

$$
\left(-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} + \frac{1}{2} kx^2 \right) \psi = E \psi
$$

esta es una ecuación cuyas soluciones a la energía y la función de onda son

$$
E_n = \left(n + \frac{1}{2}\right) h \nu
$$

$$
\psi_n = N_n H_n(\alpha x) e^{-\alpha^2x^2/2}
$$

aquí $n$ es un número cuántico tal que $n = 0,1,2,3,...$, y los términos $N_n$ y $\alpha$ están dados por

$$
N_n = \left(\frac{\alpha}{2^n n! \pi^{1/2}} \right)^{1/2}
$$

$$
\alpha = \left( \frac{mk}{\hbar^2} \right)^{1/4}
$$

$H_n(x)$ son los polinomios de Hermite

|$n$|$H_n(x)$|
|-|-|
|0|$1$|
|1|$2x$|
|2|$4x^2 -2$|
|3|$8x^3 - 12 x$|
|4|$16x^4 - 48 x^2 +12$|

**Importe las librerías numpy, math, pyplot de matplotlib y eval_hermite de scipy.special** 

# Importe librerías

import numpy as np
import math
from matplotlib import pyplot as plt
from scipy.special import eval_hermite

**Asigne valores a las constantes $m$, $k$, $\hbar$, $h$, $\alpha$. Considere $m=1$, $k=1$, $\hbar=1$.**

# Asigne valores

m=1.0
k=1.0
hbar=1.0
h = hbar*2*np.pi

alpha = np.power(m*k/(hbar**2),1.0/4.0)

**Defina un dominio de puntos para x, por ejemplo $1000$ puntos de $-5$ a $5$.**

# Defina dominio

x = np.linspace(-5,5,1000)

**Grafique la función de onda del oscilador Harmónico cuántico y su cuadrado para $n=0,2,4$.**

Compare con Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.

# Grafica

#Psi
for n in range(0,5,2):
    N_n = np.sqrt(alpha/(2**n*math.factorial(n)*np.pi**2))
    H_n = eval_hermite(n,alpha*x)
    psi_n = N_n*H_n*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_n,label=n)
plt.legend()
plt.title("$\psi$")
plt.show()

#Psi^2
for n in range(0,5,2):
    N_n = np.sqrt(alpha/(2.0**n*math.factorial(n)*np.pi**2))
    H_n = eval_hermite(n,alpha*x)
    psi_n = N_n*H_n*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_n**2,label=n)
plt.legend()
plt.title("$\psi^2$")
plt.show()

**Grafique la función de onda del oscilador Harmónico cuántico para $n=1,3,5$.**

Compare con Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.

# Gráfica

#Psi
for n in range(1,6,2):
    N_n = np.sqrt(alpha/(2**n*math.factorial(n)*np.pi**2))
    H_n = eval_hermite(n,alpha*x)
    psi_n = N_n*H_n*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_n,label=n)
plt.legend()
plt.title("$\psi$")
plt.show()

#Psi^2
for n in range(1,6,2):
    N_n = np.sqrt(alpha/(2**n*math.factorial(n)*np.pi**2))
    H_n = eval_hermite(n,alpha*x)
    psi_n = N_n*H_n*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_n**2,label=n)
plt.legend()
plt.title("$\psi^2$")
plt.show()

## Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.