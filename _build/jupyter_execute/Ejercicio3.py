# Oscilador Armónico Cuántico

El sistema es una partícula que se mueve en un potencial dado por
\begin{equation}
V = \frac{1}{2} kx^2
\end{equation}
donde $k$ es una constante. En mecánica clásica podemos encontrar sistemas con este potencial como los resortes, mientras que en química cúantica esto sirve para modelar la vibración de los enlaces.

En este caso, el Hamiltoniano contiene la energía cinética y el potencial
\begin{equation}
H = -\frac{\hbar^2}{2m} \frac{d}{dx} + \frac{1}{2} kx^2
\end{equation}

y la ecuación de Schrodinger a resolver es
\begin{equation}
-\frac{\hbar^2}{2m} \frac{d}{dx}\psi + \frac{1}{2} kx^2 \psi = E \psi
\end{equation}    

esta es una ecuación cuyas soluciones a la energía y la función de onda son
\begin{equation}
E_\nu = \left(\nu + \frac{1}{2} \right) \hbar \omega
\end{equation}    

\begin{equation}
\psi_\nu = N_\nu H_\nu(\alpha x) e^{-\alpha^2x^2/2}
\end{equation}    

aquí $\nu$ es un número cuántico tal que $\nu = 0,1,2,3,...$, y los términos $N_\nu$ y $\alpha$ están dados por
\begin{equation}
N_\nu = \left(\frac{\alpha}{2^\nu \nu! \pi^{1/2}} \right)^{1/2}
\end{equation}

\begin{equation}
\alpha = \left( \frac{mk}{\hbar^2} \right)^{1/4}
\end{equation}    

$H_\nu(x)$ son los polinomios de Hermite

|$\nu$|$H_\nu(x)$|
|-|-|
|0|$1$|
|1|$2x$|
|2|$4x^2 -2$|
|3|$8x^3 - 12 x$|
|4|$16x^4 - 48 x^2 +12$|

**Grafique la función de onda del oscilador Harmónico cuántico para $\nu=0,2,4$. Considere $m=1$, $k=1$, $\hbar=1$.**

Compare con Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.

from scipy.special import hermite as hermite
import numpy as np
from matplotlib import pyplot as plt

m=1.0
k=1.0
hbar=1.0

x = np.linspace(-5,5,100)

#Psi
for nu in range(0,5,2):
    alpha = np.power(m*k/(hbar**2),1.0/4.0)
    N_nu = np.sqrt(alpha/(np.power(2.0,nu)*np.math.factorial(nu)*np.pi**2))
    H_nu = hermite(nu)(alpha*x)
    psi_nu=N_nu*H_nu*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_nu,label=nu)
plt.legend()
plt.title("$\psi$")
plt.show()

#Psi^2
for nu in range(0,5,2):
    alpha = np.power(m*k/(hbar**2),1.0/4.0)
    N_nu = np.sqrt(alpha/(np.power(2.0,nu)*np.math.factorial(nu)*np.pi**2))
    H_nu = hermite(nu)(alpha*x)
    psi_nu=N_nu*H_nu*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_nu**2,label=nu)
plt.legend()
plt.title("$\psi^2$")
plt.show()

**Grafique la función de onda del oscilador Harmónico cuántico para $\nu=1,3,5$. Considere $m=1$, $k=1$, $\hbar=1$.**

Compare con Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.

from scipy.special import hermite as hermite
import numpy as np
from matplotlib import pyplot as plt

m=1.0
k=1.0
hbar=1.0

x = np.linspace(-5,5,100)

#Psi
for nu in range(1,6,2):
    alpha = np.power(m*k/(hbar**2),1.0/4.0)
    N_nu = np.sqrt(alpha/(np.power(2.0,nu)*np.math.factorial(nu)*np.pi**2))
    H_nu = hermite(nu)(alpha*x)
    psi_nu=N_nu*H_nu*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_nu,label=nu)
plt.legend()
plt.title("$\psi$")
plt.show()

#Psi^2
for nu in range(1,6,2):
    alpha = np.power(m*k/(hbar**2),1.0/4.0)
    N_nu = np.sqrt(alpha/(np.power(2.0,nu)*np.math.factorial(nu)*np.pi**2))
    H_nu = hermite(nu)(alpha*x)
    psi_nu=N_nu*H_nu*np.exp(-alpha**2*x**2/2)
    plt.plot(x,psi_nu**2,label=nu)
plt.legend()
plt.title("$\psi^2$")
plt.show()

## Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.