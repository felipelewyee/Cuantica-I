# Ejercicio 3

## Serie de Taylor

Sea una función f(x) continua e infinitamente diferenciable, esta puede expresarse en torno a $a$ mediante una serie de potencias
\begin{equation}
f(x) = \sum_{n=0}^\infty c_n (x-a)^n
\end{equation}

La m-ésima derivada de $f(x)$ es:
\begin{equation}
f^{(m)}(x) = \sum_{n=m}^\infty c_n (m!)(x-a)^{n-m} = \sum_{p=0}^\infty c_{p+m} (m!)(x-a)^{p}
\end{equation}

donde se ha hecho el cambio de índice $p=n-m$. Al evaluar en x=a se obtiene
\begin{equation}
f^{(m)}(a) = \sum_{p=0}^\infty c_{p+m} (m!)(a-a)^{p} = c_p (m!)
\end{equation}

El único término que sobrevive es el que tiene potencia cero, despejando $c_p$:
\begin{equation}
c_m = \frac{f^{(m)}(a)}{m!}
\end{equation}

Sustituyendo en la serie
\begin{equation}
f(x) = \sum_{n=0}^\infty \frac{f^{(n)}(a)}{n!} (x-a)^n
\end{equation}

**Haga la expansión en series de Taylor en torno a $x=0$ de $f(x)=e^{x}$ en el intervalo $x \varepsilon [0,10]$ hasta potencias de grado 10.**



## Serie de Fourier

Una funcion periodica se puede aproximar por
\begin{equation}
f(x)=\frac{a_0}{2} + \sum_{n=1}^{\infty} \left[ a_n cos \left( \frac{2n\pi}{T}x \right) + b_n sin \left( \frac{2n\pi}{T}x \right) \right]
\end{equation}

donde
\begin{equation}
a_0 = \frac{2}{T} \int\limits_{-T/2}^{T/2} f(x) dx
\end{equation}
\begin{equation}
a_n = \frac{2}{T} \int\limits_{-T/2}^{T/2} f(x) cos \left( \frac{2n\pi}{T}x \right) dx
\end{equation}
\begin{equation}
b_n = \frac{2}{T} \int\limits_{-T/2}^{T/2} f(x) sin \left( \frac{2n\pi}{T}x \right) dx
\end{equation}


**Haga la expansión en series de Fourier de la $f(x)=x$ en el intervalo $x \varepsilon [-10,10]$ con $n=10$**



**Haga la expansión en series de Fourier de la fución escalón en el intervalo $x \varepsilon [-10,10]$ con $n=10$**

\begin{equation}
f(x) = \left\{
  \begin{array}{ll}
  0      & \mathrm{si\ } x < 0\\
  1 & \mathrm{si\ } x>0 \\
  \end{array}
  \right.
\end{equation}




## Oscilador Armónico Cuántico

El sistema es una partícula en un potencial
\begin{equation}
V = \frac{1}{2} kx^2
\end{equation}
donde $k$ es una constante.

En este caso, el Hamiltoniano es:
\begin{equation}
H = -\frac{\hbar^2}{2m} \frac{d}{dx} + \frac{1}{2} kx^2
\end{equation}

y la ecuación de Schrodinger a resolver es:
En este caso, el Hamiltoniano es:
\begin{equation}
-\frac{\hbar^2}{2m} \frac{d}{dx}\psi + \frac{1}{2} kx^2 \psi = E \psi
\end{equation}    

Las soluciones a este problema son:
\begin{equation}
E_\nu = \left(\nu + \frac{1}{2} \right) \hbar \omega
\end{equation}    

\begin{equation}
\psi_\nu = N_\nu H_\nu(\alpha x) e^{-\alpha^2x^2/2}
\end{equation}    

Con:

\begin{equation}
N_\nu = \left(\frac{\alpha}{2^\nu \nu! \pi^{1/2}} \right)^{1/2}
\end{equation}

\begin{equation}
\alpha = \left( \frac{mk}{\hbar^2} \right)^{1/4}
\end{equation}    

y $H_\nu(x)$ los polinomios de Hermite

|$\nu$|$H_\nu(x)$|
|-|-|
|0|$1$|
|1|$2x$|
|2|$4x^2 -2$|
|3|$8x^3 - 12 x$|
|4|$16x^4 - 48 x^2 +12$|

**Grafique la función de onda del oscilador Harmónico cuántico para $\nu=0,2,4$. Considere $m=1$, $k=1$, $\hbar=1$.**

Compare con Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.



## Referencias

- Gersting, J. L. Technical Calculus with Analytic Geometry, Dover ed.; Dover: New York, 1992.
- Jackson, J. D. Mathematics for Quantum Mechanics: An Introductory Survey of Operators, Eigenvalues, and Linear Vector Spaces, Dover ed.; Dover books on mathematics; Dover: Mineola, N.Y, 2006.
- Tolstov, G. P. Fourier Series, Nachdr.; Dover books on mathematics; Dover: New York, 2009.
- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005, p.62.

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.