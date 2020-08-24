# Ejercicio 2

# Partícula en el anillo

Es el sistema de una partícula moviéndose en una trayectoria de radio constante tal que $x^2 + y^2 = r^2$.

El Hamiltoniano para una partícula en dos dimensiones es:
\begin{equation}
H = -\frac{\hbar^2}{2m} \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right) = -\frac{\hbar^2}{2mr^2} \left( \frac{\partial^2}{\partial r^2} + \frac{1}{r} \frac{\partial}{\partial r} + \frac{1}{r^2} \frac{\partial^2}{\partial \phi^2} \right)
\end{equation}


Como r es constante, entonces:
\begin{equation}
H = -\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2}
\end{equation}

Tal que la ecuación de Schrodinger resulta:
\begin{equation}
-\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2} \psi = E \psi
\end{equation}

Cuyas soluciones son
\begin{equation}
\psi = Ae^{im_l\phi} + Be^{-im_l\phi}
\end{equation}

con $m_l = ( 2mr^2 E/\hbar^2 )^{1/2}$, de donde se obtiene que 
\begin{equation}
E = \frac{\hbar^2 m_l^2}{2mr^2}
\end{equation}

La condición cíclica implica que $\psi(\phi) = \psi(\phi+2\pi)$, es decir:
\begin{equation}
Ae^{im_l\phi} + Be^{-im_l\phi} = Ae^{im_l\phi}e^{im_l2\pi} + Be^{-im_l2\pi}
\end{equation}

Por la cuál $m_l$ debe ser un número entero. Es decir: $m_l = {0, \pm 1, \pm 2, \cdots}$.

Tras normslizar y con B=0,
\begin{equation}
\psi = \left( \frac{1}{2\pi} \right)^{1/2} e^{i m_l \phi} = \left( \frac{1}{2\pi} \right)^{1/2} cos(m_l \phi) + i \left( \frac{1}{2\pi} \right)^{1/2} sin(m_l \phi)
\end{equation}

**Grafique la función de onda y su cuadrado para $m_l=1$**



# Particula en la esfera

Se tiene una partícula moviéndose sobre una superficie esférica de radio constante.

La ecuación de Schrodinger a resolver es:
\begin{equation}
-\frac{\hbar^2}{2m} \nabla^2 \psi(\theta,\phi) = E \psi(\theta,\phi)
\end{equation}

Donde
\begin{equation}
\nabla^2 = \frac{1}{r} \frac{\partial^2}{\partial r^2}r + \frac{1}{r^2} \Lambda^2 = \frac{1}{r} \frac{\partial^2}{\partial r^2} r + \frac{1}{r^2} \left( \frac{1}{sin^2 \theta} \frac{\partial^2}{\partial \phi^2} + \frac{1}{sin \theta} \frac{\partial}{\partial \theta} sin \theta \frac{\partial}{\partial \theta} \right)
\end{equation}

Si r es constante entonces
\begin{equation}
-\frac{\hbar^2}{2mr^2} \Lambda^2 \psi(\theta,\phi) = E \psi(\theta,\phi)
\end{equation}

Las soluciones de esta ecuación son los armónicos esféricos.

\begin{equation}
Y_l^{m_l}(\theta,\phi) = \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m_l|)!}{(l+|m_l|)!}} e^{i m_l \phi} P_l^{m_l}(cos(\theta))
\end{equation}

donde $P_l^{m_l}(cos(\phi))$ son los polinomios asociados de Legendre dados por:
\begin{equation}
P_l^{m_l}(x) = \frac{l}{2^l l!}(1-x^2)^{|m_l|/2} \frac{d^{l+|m_l|}}{dx^{l+|m_l|}} (x^2-1)^l
\end{equation}

En la tabla se muestran los primeros armónicos esféricos.


|$l$|$m_l$|Armónico esférico $Y_l^{m_l}(\theta,\phi)$|
|---|---|---|
|0|0|$\frac{1}{(4\pi)^{1/2}}$|
|1|-1|$+\frac{3}{(8\pi)^{1/2}} sin \theta e^{-i\phi}$|
|1|0|$\frac{3}{(4\pi)^{1/2}} cos \theta$|
|1|1|$-\frac{3}{(8\pi)^{1/2}} sin \theta e^{i\phi}$|
|2|-2|$+\frac{15}{(32\pi)^{1/2}} sin^2 \theta e^{-2i\phi}$|
|2|-1|$+\frac{15}{(8\pi)^{1/2}} sin \theta cos \theta e^{-i\phi}$|
|2|0|$\frac{5}{(16\pi)^{1/2}} (3cos^2 \theta - 1)$|
|2|1|$-\frac{15}{(8\pi)^{1/2}} sin \theta cos \theta e^{i\phi}$|
|2|2|$-\frac{15}{(32\pi)^{1/2}} sin^2 \theta e^{2i\phi}$|



**Grafique el armónico esférico $|Y_1^{0}|$ y su cuadrado.**



# Rotor Rígido

Se trata de un sistema de dos masas $m_1$ y $m_2$ moviéndose a una separación fija $r=r_1-r_2$. El Hamiltoniano consta de los términos de energía cinética:
\begin{equation}
H = -\frac{\hbar^2}{2m_1} \nabla^2_1 -\frac{\hbar^2}{2m_2} \nabla^2_2
\end{equation}

Este sistema es equivalente al de una partícula de masa reducida girando en torno al centro de masa de ambas partículas, con
\begin{equation}
m_T = m_1 + m_2
\end{equation}
\begin{equation}
\frac{1}{\mu} = \frac{1}{m_1} + \frac{1}{m_2}
\end{equation}
\begin{equation}
R_{cm} = \left( \frac{m_1}{m_T} \right) r_1 + \left( \frac{m_2}{m_T} \right) r_2
\end{equation}

y los Hamiltonianos son equivalentes, es decir:
\begin{equation}
H = -\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}
\end{equation}

La ecuación de Schrodinger a resolver es:
\begin{equation}
\left(-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}\right) \psi = E \psi
\end{equation}

Se propone una función de onda separable en una función de onda del centro de masa y una función de onda de la partícula de masa reducida $\psi=\psi_{cm}\psi_r$:
\begin{equation}
\left(-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}\right) \psi_{cm}\psi_r = E \psi_{cm}\psi_r
\end{equation}

\begin{equation}
-\psi_{r} \frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} \psi_{cm} - \psi_{cm}\frac{\hbar^2}{2\mu} \nabla^2_{r} \psi_r = \psi_r E_{cm} \psi_{cm} +  \psi_{cm} E_r \psi_r
\end{equation}

Separando los términos de $\psi_r$ y $\psi_{cm}$, y considrando r constante:

**Partícula libre**
\begin{equation}
-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} \psi_{cm} = E_{cm} \psi_{cm}
\end{equation}
**Partícula en la esfera**
\begin{equation}
- \frac{\hbar^2}{2\mu} \nabla^2_{r} \psi_r = E_r \psi_r
\end{equation}

# Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005.
- Zettili, N. Quantum Mechanics: Concepts and Applications, 2nd ed.; Wiley: Chichester, U.K, 2009.
- Levine, I. N. Quantum Chemistry, 5th ed.; Prentice Hall: Upper Saddle River, N.J, 2000.
- McQuarrie, D. A.; Simon, J. D. Physical Chemistry: A Molecular Approach; University Science Books: Sausalito, Calif, 1997.

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.