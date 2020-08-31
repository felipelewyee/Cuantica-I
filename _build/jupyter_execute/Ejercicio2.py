# Parícula en el anillo, Partícula en la esfera y Rotor Rígido

## Partícula en el anillo

Es el sistema de una partícula moviéndose en una trayectoria de radio constante tal que $x^2 + y^2 = r^2$.

El Hamiltoniano para una partícula en dos dimensiones, tanto en coordenadas cartesianas como en coordenadas polares es
\begin{equation}
H = -\frac{\hbar^2}{2m} \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} \right) = -\frac{\hbar^2}{2mr^2} \left( \frac{\partial^2}{\partial r^2} + \frac{1}{r} \frac{\partial}{\partial r} + \frac{1}{r^2} \frac{\partial^2}{\partial \phi^2} \right)
\end{equation}

Como r es constante en el anillo, entonces se eliminan las derivadas respecto a r, y el Hamiltoniano se vuelve más simple
\begin{equation}
H = -\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2}
\end{equation}

Al sustituir el Hamiltoniano en la ecuación de Schrodinger, se obtiene
\begin{equation}
-\frac{\hbar^2}{2mr^2} \frac{d^2}{d\phi^2} \psi = E \psi
\end{equation}

La solución a la ecuación diferencial tiene la forma
\begin{equation}
\psi = Ae^{im_l\phi} + Be^{-im_l\phi}
\end{equation}

con $m_l = ( 2mr^2 E/\hbar^2 )^{1/2}$. Despejando se obtiene que 
\begin{equation}
E = \frac{\hbar^2 m_l^2}{2mr^2}
\end{equation}

Debido a que la función de onda debe ser contínua, y a que la partícula se mueve en un anillo, debe de cumplirse la condición cíclica $\psi(\phi) = \psi(\phi+2\pi)$, es decir que al dar una vuelta, la función de onda debe terminar en el mismo punto donde comenzó. Sustituyendo la función de onda en la condición cíclica se obtiene
\begin{equation}
Ae^{im_l\phi} + Be^{-im_l\phi} = Ae^{im_l\phi}e^{im_l2\pi} + Be^{-im_l\phi}e^{-im_l2\pi}
\end{equation}

Para que la igualdad anterior pueda cumplirse, $m_l$ debe ser un número entero, tal que se cumpla $e^{im_l2\pi}=e^{-im_l2\pi}=1$. Esto origina la cuantización $m_l = {0, \pm 1, \pm 2, \cdots}$.

Tras normalizar y con B=0,
\begin{equation}
\psi = \left( \frac{1}{2\pi} \right)^{1/2} e^{i m_l \phi} = \left( \frac{1}{2\pi} \right)^{1/2} cos(m_l \phi) + i \left( \frac{1}{2\pi} \right)^{1/2} sin(m_l \phi)
\end{equation}

**Grafique la función de onda y su cuadrado para $m_l=1$**

# Grafica

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

#Numero cuantico
ml=1

#Coordenadas polares
phi = np.linspace(0, 2 * np.pi, 100)
r=1.0
psi_r = np.sqrt(1/(2*np.pi))*np.cos(ml*phi)
psi_i = np.sqrt(1/(2*np.pi))*np.sin(ml*phi)

#Grafica de la funcion de onda
fig = plt.figure()
ax = fig.gca(projection='3d')
x = r * np.sin(phi)
y = r * np.cos(phi)
ax.plot(x, y, 0, color='k') #Eje de la grafica
ax.plot(x, y, psi_r, label='psi-R')
ax.plot(x, y, psi_i, label='psi-I')
ax.legend()
plt.show()

#Grafica del cuadrado de la funcion de onda
fig = plt.figure()
ax = fig.gca(projection='3d')
x = r * np.sin(phi)
y = r * np.cos(phi)
ax.plot(x, y, 0, color='k') #Eje de la grafica
ax.plot(x, y, ((psi_r+1J*psi_i)*(psi_r-1J*psi_i)).real, label='$psi^2$')
ax.legend()
plt.show()

## Particula en la esfera

Se tiene una partícula moviéndose sobre una superficie esférica de radio constante.

La ecuación de Schrodinger a resolver es
\begin{equation}
-\frac{\hbar^2}{2m} \nabla^2 \psi(\theta,\phi) = E \psi(\theta,\phi)
\end{equation}

Donde
\begin{equation}
\nabla^2 = \frac{1}{r} \frac{\partial^2}{\partial r^2}r + \frac{1}{r^2} \Lambda^2 = \frac{1}{r} \frac{\partial^2}{\partial r^2} r + \frac{1}{r^2} \left( \frac{1}{sin^2 \theta} \frac{\partial^2}{\partial \phi^2} + \frac{1}{sin \theta} \frac{\partial}{\partial \theta} sin \theta \frac{\partial}{\partial \theta} \right)
\end{equation}

Si r es constante, entonces
\begin{equation}
-\frac{\hbar^2}{2mr^2} \Lambda^2 \psi(\theta,\phi) = E \psi(\theta,\phi)
\end{equation}

Las soluciones de esta ecuación son los armónicos esféricos, dados por

\begin{equation}
Y_l^{m_l}(\theta,\phi) = \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m_l|)!}{(l+|m_l|)!}} e^{i m_l \phi} P_l^{m_l}(cos(\theta))
\end{equation}

donde $P_l^{m_l}(cos(\phi))$ son los polinomios asociados de Legendre, dados por
\begin{equation}
P_l^{m_l}(x) = \frac{l}{2^l l!}(1-x^2)^{|m_l|/2} \frac{d^{l+|m_l|}}{dx^{l+|m_l|}} (x^2-1)^l
\end{equation}

En la tabla se muestra la forma de los primeros armónicos esféricos. En este punto han aparecido dos números cuánticos, tal que $l = 0,1,2,3,...$ y $m_l = -l, -l+1, 0, l-1, l$


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

La energía del sistema está dada por
\begin{equation}
E = -l(l+1)
\end{equation}

**Grafique el armónico esférico $|Y_1^{0}|$ y su cuadrado.**

#Grafica

import scipy.special as sp
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

ml=0
l=1

theta = np.linspace(0,np.pi,100)
phi = np.linspace(0,2*np.pi,100)
THETA,PHI=np.meshgrid(theta,phi)

R=np.abs(sp.sph_harm(ml,l,PHI,THETA))

X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)

fig = plt.figure()
ax = plt.gca(projection='3d')
ax.plot_surface(X, Y, Z,cmap='YlOrRd')
#ax.set_xlim(-0.4,0.4)
#ax.set_ylim(-0.4,0.4)
#ax.set_zlim(-0.4,0.4)
ax.set_title("$\Psi$")
plt.show()

R=np.power(R,2.0)

X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)

fig = plt.figure()
ax = plt.gca(projection='3d')
ax.plot_surface(X, Y, Z,cmap='YlOrRd')
ax.set_xlim(-0.2,0.2)
ax.set_ylim(-0.2,0.2)
ax.set_zlim(-0.2,0.2)
ax.set_title("$\Psi^2$")
plt.show()

## Rotor Rígido

Este sistema consiste en dos partículas de masa $m_1$ y $m_2$ moviéndose con una separación constante $r=r_2-r_1$.

El Hamiltoniano consta de los términos de energía cinética de ambas partículas
\begin{equation}
H = -\frac{\hbar^2}{2m_1} \nabla^2_1 -\frac{\hbar^2}{2m_2} \nabla^2_2
\end{equation}

Este sistema es equivalente al de una partícula de masa reducida girando en torno al centro de masa de ambas partículas. La masa total del sistema está dada por la suma de las masas de las partículas
\begin{equation}
m_T = m_1 + m_2
\end{equation}

La masa reducida de la nueva partícula está dada por
\begin{equation}
\frac{1}{\mu} = \frac{1}{m_1} + \frac{1}{m_2}
\end{equation}

El centro de masa del sistema se calcula mediante
\begin{equation}
R_{cm} = \left( \frac{m_1}{m_T} \right) r_1 + \left( \frac{m_2}{m_T} \right) r_2
\end{equation}

y el Hamiltoniano se calcula con la energía cinética de la masa reducida y del centro de masa, es decir
\begin{equation}
H = -\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}
\end{equation}

La ecuación de Schrodinger a resolver es
\begin{equation}
\left(-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}\right) \psi = E \psi
\end{equation}

Se propone que la función de onda se puede separar en el producto de una función de onda del centro de masa y una función de onda de la partícula de masa reducida
\begin{equation}
\psi=\psi_{cm}\psi_r
\end{equation}

Al sustituir en la ecuación de Schrodinger se obtiene
\begin{equation}
\left(-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} - \frac{\hbar^2}{2\mu} \nabla^2_{r}\right) \psi_{cm}\psi_r = E \psi_{cm}\psi_r
\end{equation}

Si consideramos que la energía está dada por $E_T = E_{cm} + E_{r}$ y distribuímos, resulta

\begin{equation}
-\psi_{r} \frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} \psi_{cm} - \psi_{cm}\frac{\hbar^2}{2\mu} \nabla^2_{r} \psi_r = \psi_r E_{cm} \psi_{cm} +  \psi_{cm} E_r \psi_r
\end{equation}

Si multiplicamos ambos lados de la ecuación anterior por $\frac{1}{\psi_{r}\psi_{cm}}$, resulta

\begin{equation}
-\frac{1}{\psi_{cm}} \left( \frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} \psi_{cm} - E_{cm} \psi_{cm} \right) = \frac{1}{\psi_{r}} \left( \frac{\hbar^2}{2\mu} \nabla^2_{r} \psi_r + E_r \psi_r \right)
\end{equation}

ya que el lado izquierdo solo depende de las coordenadas del centro de masa, y el lado derecho solo depende de las coordenadas de la masa reducida, ambos lados deben ser igual a una constante. Si elegimos esta constante como cero, y despejamos lo que está dentro de cada paréntesis se obtienen dos ecuaciones independientes. La primera ecuación corresponde al movimiento del centro de masa del sistema y la hemos estudiado previamene en el movimiento de la partícula libre
\begin{equation}
-\frac{\hbar^2}{2m_T} \nabla^2_{R_{cm}} \psi_{cm} = E_{cm} \psi_{cm}
\end{equation}
esta ecuación tiene como soluciones
\begin{equation}
\psi_{cm} = A e^{ikx} + B e^{-ikx}
\end{equation}
con $k^2=2m_TE/\hbar^2$, y simplemente nos dice que el sistema en conjunto se mueve libremente por el espacio.

La segunda ecuación corresponde a la masa reducida, y la hemos estudiado en la partícula en la esfera
\begin{equation}
-\frac{\hbar^2}{2\mu} \nabla^2_{r} \psi_r = E_r \psi_r
\end{equation}
sabemos por tanto que su solución son los armónico esféricos $Y_l^{m_l}(\theta,\phi)$ con $E = -l(l+1)$.

## Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005.
- Zettili, N. Quantum Mechanics: Concepts and Applications, 2nd ed.; Wiley: Chichester, U.K, 2009.
- Levine, I. N. Quantum Chemistry, 5th ed.; Prentice Hall: Upper Saddle River, N.J, 2000.
- McQuarrie, D. A.; Simon, J. D. Physical Chemistry: A Molecular Approach; University Science Books: Sausalito, Calif, 1997.