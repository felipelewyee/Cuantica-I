# Ejercicio 4

# Átomo de Hidrógeno

Sistema. Un electrón girándo en torno a un núcleo (protón) con radio no fijo.

El Hamiltoniano de este sistema contiene la energía cinética del electrón, la energía cinética del núcleo y la interacción coulómbica núcleo-electrón.

\begin{equation}
H=-\frac{\hbar^2}{2m_e}\nabla^2_e-\frac{\hbar^2}{2m_N}\nabla^2_N-\frac{e^2}{|\vec{r}_p-\vec{r}_e|}
\end{equation}

La ecuación de Schrodinger a resolver
\begin{equation}
\left(-\frac{\hbar^2}{2m_e}\nabla^2_e-\frac{\hbar^2}{2m_N}\nabla^2_N-\frac{e^2}{|\vec{r}_p-\vec{r}_e|}\right) \Psi = E \Psi
\end{equation}

Utilizamos coordenadas de masa reducida, y obtenemos:
\begin{equation}
\left(-\frac{\hbar^2}{2M_T}\nabla^2_{R_{cm}}-\frac{\hbar^2}{2\mu}\nabla^2_r-\frac{e^2}{|\vec{r}|}\right) \Psi = E \Psi
\end{equation}

La función de onda $\Psi$ de la ecuación anterior depende de $R$ y $r$. Se propone una solución por separación de variables, tal que $\Psi(R,r) = \Phi(R) \psi(r)$. Al sustituir se obtienen 2 ecuaciones:
\begin{equation}
\left(-\frac{\hbar^2}{2M}\nabla^2_{R_{cm}}\right) \Phi = E_R \Phi
\end{equation}
\begin{equation}
\left(-\frac{\hbar^2}{2\mu}\nabla^2_r-\frac{e^2}{|\vec{r}|}\right) \psi = E \psi
\end{equation}

La primera ecuación corresponde al movimiento de una partícula libre. Para resolver la segunda ecuación cambiamos a coordenadas polares, recordando que 
\begin{equation}
\nabla^2_r=\left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2} \right) = \left( \frac{1}{r} \frac{\partial}{\partial r^2} r + \frac{1}{r^2} \Lambda^2 \right) 
\end{equation}

Entonces
\begin{equation}
\left(-\frac{\hbar^2}{2\mu} \frac{1}{r} \frac{\partial}{\partial r^2} r -\frac{\hbar^2}{2\mu} \frac{1}{r^2} \Lambda^2 -\frac{e^2}{|\vec{r}|}\right) \psi = E \psi
\end{equation}

Se propone que la función de onda $\psi$ puede ser separada en una parte radial y una parte angular, es decir $\psi=R(r)Y(\theta,\phi)$. Esto genera la ecuación:

\begin{equation}
-\frac{\hbar^2}{2\mu} \frac{1}{r} \frac{d(rR)}{dr} -\frac{\hbar^2}{2\mu} \left[ \frac{e^2}{r} -\frac{l(l+1)}{r^2} \right] R = ER 
\end{equation}

Las soluciones a esta ecuación son:
\begin{equation}
R_{n,l}(r) = N_{n,l} \left( \frac{2r}{na_0} \right)^l e^{-r/na_0} L_{n+l}^{2l+1} \left( \frac{2r}{n a_0} \right)
\end{equation}

Donde
\begin{equation}
N_{n,l} = \left( \frac{2}{na_0} \right)^{3/2} \sqrt{\frac{(n-l-1)!}{2n[(n+l)!]^3}}
\end{equation}

Donde $L_{n+l}^{2l+1}$ son los polinomios asociados de Laguerre
\begin{equation}
L_{k}^N (r) = \frac{d^N}{dr^N} L_k(r)
\end{equation}
\begin{equation}
L_{k} (r) = e^r \frac{d^k}{dr^k} \left(r^k e^{-r}\right)
\end{equation}


En la tabla se muestran los primeros armónicos esféricos.


|$l$|$m_l$|Armónico esférico $Y_l^{m_l}(\theta,\phi)$|
|---|---|---|
|0|0|$\frac{1}{(4\pi)^{1/2}}$|
|1|-1|$+\frac{3}{(8\pi)^{1/2}} sin \theta e^{-i\phi}$|
|1|0|$\frac{3}{(4\pi)^{1/2}} cos \theta$|
|1|1|$-\frac{3}{(8\pi)^{1/2}} sin \theta e^{i\phi}$|
|2|-2|$\frac{15}{(32\pi)^{1/2}} sin^2 \theta e^{-2i\phi}$|
|2|-1|$+\frac{15}{(8\pi)^{1/2}} sin \theta cos \theta e^{-i\phi}$|
|2|0|$\frac{5}{(16\pi)^{1/2}} (3cos^2 \theta - 1)$|
|2|1|$-\frac{15}{(8\pi)^{1/2}} sin \theta cos \theta e^{i\phi}$|
|2|2|$\frac{15}{(32\pi)^{1/2}} sin^2 \theta e^{2i\phi}$|



**Realice la gráfica de R(r) (la parte radial de la función de onda) para los orbitales 1s (n=1, l=0),2s (n=2, l=0),3s (n=3, l=0) y 4s (n=4, l=0), y de $4\pi r^2R^2$.**

from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import derivative as derivative
from scipy.special import laguerre as laguerre

#Cambiar aqui para ajustar limites del eje X
r=np.linspace(0,25,100)

a0=1.0

#Cambiar aqui para usar s,p,d,etc
l=0
#Cambiar aquí para elegir limites de los numeros cuanticos n
n_min=1
n_max=4

for n in range(n_min,n_max+1):
    N=-np.power(2/(n*a0),3.0/2.0)*np.sqrt((np.math.factorial(n-l-1))/(2*n*np.power((np.math.factorial(n+l)),3.0)))
    L=laguerre(n+l)
    L=L/np.abs(L[n+l])
    assoc_L=derivative(L,2*r/(n*a0),n=2*l+1,order=2*l+3)

    R=N*(2*r/(n*a0))**l*np.exp(-r/(n*a0))*assoc_L

    plt.plot(r,R,label=n)

plt.legend()
#Cambiar aqui los titulos de los ejes
plt.xlabel("$r$")
plt.ylabel("$R$")
plt.title("$R(r)$ 1s, 2s ,3s y 4s")
plt.show()

for n in range(n_min,n_max+1):
    N=
    L=
    L=
    assoc_L=

    R=

    plt.plot(r,r**2.0*R**2.0,label=n)

plt.legend()
#Cambiar aqui los titulos de los ejes
plt.xlabel("$r$")
plt.ylabel("$4\pi r^2R^2$")
plt.title("$4\pi r^2R^2(r)$ 1s, 2s ,3s y 4s")
plt.show()

**Realice la gráfica de la parte radial de la función de onda para los orbitales 3s (n=3, l=0),3p (n=3, l=1) y 3d (n=3, l=2), y de $4\pi r^2R^2$.**

from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import derivative as derivative
from scipy.special import laguerre as laguerre

#Cambiar aqui para ajustar eje X
r=np.linspace(0,25,100)

a0=1.0

#Cambiar aqui para elegir n
n=3
#Cambiar aquí para elegir numeros cuanticos l
lmin=0
lmax=2

for l in range(lmin,lmax+1):
    N=-np.power(2/(n*a0),3.0/2.0)*np.sqrt((np.math.factorial(n-l-1))/(2*n*np.power((np.math.factorial(n+l)),3.0)))
    L=laguerre(n+l)
    L=L/np.abs(L[n+l])
    assoc_L=derivative(L,2*r/(n*a0),n=2*l+1,order=2*l+3)

    R=N*(2*r/(n*a0))**l*np.exp(-r/(n*a0))*assoc_L

    plt.plot(r,R,label=l)

plt.legend()
#Cambiar aqui los titulos de los ejes
plt.xlabel("$r$")
plt.ylabel("$R$")
plt.title("$R(r)$ 3s, 3p ,3d")
plt.show()

for l in range(lmin,lmax+1):
    N=
    L=
    L=
    assoc_L=

    R=

    plt.plot(r,r**2.0*R**2.0,label=l)

plt.legend()
#Cambiar aqui los titulos de los ejes
plt.xlabel("$r$")
plt.ylabel("$4\pi r^2R^2$")
plt.title("$4\pi r^2R^2$ 3s, 3p ,3d")
plt.show()

**A continuación se muestran los armónicos esféricos para $|Y_1^{-1}|$, $|Y_1^{0}|$, $|Y_1^{+1}|$ y su cuadrado.**

import scipy.special as sp
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Cambiar aqui para elegir l
l=1

#Cambiar aqui para elegir minimo de ml y maximo de ml en caso de que no se quieran todos los posibles
ml_min=-l
ml_max=l

theta = np.linspace(0,np.pi,100)
phi = np.linspace(0,2*np.pi,100)
THETA,PHI=np.meshgrid(theta,phi)

fig = plt.figure(figsize=plt.figaspect(0.5))

ml_n=0
    
for ml in range(ml_min,ml_max+1):
    
    ml_n=ml_n+1
    
    R=abs(sp.sph_harm(ml,l,PHI,THETA))

    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    ax = fig.add_subplot(2, l*2+1, ml_n, projection='3d')
    ax.plot_surface(X, Y, Z,cmap='YlOrRd')
    ax.set_xlim(-R.max(),R.max())
    ax.set_ylim(-R.max(),R.max())
    ax.set_zlim(-R.max(),R.max())    
    ax.set_title("$Y$"+" l="+str(l)+" ml="+str(ml))
    
    R=np.power(R,2.0)
    
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    ax = fig.add_subplot(2, l*2+1, ml_n+l*2+1, projection='3d')
    ax.plot_surface(X, Y, Z,cmap='YlOrRd')
    ax.set_xlim(-R.max(),R.max())
    ax.set_ylim(-R.max(),R.max())
    ax.set_zlim(-R.max(),R.max())
    ax.set_title("$Y^2$"+" l="+str(l)+" ml="+str(ml))

plt.show()

Recordemos que la función de onda es el producto de una parte radial y una parte angular $\psi=R(r)Y(\theta,\phi)$

Nombraremos a los orbitales p como:
\begin{eqnarray}
p_z &=& R_{n1} Y_{1}^{0} = R_{n1} \left(\frac{3}{4 \pi}\right)^{1/2} cos \theta \\
p_+ &=& R_{n1} Y_{1}^{+1} = -R_{n1} \left(\frac{3}{8 \pi}\right)^{1/2} sin \theta e^{i\phi} \\
p_- &=& R_{n1} Y_{1}^{-1} = R_{n1} \left(\frac{3}{8 \pi}\right)^{1/2} sin \theta e^{-i\phi}
\end{eqnarray}

El orbital $p_z$ es real, pero $p_+$ y $p_-$ son complejos. Hacemos la combinación lineal:
\begin{eqnarray}
p_x &=& \frac{1}{\sqrt{2}}(p_- - p_+) = R_{n1} \left(\frac{3}{4 \pi}\right)^{1/2} sin \theta cos \phi \\
p_y &=& \frac{i}{\sqrt{2}}(p_- + p_+) = R_{n1} \left(\frac{3}{4 \pi}\right)^{1/2} sin \theta sin \phi 
\end{eqnarray}

Represente la parte angular de las combinaciones lineales de $p_-$ y $p_+$ para formar $p_x$ y $p_y$.

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Seleccione numero cuantico l
l=1

theta = 
phi = 
THETA,PHI=np.meshgrid(theta,phi)

#Estas lineas seleccionan el harmonico esferico, en este ejemplo es ml=0, ml=-1, ml=1
Rz=sp.sph_harm(0,l,PHI,THETA)
R_m=sp.sph_harm(-1,l,PHI,THETA)
R_p=sp.sph_harm(1,l,PHI,THETA)

#Estas lineas hacen la combinacion lineal de harmonicos esfericos
Rx=(R_m-R_p)/np.sqrt(2)
Ry=1j*(R_m+R_p)/np.sqrt(2)

fig = plt.figure(figsize=plt.figaspect(0.5))

ml_n=0
    
for ml in range(-1,2):
    
    ml_n=ml_n+1

    #Estas lineas seleccionan la funion a graficar segun el ml    
    if(ml==-1):
        R=abs(Rx)
        orb_name="px"
    elif(ml==0):
        R=abs(Rz)
        orb_name="pz"        
    elif(ml==1):
        R=abs(Ry)
        orb_name="py"        
        
    
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    ax = fig.add_subplot(2, 2*l+1, ml_n, projection='3d')
    ax.plot_surface(X, Y, Z,cmap='YlOrRd')
    ax.set_xlim(-R.max(),R.max())
    ax.set_ylim(-R.max(),R.max())
    ax.set_zlim(-R.max(),R.max())
    ax.set_title("$Y$ "+orb_name)
    
    R=np.power(R,2.0)
    
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    ax = fig.add_subplot(2, 2*l+1, ml_n+(2*l+1), projection='3d')
    ax.plot_surface(X, Y, Z,cmap='YlOrRd')
    ax.set_xlim(-R.max(),R.max())
    ax.set_ylim(-R.max(),R.max())
    ax.set_zlim(-R.max(),R.max())
    ax.set_title("$Y^2$ "+orb_name)

plt.show()

A continuación se dan las expresiones de algunos orbitales

|Orbital|Función de onda|
|-------|---------------|
|$1s$|$N_1 e^{-r}$|
|$2s$|$N_2 (2-r)e^{-r/2}$|
|$2p_x$|$N_2 r sin\theta cos\phi e^{-r/2}$|
|$2p_y$|$N_2 r sin\theta sin\phi e^{-r/2}$|
|$2p_z$|$N_2 r cos\theta e^{-r/2}$|


**Haga la gráfica del orbital 1s para z=0**

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Dimensiones del dibujo
R=5.0

x=np.linspace(-5,5,100)
y=np.linspace(-5,5,100)
X,Y=np.meshgrid(x,y)

#Aqui se crea r, note que Z=0 (no hay Z)
r=np.sqrt(X**2.0+Y**2.0)

#Aqui va la formula del orbital
psi=np.exp(-r)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi,cmap='RdBu')
ax.set_title("$\Psi$")


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi**2,cmap='RdBu')
ax.set_title("$\Psi^2$")

plt.show()

**Haga la gráfica del orbital 2s para z=0**

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

R=5.0

x=np.linspace(-10,10,100)
y=np.linspace(-10,10,100)
X,Y=np.meshgrid(x,y)

#Aqui se crea r, note que Z=0 (no hay Z)
r=np.sqrt(X**2.0+Y**2.0)

#Aqui va la formula del orbital
psi=

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi,cmap='RdBu')
ax.set_title("$\Psi$")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi**2,cmap='RdBu')
ax.set_title("$\Psi^2$")

plt.show()

**Haga la gráfica del orbital $2p_z$ para z=0**

import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Dimensiones del dibujo
R=5.0

x=np.linspace(-15,15,100)
y=np.linspace(-15,15,100)

X,Y=np.meshgrid(x,y)

#Aqui se crea r, note que Z=0 (no hay Z)
r=np.sqrt(X**2.0+Y**2.0)

#Aqui va la formula del orbital
psi=np.exp(-r/2)*x


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi,cmap='YlOrRd')
ax.set_title("$\Psi$")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, psi**2,cmap='YlOrRd')
ax.set_title("$\Psi^2$")
    
plt.show()

# Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005.
- Pilar, F. L. Elementary Quantum Chemistry; 2001.
- Zettili, N. Quantum Mechanics: Concepts and Applications, 2nd ed.; Wiley: Chichester, U.K, 2009.
- Levine, I. N. Quantum Chemistry, 5th ed.; Prentice Hall: Upper Saddle River, N.J, 2000.
- McQuarrie, D. A.; Simon, J. D. Physical Chemistry: A Molecular Approach; University Science Books: Sausalito, Calif, 1997.

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.