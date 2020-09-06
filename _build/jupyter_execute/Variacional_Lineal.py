# Método variacional lineal

El método variacional lineal permite resolver el problema

$$
H\psi = E \psi
$$

al expresar la función de onda como una combinación lineal de funciones. Cuando las funciones $\psi_i$ forman un conjunto completo, la combinación lineal es exacta, sin embargo, en la práctica se usan solo algunas funciones, por lo que la función de onda resultante es aproximada, es decir

$$
\psi_{prueba} = \sum_{i=1} c_i \psi_i  
$$

Al sustituir la expansión de la función de onda en la ecuación de Schrodinger, se obtiene la ecuación

$$
\mathcal{H} \mathcal{C} = \mathcal{S} \mathcal{C} \mathcal{\varepsilon}
$$

En general uno raliza los siguientes pasos:

1. Se seleccionan las funciones $\psi_i$ que se usarán para expandir la función de onda, es común elegir funciones exponenciales o gaussianas.
2. Se evalúan las matrices $\mathcal{H}$, y $\mathcal{S}$.
3. Se resuelve el problema de valores propios $\mathcal{H}\mathcal{C} = \mathcal{S}\mathcal{C} \mathcal{\epsilon}$.
4. Se construye la función de onda utilizando los coeficientes obtenidos.

Una consecuencia de este método es que sin importar las funciones $\psi_i$ que se usen, siempre que cumplan con las restricciones del problema, la solución de prueba siempre tiene una energía mayor o igual a la solución exacta. Por lo tanto, podemos construir varias funciones de prueba y tomar la que de la energía más baja.

Para ejemplificar este procedimiento, resolveremos el átomo de hidrógeno utilizando el método variacional lineal. Ya que sí conocemos las soluciones exactas del átomo de hidrógeno, podremos comparar las funciones de onda aproximadas con las funciones de onda exactas.

**Importe las siguientes librerías**

- numpy
- sympy

# Librerías

import numpy as np
import sympy as sp

## Átomo de hidrógeno con dos gausianas

**Paso 1.** Seleccionar funciones $\psi_i$ para construir $\psi_{prueba}$

$$
\psi_{prueba} = \sum_{i=1} c_i \psi_i  
$$

Para este ejemplo tomaremos dos funciones gaussianas:

$$
\psi_1 = \left( \frac{2(1.309756377)}{\pi} \right)^{\frac{3}{4}} e^{-1.309756377r^2}
$$

$$
\psi_2 = \left( \frac{2(0.233135974)}{\pi} \right)^{\frac{3}{4}} e^{-0.233135974r^2}
$$

**Defina las funciones gaussianas usando álgebra simbólica**

# Defina funciones



r = sp.Symbol("r")
psi_1 = (2*1.309756377/sp.pi)**(3/4)*sp.exp(-1.309756377*r**2)
psi_2 = (2*0.233135974/sp.pi)**(3/4)*sp.exp(-0.233135974*r**2)
print("psi_1")
sp.pprint(psi_1)
print("psi_2")
sp.pprint(psi_2)

**Paso 2.** Evaluar las matrices $\mathcal{H}$ y $\mathcal{S}$.

El Hamiltoniano del átomo de hidrógeno es

$$
\hat{H} = -\frac{1}{2} \nabla^2 - \frac{1}{r} = -\frac{1}{2} \frac{1}{r} \frac{d^2}{dr^2} r -\frac{1}{r}
$$

Por tanto, la matriz $\mathcal{H}$ es

$$
\mathcal{H} = \begin{pmatrix} H_{11} & H_{12}\\
H_{21} & H_{22}
\end{pmatrix} = \begin{pmatrix} 4 \pi\int_0^\infty \psi_1^* \hat{H} \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_1^* \hat{H} \psi_2 r^2dr\\
4 \pi\int_0^\infty \psi_2^* \hat{H} \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_2^* \hat{H} \psi_2 r^2dr
\end{pmatrix}
$$

**Genere una matriz de 2x2 y evalúe las integrales**

# Matriz H

H_1 = -1/2*1/r*sp.diff(r*psi_1,r,r) - 1/r*psi_1
H_2 = -1/2*1/r*sp.diff(r*psi_2,r,r) - 1/r*psi_2

H = sp.zeros(2)
H[0,0] = 4*sp.pi*sp.integrate(psi_1*H_1*r**2,(r,0,sp.oo))
H[0,1] = 4*sp.pi*sp.integrate(psi_1*H_2*r**2,(r,0,sp.oo))
H[1,0] = 4*sp.pi*sp.integrate(psi_2*H_1*r**2,(r,0,sp.oo))
H[1,1] = 4*sp.pi*sp.integrate(psi_2*H_2*r**2,(r,0,sp.oo))
H=H.evalf()
H

La matriz S es

$$
\mathcal{S} = \begin{pmatrix} S_{11} & S_{12}\\
S_{21} & S_{22}
\end{pmatrix} = \begin{pmatrix} 4 \pi\int_0^\infty \psi_1^* \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_1^* \psi_2 r^2dr\\
4 \pi\int_0^\infty \psi_2^* \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_2^* \psi_2 r^2dr
\end{pmatrix}
$$

**Genere una matrix de 2x2 y evalúe las integrales**

# Matriz S

S = sp.zeros(2)
S[0,0] = 4*sp.pi*sp.integrate(psi_1*psi_1*r**2,(r,0,sp.oo))
S[0,1] = 4*sp.pi*sp.integrate(psi_1*psi_2*r**2,(r,0,sp.oo))
S[1,0] = 4*sp.pi*sp.integrate(psi_2*psi_1*r**2,(r,0,sp.oo))
S[1,1] = 4*sp.pi*sp.integrate(psi_2*psi_2*r**2,(r,0,sp.oo))
S

**Paso 3.** Resolver $\color{red}{\mathcal{H}\mathcal{C} = \mathcal{S}\mathcal{C} \mathcal{\epsilon}}$

Para ello utilizaremos la instrucción
~~~python
E,C = LA.eigh(H,S)
~~~
la cual resuelve directamente el problema $\mathcal{H}\mathcal{C} = \mathcal{S}\mathcal{C} \mathcal{\epsilon}$. La columna de $\mathcal{C}$ con la energía más baja nos indica los coeficientes de la combinación lineal de la función de onda.

# Resuelva HC = SCe

from scipy import linalg as LA
import numpy

H = np.array(H).astype(np.float64)
S = np.array(S).astype(np.float64)

E,C = LA.eigh(H,S)
print(E)
print(C)

**Paso 4.** Susituir los coeficientes en $\psi_{prueba} = \sum_{i=1} c_i \psi_i$ para construir la función de onda.

# Genere función de prueba

psi_p = C[0][0]*psi_1 + C[1][0]*psi_2
psi_p

Adicionalmente, se puede comprobar que se cumple la condición de normalización. **Evalúe la integral**

$$
4\pi \int_0^{\infty} r^2 |\psi_{prueba}|^2 dr= 1
$$

# Integral

4*sp.pi*sp.integrate(psi_p*psi_p*r**2,(r,0,sp.oo))

**Almacene la función de prueba, así como cada una de las gausianas que la componen para comparar con la función de onda exacta**.

# Código

psi_2g = psi_p
psi_2g_1 = psi_1
psi_2g_2 = psi_2

**Genere la gráfica de las dos funciones $\psi_1$ y $\psi_2$, así como la $\psi_{prueba}$ y la solución exacta para el átomo de hidrógeno (1s)**

```{tip}
Recuerde

$$
1s = \pi^{-1/2} e^{-|r|}
$$
```

# Gráfica

from matplotlib import pyplot as plt

s=1/sp.sqrt(sp.pi)*sp.exp(-sp.Abs(r))
lam_s = sp.lambdify(r,s,modules=['numpy'])
lam_psi_2g = sp.lambdify(r,psi_2g,modules=['numpy'])
lam_psi_2g_1 = sp.lambdify(r,psi_2g_1,modules=['numpy'])
lam_psi_2g_2 = sp.lambdify(r,psi_2g_2,modules=['numpy'])

r1 = np.linspace(-7,7,100)
psi_s = lam_s(r1)
psi_2g1 = lam_psi_2g(r1)
psi_2g_1 = lam_psi_2g_1(r1)
psi_2g_2 = lam_psi_2g_2(r1)

plt.plot(r1,psi_s,label="1s")
plt.plot(r1,-psi_2g1,label="2g")
plt.plot(r1,psi_2g_1,label="psi_1",linestyle=':')
plt.plot(r1,psi_2g_2,label="psi_2",linestyle=':')
plt.legend()
plt.show()

## Átomo de hidrógeno con tres gausianas

**Paso 1.** Seleccionar funciones $\psi_i$ para construir $\psi_{prueba}$

$$
\psi_{prueba} = \sum_{i=1} c_i \psi_i  
$$

Para este ejemplo **declare tres funciones gaussianas**

$$
\psi_1 = \left( \frac{2(3.42525091)}{\pi} \right)^{\frac{3}{4}} e^{-3.42525091r^2}
$$

$$
\psi_2 = \left( \frac{2(0.62391373)}{\pi} \right)^{\frac{3}{4}} e^{-0.62391373r^2}
$$

$$
\psi_3 = \left( \frac{2(0.16885540)}{\pi} \right)^{\frac{3}{4}} e^{-0.16885540r^2}
$$

# Declare funciones

sp.init_printing()

r = sp.Symbol("r")
psi_1 = (2*3.42525091/sp.pi)**(3/4)*sp.exp(-3.42525091*r**2)
psi_2 = (2*0.62391373/sp.pi)**(3/4)*sp.exp(-0.62391373*r**2)
psi_3 = (2*0.16885540/sp.pi)**(3/4)*sp.exp(-0.16885540*r**2)
print("psi_1")
sp.pprint(psi_1)
print("psi_2")
sp.pprint(psi_2)
print("psi_3")
sp.pprint(psi_3)

El Hamiltoniano del átomo de hidrógeno es

$$
\hat{H} = -\frac{1}{2} \nabla^2 - \frac{1}{r} = -\frac{1}{2} \frac{1}{r} \frac{d^2}{dr^2} r -\frac{1}{r}
$$

Por tanto, la matriz $\mathcal{H}$ es

$$
\mathcal{H} = \begin{pmatrix} H_{11} & H_{12} & H_{13}\\
H_{21} & H_{22} & H_{23}\\
H_{31} & H_{32} & H_{33}\\
\end{pmatrix} = \begin{pmatrix} 4 \pi\int_0^\infty \psi_1^* \hat{H} \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_1^* \hat{H} \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_1^* \hat{H} \psi_3 r^2dr\\
4 \pi\int_0^\infty \psi_2^* \hat{H} \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_2^* \hat{H} \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_2^* \hat{H} \psi_3 r^2dr\\
4 \pi\int_0^\infty \psi_3^* \hat{H} \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_3^* \hat{H} \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_3^* \hat{H} \psi_3 r^2dr\\
\end{pmatrix}
$$

**Evalúe la matriz H**

# Matriz H

H_1 = -1/2*1/r*sp.diff(r*psi_1,r,r) - 1/r*psi_1
H_2 = -1/2*1/r*sp.diff(r*psi_2,r,r) - 1/r*psi_2
H_3 = -1/2*1/r*sp.diff(r*psi_3,r,r) - 1/r*psi_3

H = sp.zeros(3)
H[0,0] = 4*sp.pi*sp.integrate(psi_1*H_1*r**2,(r,0,sp.oo))
H[0,1] = 4*sp.pi*sp.integrate(psi_1*H_2*r**2,(r,0,sp.oo))
H[0,2] = 4*sp.pi*sp.integrate(psi_1*H_3*r**2,(r,0,sp.oo))
H[1,0] = 4*sp.pi*sp.integrate(psi_2*H_1*r**2,(r,0,sp.oo))
H[1,1] = 4*sp.pi*sp.integrate(psi_2*H_2*r**2,(r,0,sp.oo))
H[1,2] = 4*sp.pi*sp.integrate(psi_2*H_3*r**2,(r,0,sp.oo))
H[2,0] = 4*sp.pi*sp.integrate(psi_3*H_1*r**2,(r,0,sp.oo))
H[2,1] = 4*sp.pi*sp.integrate(psi_3*H_2*r**2,(r,0,sp.oo))
H[2,2] = 4*sp.pi*sp.integrate(psi_3*H_3*r**2,(r,0,sp.oo))

H=H.evalf()
H

**Evalúe la matriz S**

$$
\mathcal{S} = \begin{pmatrix} S_{11} & S_{12} & S_{13}\\
S_{21} & S_{22} & S_{23} \\
S_{31} & S_{32} & S_{33} \\
\end{pmatrix} = \begin{pmatrix} 4 \pi\int_0^\infty \psi_1^* \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_1^* \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_1^* \psi_3 r^2dr\\
4 \pi\int_0^\infty \psi_2^* \psi_1 r^2dr & 4 \pi\int_0^\infty \psi_2^* \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_2^* \psi_3 r^2dr\\
4 \pi\int_0^\infty \psi_2^* \psi_3 r^2dr & 4 \pi\int_0^\infty \psi_3^* \psi_2 r^2dr & 4 \pi\int_0^\infty \psi_3^* \psi_3 r^2dr
\end{pmatrix}
$$

# Matriz S

S = sp.zeros(3)
S[0,0] = 4*sp.pi*sp.integrate(psi_1*psi_1*r**2,(r,0,sp.oo))
S[0,1] = 4*sp.pi*sp.integrate(psi_1*psi_2*r**2,(r,0,sp.oo))
S[0,2] = 4*sp.pi*sp.integrate(psi_1*psi_3*r**2,(r,0,sp.oo))
S[1,0] = 4*sp.pi*sp.integrate(psi_2*psi_1*r**2,(r,0,sp.oo))
S[1,1] = 4*sp.pi*sp.integrate(psi_2*psi_2*r**2,(r,0,sp.oo))
S[1,2] = 4*sp.pi*sp.integrate(psi_2*psi_3*r**2,(r,0,sp.oo))
S[2,0] = 4*sp.pi*sp.integrate(psi_3*psi_1*r**2,(r,0,sp.oo))
S[2,1] = 4*sp.pi*sp.integrate(psi_3*psi_2*r**2,(r,0,sp.oo))
S[2,2] = 4*sp.pi*sp.integrate(psi_3*psi_3*r**2,(r,0,sp.oo))
S

**Paso 3.** Resolver $\color{red}{\mathcal{H}\mathcal{C} = \mathcal{S}\mathcal{C} \mathcal{\epsilon}}$

Para ello utilizaremos la instrucción
~~~python
E,C = LA.eigh(H,S)
~~~
la cual resuelve directamente el problema $\mathcal{H}\mathcal{C} = \mathcal{S}\mathcal{C} \mathcal{\epsilon}$. La columna de $\mathcal{C}$ con la energía más baja nos indica los coeficientes de la combinación lineal de la función de onda.

# Resuelva HC=SCe

from scipy import linalg as LA

H = np.array(H).astype(np.float64)
S = np.array(S).astype(np.float64)

E,C = LA.eigh(H,S)
print(E)
print(C)

**Paso 4.** Sustituir los coeficientes en $\psi_{prueba} = \sum_{i=1} c_i \psi_i$ para construir la función de onda.

# Genere la función de onda

psi_p = C[0][0]*psi_1 + C[1][0]*psi_2 + C[2][0]*psi_3
psi_p

**Compruebe la normalización**

# Normalizacióm

4*sp.pi*sp.integrate(psi_p*psi_p*r**2,(r,0,sp.oo))

**Guarde sus funciones para graficar**

# Guarde funciones

psi_3g = psi_p
psi_3g_1 = psi_1
psi_3g_2 = psi_2
psi_3g_3 = psi_3

**Genere la gráfica con tres gaussianas**

# Gráfica

from matplotlib import pyplot as plt

s=1/sp.sqrt(sp.pi)*sp.exp(-sp.Abs(r))
lam_s = sp.lambdify(r,s,modules=['numpy'])
lam_psi_3g = sp.lambdify(r,psi_3g,modules=['numpy'])
lam_psi_3g_1 = sp.lambdify(r,psi_3g_1,modules=['numpy'])
lam_psi_3g_2 = sp.lambdify(r,psi_3g_2,modules=['numpy'])
lam_psi_3g_3 = sp.lambdify(r,psi_3g_3,modules=['numpy'])

r1 = np.linspace(-7,7,100)
psi_s = lam_s(r1)
psi_3g1 = lam_psi_3g(r1)
psi_3g_1 = lam_psi_3g_1(r1)
psi_3g_2 = lam_psi_3g_2(r1)
psi_3g_3 = lam_psi_3g_3(r1)

plt.plot(r1,psi_s,label="1s")
plt.plot(r1,psi_3g1,label="3g")
plt.plot(r1,psi_3g_1,label="psi_1",linestyle=':')
plt.plot(r1,psi_3g_2,label="psi_2",linestyle=':')
plt.plot(r1,psi_3g_3,label="psi_3",linestyle=':')
plt.legend()
plt.show()

## Átomo de hidrógeno dos gaussianas vs tres gaussianas

En la gráfica se compara la función de onda aproximada por el método variacional lineal con dos gausianas y con tres gausianas contra la solución exacta del átomo de hidrógeno del orbital 1s. 

from matplotlib import pyplot as plt

s=1/sp.sqrt(sp.pi)*sp.exp(-sp.Abs(r))
lam_s = sp.lambdify(r,s,modules=['numpy'])
lam_psi_2g = sp.lambdify(r,psi_2g,modules=['numpy'])
lam_psi_3g = sp.lambdify(r,psi_3g,modules=['numpy'])

r1 = np.linspace(-7,7,100)
psi_s = lam_s(r1)
psi_2g1 = lam_psi_2g(r1)
psi_3g1 = lam_psi_3g(r1)

plt.plot(r1,psi_s,label="1s")
plt.plot(r1,-psi_2g1,label="2g")
plt.plot(r1,psi_3g1,label="3g")
plt.legend()
plt.show()

## Referencias

- Atkins, P. W.; Friedman, R. Molecular Quantum Mechanics, 4th ed.; Oxford University Press: New York, 2005.
- Pilar, F. L. Elementary Quantum Chemistry; 2001.
- Zettili, N. Quantum Mechanics: Concepts and Applications, 2nd ed.; Wiley: Chichester, U.K, 2009.
- Levine, I. N. Quantum Chemistry, 5th ed.; Prentice Hall: Upper Saddle River, N.J, 2000.
- McQuarrie, D. A.; Simon, J. D. Physical Chemistry: A Molecular Approach; University Science Books: Sausalito, Calif, 1997.