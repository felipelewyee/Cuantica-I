# Ejercicio 10

## Catenaria: Un problema de cálculo de variaciones

Considere una cadena que cuelga de los puntos (-a,0) y (a,0), utilizaremos el cálculo de variaciones para obtener la función de mínima energía que describe a la cadena en un campo gravitacional constante.

### Energía de la cadena

Si la cadena tiene una longitud de $2L$ y una masa de $M$, la masa de un segmento $ds$ de la cadena está dada por
\begin{equation}
dm = \frac{M}{2L} ds
\end{equation}
y la energía de una porción de la cadena está dada por
\begin{equation}
dU = gydm
\end{equation}
Sustituyendo
\begin{equation}
dU = gy \frac{M}{2L} ds = gy \frac{M}{2L} \sqrt{(dx)^2 + (dy)^2} = gy \frac{M}{2L} \sqrt{1 + (y_x)^2} dx
\end{equation}
Por tanto, la energía está dada por
\begin{equation}
U = \int_{-a}^a gy \frac{M}{2L} \sqrt{1 + (y_x)^2} dx
\end{equation}

### Longitud de la cadena

El segmento de la cadena debe ser la longitud total
\begin{equation}
\int_{-a}^a ds = \int_{-a}^a \sqrt{1 + (y_x)^2} dx = 2L
\end{equation}
Esto es equivalente a
\begin{equation}
\int_{-a}^a \left( \sqrt{1 + (y_x)^2} - \frac{L}{a} \right) dx = 0
\end{equation}

### Multiplicador de Lagrange

A la función de energía se le puede agregar la restricción de la longitud
\begin{equation}
U = \int_{-a}^a gy \frac{M}{2L} \sqrt{1 + (y_x)^2} dx + \lambda \int_{-a}^a \left( \sqrt{1 + (y_x)^2} - \frac{L}{a} \right) dx
\end{equation}
donde $\lambda$ es un multiplicador indeterminado de Lagrange, nótese que formalente solo le hemos agregado un cero a la energía $U$.

Otra forma de escribir $U$ es
\begin{equation}
U = g \frac{M}{2L} \int_{-a}^a \left[ y \sqrt{1 + (y_x)^2} + \lambda' \left( \sqrt{1 + (y_x)^2} - \frac{L}{a} \right) \right ] dx = g \frac{M}{2L} \int_{-a}^a f dx
\end{equation}
con
\begin{equation}
f = y \sqrt{1 + (y_x)^2} + \lambda' \left( \sqrt{1 + (y_x)^2} - \frac{L}{a} \right)
\end{equation}


### Ecuación de Euler

La ecuación de Euler junto con la identidad de Beltrami permite obtener 
\begin{equation}
f - y_x \frac{\partial f}{\partial y_x} = C
\end{equation}

que al sustituir el valor de f resulta
\begin{equation}
y\sqrt{1+(y_x)^2} + \lambda'\left( \sqrt{1 + (y_x)^2} - \frac{L}{a} \right) - y_x \left( \frac{yy_x}{\sqrt{1+(y_x)^2}} + \lambda' \frac{y}{\sqrt{1+(y_x)^2}} \right) = C
\end{equation}

al despejar $y_x$
\begin{equation}
\frac{dy}{dx} = y_x = \sqrt{\left( \frac{a(y+\lambda')}{L\lambda' + aC} \right)^2 - 1} = \sqrt{\left( \frac{y+D}{k} \right)^2 - 1} 
\end{equation}

Al poner $x$ en función de $y$
\begin{equation}
dx = \frac{k}{\sqrt{\left( y+D \right)^2 - k^2 }} dy
\end{equation}

Finalmente, integramos
\begin{equation}
x = k arccosh \left(\frac{y+D}{k} \right) + F
\end{equation}

Despejamos $y$ en función de $x$
\begin{equation}
y = k cosh \left(\frac{x-F}{k} \right) - D
\end{equation}

### Condiciones de frontera

#### De donde se sostiene la cadena

Se debe cumplir $y(-a) = 0$
\begin{equation}
y(-a) = kcosh \left(\frac{-a-F}{k}\right)- D = kcosh \left(\frac{a+F}{k}\right)- D = 0
\end{equation}

Se debe cumplir $y(a) = 0$
\begin{equation}
y(a) = kcosh \left(\frac{a-F}{k}\right)- D = 0
\end{equation}

Por tanto
\begin{equation}
kcosh \left(\frac{a+F}{k}\right)- D = kcosh \left(\frac{a-F}{k}\right)- D
\end{equation}

Lo que solo se cumple si $F = 0$. Sustituyendo en $y(a)$
\begin{equation}
D = kcosh \left(\frac{a}{k}\right)
\end{equation}

#### De la longitud de la cadena

Ya podemos sustituir la forma de $y$ para ajustarla a la longitud de la cadena

\begin{equation}
2L = \int_{-a}^a \sqrt{1 + (y_x)^2} dx = \int_{-a}^a \sqrt{1 + sinh^2 \left(\frac{x}{k} \right)} dx = 2k sinh\left(\frac{a}{k}\right)
\end{equation}

Finalmente
\begin{equation}
k sinh\left(\frac{a}{k} \right) = L
\end{equation}

## Gráfica

La forma de la cadena está dada por
\begin{equation}
y = k cosh \left(\frac{x}{k} \right) - D
\end{equation}
con
\begin{equation}
D = kcosh \left(\frac{a}{k}\right)
\end{equation}
\begin{equation}
k sinh\left(\frac{a}{k} \right) = L
\end{equation}
donde la constante $D$ y $k$ se determinan únicamente por la longitud $L$ y los puntos $-a$ y $a$.

Damos un valor a las constantes $a$ y $L$. 

a = 10
L = 20

Importamos librerías.

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

Definimos la función:
\begin{equation}
f(k) = k sinh \left(\frac{a}{k} \right) - L
\end{equation}
El valor adecuado de $k$ es una raíz de $f(k)$, es decir
\begin{equation}
f(k) = 0
\end{equation}

def f(k):
    return k*np.sinh(a/k) - L

Buscamos la raíz de la función f(k).

k = brentq(f, 0.02, 100)
print(k)

Generamos la función

def D(k):
    return k*np.cosh(a/k)

def y(x):
    return k*np.cosh(x/k) - D(k)

Se realiza la gráfica

x = np.linspace(-a,a,100)
plt.scatter(-a,0)
plt.scatter(a,0)
plt.plot(x,y(x))
plt.ylim(-40,10)

## Referencias

- George, A. H. W. and F. H. Mathematical Methods for Physicist Seventh Edition; 2013, pp. 1081-1092. 
- Equilibrium Shape of a Rope Hanging from Two Endpoints http://semmat.dmf.unicatt.it/~paolini/divulgazione/mateott/catenaria/catenary/catenary.htm?fbclid=IwAR2VWhmP33C48hPv0GgizjbyP4t94EAKuStE_gJjoKKGtYNmE6YMGhlU8ig (accessed Feb 10, 2020).    

## Autores

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.