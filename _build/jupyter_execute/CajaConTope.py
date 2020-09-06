# Caja de potencial con un tope en el medio

Este problema es continuación de la partícula en la caja. Para ello, planteemos una caja de $-L$ a $L$, con el potencial definido por

$$
V(x) = \left\{
  \begin{array}{lll}
  \infty      & \mathrm{si\ } x < -L & \\
  0      & \mathrm{si\ } -L \le x < -a & II\\
  V & \mathrm{si\ } -a \le x \le a & III \\
  0      & \mathrm{si\ } a < x \le L & IV\\
  \infty     & \mathrm{si\ } x > L & 
  \end{array}
  \right.
$$

Es decir, el potencial vale infinito fuera de la caja, cero en las zonas de la izquierda y la derecha (de $-L$ a $-a$ y de $+a$ a $+L$) y vale $V$ en el centro (de $-a$ a $+a$).

```{admonition} Para pensar
:class: tip
Considere que tiene una partícula moviéndose dentro de la caja con una energía menor que V. De manera clásica, la partícula no podría pasar de un lado de la caja al otro porque no tiene suficiente energía para atravesar el potencial. ¿Qué pasará cuánticamente?
```

La función de onda se obtiene resolviendo la ecuación de Schrödinger

$$
\left( -\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V(x) \right) \psi(x) = E \psi(x)
$$

Como se vio antes, la función de onda vale cero afuera de la caja. Por lo que se puede plantear la ecuación de Schrödinger por secciones.

```{admonition} Inserto matemático: Hamiltoniano por secciones
:class: dropdown

Si analizamos la ecuación de Schrodiger por zonas se tiene:

| Zona      | Hamiltoniano | Función de onda | Constantes |
|:----------------:|:---------:|:--------:|:--------:|
| I | $-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi_I(x) = E \psi_I(x)$ | $\psi_I(x) = A sin(k_1 x) + Bcos(k_1x)$ | $k_1^2 = \frac{2mE}{\hbar^2}$ |
| II | $-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi_{II}(x) + U\psi_{II}(x)= E \psi_{II}(x)$ | $\psi_{II}(x) = C e^{k_2 x} + De^{-k_2x}$ | $k_2^2 = -\frac{2m(E-U)}{\hbar^2} = \frac{2m(U-E)}{\hbar^2}$ |
| III | $-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi_{III}(x) = E \psi_{III}(x)$ | $\psi_{III}(x) = A' sin(k_1 x) + B' cos(k_1x)$ | $k_1^2 = \frac{2mE}{\hbar^2}$ |

```

**Importe las siguientes librerías**
- numpy
- pyplot de matplotlib
- optimize de scipy
- integrate de scipy

# Importe librerías

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
from scipy import integrate

**Establezca valores para las constantes** $\hbar$, $m$, $V$, $a$, $L$.

# De valor a las constantes

hbar = 1
m = 1.0
V = 50.0
a = 0.2
L = 1.2

La función de onda debe ser contínua, por lo que podemos igualar la función de onda en el punto donde se unen las zonas, y obtener nuevas ecuaciones.

```{admonition} Inserto matemático: Condiciones de Frontera
:class: dropdown
| Zonas | Condición | Ecuación |
|:---: |:---: | :---:    |
| Inicio y I | $\psi_I(-L) = 0$ | $B = A tan(k_1 L)$ |
| I y II | $\psi_{I}(-a) = \psi_{II}(-a)$ | $-A sin(k_1 a) + Bcos(k_1 a) = C e^{-k_2a} + D e^{k_2a}$ |
|I y II | $\psi'_{I}(-a) = \psi'_{II}(-a)$ | $k_1(A cos(k_1 a) + Bsin(k_1 a)) = k_2 (C e^{-k_2a} - D e^{k_2a})$|
| II y III | $\psi_{II}(a) = \psi_{III}(a)$ | $C e^{-k_2a} + D e^{k_2a} = A' sin(k_1 a) + B'cos(k_1 a)$|
| III y Final | $\psi_{III}(L) = 0$ | $B' = -A' tan(k_1 L)$|
```

A partir de aquí podemos ayudarnos de la simetría del problema.

## Simetría Par

Empezaremos asumiendo que lo que esta del lado izquierdo del potencial es simétrico respecto a lo que esta del lado derecho, es decir el problema tiene simetría par ($\psi_{II}(x) = \psi_{II}(-x)$ y $C = D$). Al considerar esta condición en las ecuaciones de continuidad, se obtiene las siguiente ecuación

$$
tanh \left(\sqrt{\frac{2m(V-E)}{\hbar^2}} a \right) tan \left(\sqrt{\frac{2mE}{\hbar^2}}(a-L) \right) = \sqrt{\frac{E}{V-E}}
$$

En esta ecuación no es trivial despejar E, sin embargo, la igualdad sólo se cumplirá con la E correcta. Una forma más simple es elevar al cuadrado y pasar todo a la derecha, tal que definamos $f(E)$

$$
f(E) = tanh^2 \left(\sqrt{\frac{2m(V-E)}{\hbar^2}} a \right) tan^2 \left(\sqrt{\frac{2mE}{\hbar^2}}(a-L) \right) - \frac{E}{V-E}
$$

Cuando se tenga el $E$ correcto se cumplirá $f(E) = 0$, así que sólo tenemos que buscar los ceros (o raíces) de la función.

**Defina la función $f(E)$**

# f(E)

def f(E): 
    
    arg1 = np.sqrt(2*m*(V-E)/hbar**2)*a
    arg2 = np.sqrt(2*m*E/hbar**2)*(a-L)
    
    return np.tanh(arg1)**2*np.tan(arg2)**2 - E/(V-E)

Para encontrar los valores de E que hacen que f(E) se vuelva cero, cree un conjunto de puntos de E con muchos puntos entre 0 y V, puede usar la instrucción
```
E_dominio = np.linspace(0,V,10000)
```

# E_dominio

E_dominio = np.linspace(0,V,10000)

**Para cada uno de estos puntos evalúe si f(E) es menos a $10^{-2}$**, si la condición se cumple el valor de E es un buen candidato para ser una raíz de f(E). Haga una lista con los valores de E que cumplieron el criterio, este será su `primer guess`.

E_primerguess = []
for E_i in E_dominio:
    if (abs(f(E_i))<1e-2):
        E_primerguess.append(E_i)

Python tiene funciones especiales para buscar raíces partiendo de cierto punto. La siguiente línea busca la raíz de f(E) más cercana a un punto E_i
```
E = newton(f,x0=E_i)
```

Para cada valor de energía de su primer guess, utilice el método de Newton para encontrar la raíz más cercana y guárdela en una lista si la diferencia con la última raíz es mayor a 0.1. Este será su `segundo guess`.

E_segundoguess = [0]
for E_i in E_primerguess:
    E = optimize.newton(f,x0=E_i)
    if (abs(E_segundoguess[-1] - E) > 0.1 ):
        E_segundoguess.append(E)

Imprima su segundo guess

# Impresión

print(E_segundoguess)

Defina funciones para

$$
k_1 = \frac{2mE}{\hbar^2}
$$

$$
k_2 = \frac{2m(V-E)}{\hbar^2}
$$

# Defina funciones

def k1(E): return np.sqrt(2*m*E/hbar**2)
def k2(E): return np.sqrt(2*m*(V-E)/hbar**2)

Defina funciones para

$$
\psi_{I}(x) = \frac{e^{-k_2Ea}+e^{k_2Ea}}{-sin(k_1Ea)+tan(k_1EL)cos(k_1Ea)} sin(k_1Ex)+tan(k_1EL)cos(k_1Ex)
$$

$$
\psi_{II}(x) = e^{k_2Ex}+e^{-k_2Ex}
$$

$$
\psi_{III}(x) = \frac{e^{k_2Ea}+e^{-k_2Ea}}{sin(k_1Ea)-tan(k_1EL)cos(k_1Ea)}(sin(k_1Ex)-tan(k_1EL)cos(k_1Ex))
$$

# Defina funciones

def psi_I(x): return (np.exp(-k2(E)*a)+np.exp(k2(E)*a))/(-np.sin(k1(E)*a)+np.tan(k1(E)*L)*np.cos(k1(E)*a))*(np.sin(k1(E)*x)+np.tan(k1(E)*L)*np.cos(k1(E)*x))
def psi_II(x): return np.exp(k2(E)*x)+np.exp(-k2(E)*x)
def psi_III(x): return (np.exp(k2(E)*a)+np.exp(-k2(E)*a))/(np.sin(k1(E)*a)-np.tan(k1(E)*L)*np.cos(k1(E)*a))*(np.sin(k1(E)*x)-np.tan(k1(E)*L)*np.cos(k1(E)*x))

Cree tres dominios de 1000 puntos para la función de onda, tal que

$$
x_1 \in [-L,-a]
$$

$$
x_2 \in [-a,a]
$$

$$
x_3 \in [a,L]
$$

# x_1, x_2 y x_3

x1 = np.linspace(-L,-a,10000)
x2 = np.linspace(-a,a,10000)
x3 = np.linspace(a,L,10000)

Utilice las energías del segundo guess para graficar las funciones de onda.

for E in E_segundoguess:
    if(E>0):
        norm = 0.0
        
        norm = norm + integrate.quad(lambda x: psi_I(x)*psi_I(x), -L, -a)[0]
        norm = norm + integrate.quad(lambda x: psi_II(x)*psi_II(x), -a, a)[0]
        norm = norm + integrate.quad(lambda x: psi_III(x)*psi_III(x), a, L)[0]
                
        plt.plot(x1,psi_I(x1)/np.sqrt(norm))
        plt.plot(x2,psi_II(x2)/np.sqrt(norm))
        plt.plot(x3,psi_III(x3)/np.sqrt(norm)) 
        
        prob = integrate.quad(lambda x: psi_II(x)*psi_II(x), -a, a)[0]/norm
        print("E: " + str(E) + " Probabilidad de [-a,a]: " + str(prob))
        
        plt.plot(x1,psi_I(x1)*psi_I(x1)/norm)
        plt.plot(x2,psi_II(x2)*psi_II(x2)/norm)
        plt.plot(x3,psi_III(x3)*psi_III(x3)/norm)

        plt.show()
    else:
        print("Zero")

## Simetría impar

Ahora asumirémos que lo que esta del lado izquierdo del potencial es antisimétrico respecto a lo que esta del lado derecho, es decir el problema tiene simetría par ($\psi_{II}(x) = -\psi_{II}(-x)$ y $C = -D$). Al considerar esta condición en las ecuaciones de continuidad, se obtienen las siguiente ecuación

$$
tanh^{-1} \left(\sqrt{\frac{2m(V-E)}{\hbar^2}} a \right) tan \left(\sqrt{\frac{2mE}{\hbar^2}}(a-L) \right) = \sqrt{\frac{E}{V-E}}
$$

En esta ecuación no es trivial despear E, pero la igualdad sólo se cumplirá con la E correcta. Una forma más simple es elevar al cuadrado y pasar todo a la derecha, tal que definamos $f(E)$

$$
f(E) = \left( tanh^{-1} \left(\sqrt{\frac{2m(V-E)}{\hbar^2}} a \right) \right)^2 tan^2 \left(\sqrt{\frac{2mE}{\hbar^2}}(a-L) \right) - \frac{E}{V-E}
$$

Cuando se tenga el $E$ correcto se cumplirá $f(E) = 0$, así que sólo tenemos que buscar los ceros (o raíces) de la función.

**Defina la función $f(E)$**

def f(E): 
    arg1 = np.sqrt(2*m*(V-E)/hbar**2)*a
    arg2 = np.sqrt(2*m*(E)/hbar**2)*(a-L)
    
    return (1/np.tanh(arg1))**2*np.tan(arg2)**2 - E/(V-E)

Obtenga su guess de valores de energía

```{tip}
- Genere un conunto de 1000 puntos de E de 0 a V
- Seleccione aquellos para los que f(E) es menor que $10^{-2}$
- Utilice el método de Newton para obtener valores únicos de energía
```

E_dominio = np.linspace(0,V,10000)

E_primerguess = []
for E in E_dominio:
    if (abs(f(E))<1e-2):
        E_primerguess.append(E)
    
E_segundoguess = [0]
for E in E_primerguess:
    E_i = optimize.newton(f,x0=E)
    if (abs(E_segundoguess[-1] - E_i) > 0.1 ):
        E_segundoguess.append(E_i)

Defina funciones para:

$$
\psi_{I}(x) = \frac{e^{-k_2E*a}-e^{k_2Ea}}{-sin(k_1Ea)+tan(k_1EL)cos(k_1Ea)}(sin(k_1Ex)+tan(k_1EL)np.cos(k_1Ex))
$$

$$
\psi_{II}(x) = e^{k_2Ex}-e^{-k_2Ex}
$$

$$
\psi_{III}(x) = \frac{e^{k_2Ea}-e^{-k_2Ea}}{sin(k_1Ea)-tan(k_1EL)cos(k_1Ea)}(sin(k_1Ex)-tan(k_1EL)cos(k_1Ex))
$$

# Funciones

def psi_I(x): return (np.exp(-k2(E)*a)-np.exp(k2(E)*a))/(-np.sin(k1(E)*a)+np.tan(k1(E)*L)*np.cos(k1(E)*a))*(np.sin(k1(E)*x)+np.tan(k1(E)*L)*np.cos(k1(E)*x))
def psi_II(x): return np.exp(k2(E)*x)-np.exp(-k2(E)*x)
def psi_III(x): return (np.exp(k2(E)*a)-np.exp(-k2(E)*a))/(np.sin(k1(E)*a)-np.tan(k1(E)*L)*np.cos(k1(E)*a))*(np.sin(k1(E)*x)-np.tan(k1(E)*L)*np.cos(k1(E)*x))

Realice las gráficas de la función de onda

# Graficas

for E in E_segundoguess:
    if(E>0):
        norm = 0.0
        
        norm = norm + integrate.quad(lambda x: psi_I(x)*psi_I(x), -L, -a)[0]
        norm = norm + integrate.quad(lambda x: psi_II(x)*psi_II(x), -a, a)[0]
        norm = norm + integrate.quad(lambda x: psi_III(x)*psi_III(x), a, L)[0]
                
        plt.plot(x1,psi_I(x1)/np.sqrt(norm))
        plt.plot(x2,psi_II(x2)/np.sqrt(norm))
        plt.plot(x3,psi_III(x3)/np.sqrt(norm)) 
        
        prob = integrate.quad(lambda x: psi_II(x)*psi_II(x), -a, a)[0]/norm
        print("E: " + str(E) + " Probabilidad de [-a,a]: " + str(prob))
        
        plt.plot(x1,psi_I(x1)*psi_I(x1)/norm)
        plt.plot(x2,psi_II(x2)*psi_II(x2)/norm)
        plt.plot(x3,psi_III(x3)*psi_III(x3)/norm)

        plt.show()
    else:
        print("Zero")

Con base en lo anteriormente visto, responda a la siguiente frase con verdadero o falso.

**Si la partícula tiene energía menor que V, es imposible encontrarla en el intevalo [-a,a] (Cierto/Falso)**

**Respuesta**

## Autoría

Esta es una contribución del Dr. Carlos Amador Bedolla.