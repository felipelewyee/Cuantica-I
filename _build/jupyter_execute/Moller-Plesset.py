# Moller-Plesset (MPn)

Las ecuaciones de Moller-Plesset surgen de aplicar la teoría de perturbaciones a Hartree-Fock. Para ello se toma como sistema conocido el Hamiltoniano de Hartree-Fock ($\mathcal{H}_0$) y se le aplica  una perturbación ($\mathcal{V}$) para convertirlo en el Hamiltoniano del sistema con correlación electrónica ($\mathcal{H}$).

$$
\mathcal{H} = \mathcal{H}_0 + \mathcal{V}
$$

Para comenzar, **importe las librerías numpy y psi4**.

# Importe librerías

import numpy as np
import psi4

**Defina una molécula de su interés.** A continuación se proporciona la geometría de la molécula de hidrógeno, aunque puede usar cualquier otro sistema. Se recomienda que sea pequeño por tiempo de ejecución.
```
H 0.0000  0.0000 0.0000
H 0.0000  0.0000 0.7414 
```

# Defina molécula

mol = psi4.geometry("""
0 1
H 0.0000  0.0000 0.0000
H 0.0000  0.0000 0.7414 
units angstrom
symmetry C1
""")

Ya que la teoría de Moller-Plesset corrige el Hamiltoniano de Hartree-Fock, **realice un cálculo de Hartree-Fock y recupere la energía y la función de onda**. La energía resultante contiene tanto energía electrónica como nuclear. Seleccione la base que desee, por ejemplo, 6-31G

# Calcule E_HF y wfn

E_HF, wfn = psi4.energy('HF/6-31G',return_wfn=True)
print("E_HF =",E_HF)

**Obtenga la energía nuclear**

# E_nuc

E_nuc = mol.nuclear_repulsion_energy() 
print("E_nuc =",E_nuc)

**Obtenga los coeficientes de los orbitales moleculares a partir de la función de onda**

# Coeficientes

C = np.asarray(wfn.Ca())

**Obtenga el número de funciones de base ($N_{bf}$), de orbitales moleculares ($N_{mo})$ y de electrones ($N_{e}$)**

# nbf, nmo y ne

nbf = len(C)
nmo = wfn.nmo()
ne = wfn.nalpha() + wfn.nbeta()
print("nbf =",nbf," nmo =",nmo," ne =",ne)

Vamos a necesitar integrales de repulsión electónica

$$
[\mu \nu | \sigma \lambda] = \int \int \frac{\mu(r_1) \nu(r_1) \sigma(r_2) \lambda(r_2)}{ |r_1 - r_2| } dr_1 dr_2
$$
donde $\mu$, $\nu$, $\sigma$ y $\lambda$ se refieren a orbitales atómicos.

**Obtenga las las integrales $[\mu\nu|\sigma\lambda]$ y guardarlas en la variable I_AO**
`````{tip}
Puede usar las siguientes líneas
~~~
mints = psi4.core.MintsHelper(wfn.basisset())
I_AO = np.asarray(mints.ao_eri())
~~~
`````

# ERIs I_AO

mints = psi4.core.MintsHelper(wfn.basisset())
I_AO = np.asarray(mints.ao_eri())

```{Note}
Los software de estructura electrónica usan formas optimizados de las ecuaciones aquí mostradas. Este notebook ha sido creado sólo con fines didácticos.
```

La teoría de Moller-Plesset requiere que estas integrales sean transformadas a orbital molecular. Recuerde que un orbital molecular es una combinación lineal de orbitales atómicos

$$
p(r) = \sum_\mu C_{\mu p} \mu (r)
$$

En el siguiente recuadro **declare una variable I_MO de dimensión ($N_{MO}$,$N_{MO}$,$N_{MO}$,$N_{MO}$), y lleve a cabo la transformación**

$$
[pq|rt] = \sum_{\mu\nu\sigma\lambda}^N C_{\mu p} C_{\nu q} C_{\sigma r} C_{\lambda s} [\mu \nu | \sigma \lambda]
$$

```{tip}
Utilice 4 for para recorrer los orbitales moleculares, y 4 for para recorrer las funciones de base.
```

`````{margin}
Alternativamente puede usar la instrucción
~~~
I_MO = np.einsum('mp,nq,sr,lt,mnsl->pqrt',C,C,C,C,I_AO,optimize=True)
~~~
La cual es sustancialmente más rápida que los 8 for. Esto puede hacer su programa más rápido si su sistema no es tan pequeño.
`````

# ERIs I_MO

I_MO = np.zeros((nmo,nmo,nmo,nmo))
for p in range(nmo):
    for q in range(nmo):
        for r in range(nmo):
            for t in range(nmo):
                for m in range(nbf):
                    for n in range(nbf):
                        for s in range(nbf):
                            for l in range(nbf):
                                I_MO[p][q][r][t] = I_MO[p][q][r][t] + C[m][p]*C[n][q]*C[s][r]*C[l][t]*I_AO[m][n][s][l]

**Obtenga las energías de los orbitales moleculares ($\varepsilon_a$)**
`````{tip}
Puede usar el siguiente código
~~~
epsilon = wfn.epsilon_a()
~~~
`````

# epsilon

epsilon = np.array(wfn.epsilon_a())

Recordando teoría de perturbaciones, es posible realizar correcciones de n-orden a la energía $E_i^{(n)}$, tal que la energía electrónica total del estado basal ($i=0$) es

$$
E = E_0^{(0)} + E_0^{(1)} + E_0^{(2)} + E_0^{(3)} + E_0^{(4)} ...
$$

dependiendo del n-orden hasta el que se haga la corrección sobre la energía al cálculo se le denomina MPn, siendo los más comunes son MP2, MP3 y MP4.

En Moller-Plesset estos términos son:

$$
E_0^{(0)} = 2\sum_{a}^{N_e/2} {\varepsilon_a}
$$

$$
E_0^{(1)} = -2 \sum_{a}^{N_e/2}\sum_{b}^{N_e/2} [aa|bb] - [ab|ba]
$$

donde los índices $a$, $b$ refieren a orbitales moleculares ocupados, y los términos $\varepsilon_a$, $\varepsilon_b$, a sus energías.

**Calcule $E_0^{(0)}$ y $E_0^{(1)}$**.

```{Caution}
Aunque las sumas empiezan en el primer orbital molecular ocupado, recuerde que en Python los índices empiezan en cero.
```

# E_0 y E_1

E_0 = 0
for a in range(int(ne/2)):
    E_0 = E_0 + 2*epsilon[a]
    
E_1 = 0
for a in range(int(ne/2)):
    for b in range(int(ne/2)):
        E_1 = E_1 - 2 * I_MO[a][a][b][b] + I_MO[a][b][b][a]

**Calcule la energía total de MP1**, recuerde sumar la energía nuclear, es decir

$$
E_{MP1} = E_{nuc} + E_0^{(0)} + E_0^{(1)}
$$

# E_MP1

E_MP1 = E_nuc + E_0 + E_1
print("E_MP1 =",E_MP1)

```{admonition} Pregunta
:class: warning

Calcule la diferencia entre $E_{MP1}$ y la energía de Hartree-Fock que calculó al inicio, **¿Cuál es la energía de correlación? ¿Cómo se relaciona $E_{MP1}$ con $E_{HF}$?
```

La energía de Hartree-Fock se calcula como

$$
E_{HF} = E_{nuc} + 2\sum_a^{N_e/2} {\varepsilon_a} -2 \sum_{a}^{N/2}\sum_{b}^{N/2} [aa|bb] - [ab|ba]
$$

esta es exactamente la misma expresion que $E_{MP1}$. La primera corrección a la energía aparece en $E_{MP2}$. La corrección a segundo orden es

$$
E_0^{(2)} = 2 \sum_{a}^{N_e/2}\sum_{b}^{N_e/2}\sum_{r=N_e+1}^{N_{bf}}\sum_{s=N_e+1}^{N_{bf}} \frac{[ar|bs][ra|sb]}{\varepsilon_{a}+\varepsilon_{b}-\varepsilon_{r}-\varepsilon_{s}} - \sum_{abrs}^{N_e/2} \frac{[ar|bs][rb|sa]}{\varepsilon_{a}+\varepsilon_{b}-\varepsilon_{r}-\varepsilon_{s}}
$$

donde $r$, $s$ son los orbitales moleculares desocupados.

**Calcule $E_0^{(2)}$, la energía de $MP2$** dada por

$$
E_{MP2} = E_{nuc} + E_0^{(0)} + E_0^{(1)} + E_0^{(2)}
$$

**y la energía de correlación**

$$
E_{corr} = E_{MP2} - E_{HF}
$$

# E_2, E_MP2 y E_corr

E_2 = 0
for a in range(int(ne/2)):
    for b in range(int(ne/2)):
        for r in range(int(ne/2),nbf):
            for s in range(int(ne/2),nbf):
                E_2 = E_2 + 2*(I_MO[a][r][b][s]*I_MO[r][a][s][b])/(epsilon[a] + epsilon[b] - epsilon[r]- epsilon[s])
                E_2 = E_2 -   (I_MO[a][r][b][s]*I_MO[r][b][s][a])/(epsilon[a] + epsilon[b] - epsilon[r]- epsilon[s])
                
E_MP2 = E_nuc + E_0 + E_1 + E_2
print("E_MP2 =",E_MP2)

print("E_corr =",E_MP2-E_HF)

```{important}
Generalizando, la energía total de MPn es

$$
E_{MPn} = E_{nuc} + E_0^{(0)} + E_0^{(1)} + E_0^{(2)} + ... + E_0^{(n)}
$$

Al restarle la energía de Hartree-Fock se obtiene

$$
E_{corr} = E_{MPn} - E_{HF} = E_0^{(2)} + ... + E_0^{(n)}
$$
```

**Adapte la instrucción a MP2 junto con la base que usó para comprobar sus resultados**.
```
E_MPn = psi4.energy('MPn/Base')
print("E_MPn =",E_MPn)
```

# E_MP2 psi4

E_MP2 = psi4.energy('MP2/6-31G')
print("E_MP2 =",E_MP2)