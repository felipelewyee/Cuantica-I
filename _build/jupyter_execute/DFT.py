# Teoría de funcionales de la densidad (DFT)

La teoría de los funcionales de la densidad se basa en que la densidad electrónica contiene toda la información del sistema, por lo que se puede extraer de ella la energía mediante un funcional

$$
E[\rho] = T_s[\rho] + U[\rho] + V_{nuc}[\rho] + E_{xc}[\rho]
$$

Donde
- $E[\rho]$ es la energía del sistema 
- $T_s[\rho]$ es la energía cinética de Kohn-Sham
- $U[\rho]$ es la energía de Hartree
- $V[\rho]$ es la interacción núcleo electrón
- $E_{xc} [\rho]$ es la energía de intercambio correlación

En la formulación de Kohn y Sham, la densidad electrónica ($\rho(r)$) se calcula a partir de los orbitales de Kohn-Sham ($\psi^{KS}_i(r)$)

$$
\rho(r) = \sum_i^{N} |\psi^{KS}_i(r)|^2
$$


```{warning}
Los teoremas de Hohenberg y Kohn, y Kohn y Sham prueban que existe un funcional universal que conecta a la densidad electrónica con la energía del sistema, pero este funcional no se conoce, particularmente por la parte de $E_{xc}[\rho]$. Se han realizan diversas propuestas de como construir este funcional, por lo que en la práctica existen cientos de funcionales de DFT.
```

```{note}
No existe el mejor funcional, sino que depende del sistema químico que se esté estudiando.
```

En general los funcionales aproximan de diferente manera el intercambio y la correlación, dependiendo de como se aproxime el funcional, estos se pueden clasificar en diferentes categorías. La clasificación fue propuesta por Perdew y se conoce como la `escalera de Jacob`.

```{warning}
Evaluar la contribución del funcional de intercambio-correlación requiere de realizar integración numérica. Esto se hace con un mallado entorno a la molécula. Existen diversos esquemas para colocar los puntos en el mallado, y para seleccionar cuantos puntos poner, estos puede cambiar de un software a otro, e incluso entre diferentes versiones de un mismo software.
```

**Importe psi4**

#Importe psi4

import psi4

Para ejemplificar el uso de estos funcionales, declare la molécula de agua.

```
h2o = psi4.geometry("""
0 1
O    0.0000    0.0000    0.1173
H    0.0000    0.7572   -0.4692
H    0.0000   -0.7572   -0.4692 
units angstrom
""")
```

# h2o

h2o = psi4.geometry("""
0 1
O    0.0000    0.0000    0.1173
H    0.0000    0.7572   -0.4692
H    0.0000   -0.7572   -0.4692 
units angstrom
""")

## Aproximación Local de la Densidad (LDA)

Fue una de las primeros aproximaciones y ya ha sido superada. En este caso, el intercambio se calcula mediante

$$
E_x^{LDA} = -\frac{3}{4} \left( \frac{3}{\pi} \right)^{1/3} \int \rho^{4/3} (r) dr 
$$

y la correlación se calcula mediante el funcional de Vosko, Wilk y Nusair
$$
E_c^{LDA} = \int \varepsilon_c^{VWN} dr
$$

Haga un cálculo de energía con SVWN y la base 6-311G con la siguiente instrucción
```
psi4.energy('SVWN/6-311G')
```

# SVWN

psi4.energy('SVWN/6-311G')

## Aproximación de Gradientes Generalizados (GGA)

Esto calcula la energía con base en la densidad electrónica y su gradiente

$$
E_{xc}^{GGA} = -\int \varepsilon_{xc}^{GGA} (\rho,\nabla \rho) dr
$$

Para ello separa los funcionales en intercambio y correlación.

$$
\varepsilon_{xc}^{GGA} = \varepsilon_{x}^{GGA} + \varepsilon_{c}^{GGA}
$$

Algunos funcionales de intercambio GGA son
- PWx86: Perdew-Wang 1986
- B88: Becke 1988
- PWx91: Perdew-Wang 1991
- PBE: Perdew-Burke-Ernzerhof
    
Algunos funcionales de correlación GGA son
- LYP: Lee-Yang-Parr
- Pc86: Perdew 1986
- PWc91: Perdew-Wang 1991
- PBE: Perdew-Burke-Ernzerhof
    
La combinación de estos funcionales genera los funcionales GGA. **Haga un cálculo de energía con PBE y la base 6-311G con la siguiente instrucción**
```
psi4.energy('PBE/6-311G')
```

psi4.energy('pbe/6-311G')

## Aproximación meta-GGA

Los funcionales meta-GGA usan la densidad electrónica, el gradiente de la densidad electrónica, y el laplaciano de la densidad electrónica.

$$
E_{xc}^{meta-GGA} = -\int \varepsilon_{xc}^{meta-GGA} (\rho,\nabla \rho,\nabla^2 \rho) dr
$$

Algunos ejemplos son:
- B95: Becke 1995
- TPSS: Tau-Perdew-Staroverov-Scuseria

**Haga un cálculo de energía con TPSS y la base 6-311G con la siguiente instrucción**
```
psi4.energy('TPSS/6-311G')
```

# TPSS

psi4.energy('TPSS/6-311G')

## Funcionales Híbridos

Mezclan un funcional de intercambio con el `intercambio de Hartree-Fock` en alguna proporción.

Algunos ejemplos de estos funcionales son:

- B3LYP
- PBE0
- M05-2X y M06-2X
- TPSSh

**Haga un cálculo de energía con M06-2X y la base 6-311G con la siguiente instrucción**
```
psi4.energy('M062X/6-311G')
```

#M062X

psi4.energy('M062X/6-311G')

## Referencias

- P. Hohenberg and W. Kohn, Phys. Rev. 136, B864 (1964).
- P. Atkins, Molecular Quantum Mechanics (Oxford University Press, 2005).
- D. Rappoport, N. R. M. Crawford, F. Furche, y K. Burke, Which functional should I choose?, (2008).
- W. Koch y M.C. Holthausen, A Chemist’s Guide to Density Functional Theory, (2001).
- K. Burke y L.O. Wagner, Int. J. Quantum Chem. 113, 96 (2013).
- K. Burke, J. Chem. Phys. 136, (2012).