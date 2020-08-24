# Ejercicio 9.

El $N_2O_4$ es un compuesto que se encuentra en la atmósfera. Con los cambios de temperatura se descompone en el radical $NO_2$ mediante la reacción:

$N_2O_4 <=> 2NO_2$

Las geometrías de $N_2O_4$ y $NO_2$ se dan a continuación:

Molécula: $N_2O_4$ Carga: 0 Multiplicidad: 1

|Átomo |x (Å)  |y (Å)  |z (Å) |
|------|-------|-------|------|
|N     |0.0000 |0.0000 |0.0000|
|N     |-1.7820|0.0000 |0.0000|
|O     |0.4516 |1.1010 |0.0000|
|O     |0.4516 |-1.1010|0.0000|
|O     |-2.2336|1.1010 |0.0000|
|O     |-2.2336|-1.1010|0.0000|

Molécula: $NO_2$ Carga: 0 Multiplicidad: 2

|Átomo |x (Å)  |y (Å)  |z (Å) |
|------|-------|-------|------|
|N     | 0.0000| 0.0000|0.0000|
|O     | 0.0000| 1.0989|0.4653|
|O     | 0.0000|-1.0989|0.4653|


**Pregunta 1.** Calule la energía de la molécula de $N_2O_4$ con HF y la base aug-cc-pvdz. 

import psi4
psi4.set_memory("2 gb")
psi4.geometry("""
0 1
   N       -4.84638        1.76109        0.00000
   N       -3.46888        1.78415        0.00000
   O       -2.82385        2.93169       -0.00000
   O       -2.85055        0.76276        0.00000
   O       -5.46471        2.78248        0.00000
   O       -5.49141        0.61355        0.00000
""")
n2o4,wfn=psi4.opt("HF/aug-cc-pvdz", return_wfn=True)
print(n2o4)

**Pregunta 2.** Calule la energía de la molécula de $NO_2$ con HF y la base aug-cc-pvdz. [Complete donde haga falta - Reemplace las X]

import psi4
NO2 = psi4.geometry("""
0 2
   N        XXXXXXX        XXXXXXX        XXXXXXX
   O        XXXXXXX        XXXXXXX        XXXXXXX
   O        XXXXXXX        XXXXXXX        XXXXXXX
units angstrom
""")
psi4.set_options({'reference': 'uhf'})
no2=psi4.optimize("XX/XXX-XX-XXXX")
print(no2)

**Pregunta a.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según HF.



**Pregunta 3.** Calule la energía de la molécula de $N_2O_4$ con DFT B3LYP y la base aug-cc-pvdz.



**Pregunta 4.** Calule la energía de la molécula de $NO_2$ con DFT B3LYP y la base aug-cc-pvdz.



**Preguna b.** Calcule el $\Delta U$ de la reacción $N_2O_4 <=> 2NO_2$ según DFT B3LYP.



# Referencias

- Parrish, R. M.; Burns, L. A.; Smith, D. G. A.; Simmonett, A. C.; DePrince, A. E.; Hohenstein, E. G.; Bozkaya, U.; Sokolov, A. Y.; Di Remigio, R.; Richard, R. M.; et al. **Psi4 1.1: An Open-Source Electronic Structure Program Emphasizing Automation, Advanced Libraries, and Interoperability.** Journal of Chemical Theory and Computation 2017, 13 (7), 3185–3197.


Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.