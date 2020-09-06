# Software

Los `software de estructura electrónica` se dedican a calcular diversas propiedades de las moléculas utilizando la teoría que se ve en `química cuántica`. En general, esto nos sirve para `predecir la reactividad química`, desde poder predecir si una reacción procederá o no, hasta cosas más avanzadas como cinéticas de reacción y pKas o potenciales redox. En este notebook estaremos usando el software `psi4`, sin embargo, existe una gran cantidad de software, desde los libres y gratuitos hasta los privados y de pago.

```{note}
Algunos programas de estructura elecrónica son:
- Gaussian (https://gaussian.com/)
- Psi4 (http://www.psicode.org/)
- NWChem (http://www.nwchem-sw.org)
- QChem (http://www.q-chem.com/)
- TeraChem (http://www.petachem.com)
- deMon (http://www.demon-software.com)
- Orca (https://orcaforum.cec.mpg.de/)
- Molcas (http://www.molcas.org/)
- ADF (https://www.scm.com/)
- GAMESS (http://www.msg.chem.iastate.edu/gamess/)
- Quantum Espresso (https://www.quantum-espresso.org/)
```

## Aprendiendo a usar psi4

Para usar psi4, puede importarlo como si de una librería se tratase, es decir
~~~python
import psi4
~~~
**Importe psi4 en la siguiente celda**

# importe psi4

# importe psi4
import psi4

En la mayoría de los software es común (pero no obligatorio) que antes de mandar el cálculo de una molécula se asigne una cantidad de memoria RAM, por ejemplo 2 gb. 

En psi4 esto se hace mediante la instrucción
~~~python
psi4.set_memory("2 gb")
~~~

**Asigne memoria a su cálculo**

# Establezca memoria

# Establezca memoria
psi4.set_memory("2 gb")

```{note}
Los software usan la RAM asignada para guardar vectores y matrices como lo ha hecho en las prácticas anteriores. Si la memoria es suficiente, el programa guardará todo y el cálculo será más rápido, si no la hay, el cálculo será más lento. `A más memoria los cálculos tienden a ser igual o más rápidos`
```

```{warning}
La cantidad de memoria que puede asignar al cálculo depende de la cantidad de RAM que tenga su computadora. Recomendamos asignar menos memoria del total disponible ya que la memoria se reparte con los demás programas de su computadora. 
```

El siguiente paso es **declarar las coordenadas de los átomos que forman la molécula**. Para ello se pueden usar visualizadores como `Avogadro` o `IQmol`. También es posible obtener valores experimentales o calculados de https://cccbdb.nist.gov/ . En este caso utilizaremos los resultados experimentales de benceno.

Use las siguientes líneas para declarar la geometría
```{margin}
En este caso el 0 y el 1 indican la `carga` y `multiplicidad`, posteriormente viene el `X`, `Y`, `Z` y las unidades en las que se expresan las coordenadas.
```

```
benzene = psi4.geometry("""
0 1
C 0.0000 1.3970 0.0000
C 1.2098 0.6985 0.0000
C 1.2098 -0.6985 0.0000
C 0.0000 -1.3970 0.0000
C -1.2098 -0.6985 0.0000
C -1.2098 0.6985 0.0000
H 0.0000 2.4810 0.0000
H 2.1486 1.2405 0.0000
H 2.1486 -1.2405 0.0000
H 0.0000 -2.4810 0.0000
H -2.1486 -1.2405 0.0000
H -2.1486 1.2405 0.0000
units angstrom
""")
```

#Geometría

benzene = psi4.geometry("""
0 1
C 0.0000 1.3970 0.0000
C 1.2098 0.6985 0.0000
C 1.2098 -0.6985 0.0000
C 0.0000 -1.3970 0.0000
C -1.2098 -0.6985 0.0000
C -1.2098 0.6985 0.0000
H 0.0000 2.4810 0.0000
H 2.1486 1.2405 0.0000
H 2.1486 -1.2405 0.0000
H 0.0000 -2.4810 0.0000
H -2.1486 -1.2405 0.0000
H -2.1486 1.2405 0.0000
units angstrom
""")

Para realizar un cálculo de `energía` de una molécula `con la geometría` especificada arriba, es necesario especificar  un `método` y una `base` en la siguiente instrucción
~~~python
psi4.energy('método/base')
~~~

**Realice un cálculo con el método HF y la base 6-31G**

# Benzeno HF/6-31G

# Benzeno HF/6-31G
psi4.energy('HF/6-31G')

**Calcule la energía de benceno con MP2 y 6-31G**

# Benzeno MP2/6-31G

# Benzeno MP2/6-31G
psi4.energy('MP2/6-31G')

**Calcule la energía de benceno con el funcional M062X y la base 6-31G**

# Benzeno M062X/6-31G

# Benzeno M062X/6-31G
psi4.energy('M062X/6-31G')

Usualmente la geometría especificada no es necesariamente la geometría real. Es posible pedir al software que mueva los átomos hasta encontrar las coordenadas que representen un mínimo de energía con el método y base usados. Esto se llama `optimización de geometrías` y se hace con la línea

```
psi4.opt('método/base')
```

**Optimice la geometría de benceno con el método M062X y base 6-31G e imprima su energía**.

```{warning}
Este cálculo puede tardar entre 1 y 10 minutos dependiendo del procesador de cada computadora.
```

psi4.opt('M062X/6-31G')

## Referencias

- Smith, D. G. A.; Burns, L. A.; Sirianni, D. A.; Nascimento, D. R.; Kumar, A.; James, A. M.; Schriber, J. B.; Zhang, T.; Zhang, B.; Abbott, A. S.; et al. Psi4NumPy : An Interactive Quantum Chemistry Programming Environment for Reference Implementations and Rapid Development. Journal of Chemical Theory and Computation 2018, 14 (7), 3504–3511.
- https://github.com/psi4/psi4numpy/tree/master/Tutorials