# Software

Lista de algunos programas de química cuántica:
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

## Aprendiendo a usar psi4

Importamos psi4 con
~~~python
import psi4
~~~

# importe psi4
import psi4

Asignamos una cantidad de memoria, por ejemplo 2 gb, con
~~~python
psi4.set_memory("2 gb")
~~~

# Establezca memoria
psi4.set_memory("2 gb")

Declaramos las coordenadas de la molécula. Para ello usamos visualizadores como Avogadro. También podemos obtener valores experimentales de https://cccbdb.nist.gov/ . En este caso utilizaremos los resultados experimentales de benceno.

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

Realizamos el cálculo de energía con un **método** y una **base**. Por ejemplo, un cálculo de MP2 con base 6-31G sería
~~~python
psi4.energy('mp2/6-31G')
~~~

# Benzeno mp2/6-31G
psi4.energy('mp2/6-31G')

Calcule la energía de benceno con MP3 y 6-31G

# Benzeno mp3/6-31G

Calcule la energía de benceno con MP4 y 6-31G

# Benzeno mp4/6-31G

## Referencias

- Smith, D. G. A.; Burns, L. A.; Sirianni, D. A.; Nascimento, D. R.; Kumar, A.; James, A. M.; Schriber, J. B.; Zhang, T.; Zhang, B.; Abbott, A. S.; et al. Psi4NumPy : An Interactive Quantum Chemistry Programming Environment for Reference Implementations and Rapid Development. Journal of Chemical Theory and Computation 2018, 14 (7), 3504–3511.
- https://github.com/psi4/psi4numpy/tree/master/Tutorials

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.