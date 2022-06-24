#!/usr/bin/env python
# coding: utf-8

# # Software

# Los `software de estructura electrónica` se dedican a calcular diversas propiedades de las moléculas utilizando la teoría que se ve en `química cuántica`. En general, esto nos sirve para `predecir la reactividad química`, desde poder predecir si una reacción procederá o no, hasta cosas más avanzadas como cinéticas de reacción y pKas o potenciales redox. En este notebook estaremos usando el software `PySCF`, sin embargo, existe una gran cantidad de software, desde los libres y gratuitos hasta los privados y de pago.
# 
# ```{note}
# Algunos programas de estructura electrónica son:
# - Gaussian (https://gaussian.com)
# - PySCF (https://pyscf.org)
# - Psi4 (http://www.psicode.org)
# - NWChem (http://www.nwchem-sw.org)
# - QChem (http://www.q-chem.com)
# - TeraChem (http://www.petachem.com)
# - deMon2k (http://www.demon-software.com)
# - Orca (https://orcaforum.cec.mpg.de)
# - Molcas (http://www.molcas.org)
# - ADF (https://www.scm.com)
# - GAMESS (http://www.msg.chem.iastate.edu/gamess)
# - Quantum Espresso (https://www.quantum-espresso.org)
# - GPAW (https://wiki.fysik.dtu.dk/gpaw)
# ```

# ## Aprendiendo a usar PySCF

# ```{warning}
# Si está utilizando Google Colab o la ejecución en línea, debe de ejecutar al inicio el siguiente código
# ~~~
# !pip install pyscf
# !pip install geomeTRIC
# ~~~
# ```

# Para usar PySCF, puede importarlo como si de una librería se tratase, es decir
# 
# ~~~python
# import pyscf
# ~~~
# **Importe PySCF en la siguiente celda**

# In[1]:


# importe pyscf


# In[2]:


# Descomentar estas líneas si está en modo online

#!pip install pyscf
#!pip install geomeTRIC

import pyscf


# El siguiente paso es **declarar las coordenadas de los átomos que forman la molécula**. Para ello se pueden usar visualizadores como `Avogadro` o `IQmol`. También es posible obtener valores experimentales o calculados de https://cccbdb.nist.gov/ . En este caso utilizaremos los resultados experimentales de benceno.
# 
# Use las siguientes líneas para declarar la geometría y seleccionar una base
# ```{margin}
# Se asume que la molécula es neutra. Si desea editar la carga y/o el espín (note que esto es 2S, no la multiplicidad, 2S+1) puede usar instrucciones como
# 
# ~~~python
# mol.charge = 1
# mol.spin = 1
# ~~~
# 
# Recuerde volver a construir la molécula después de cualquier modificación.
# ```
# 
# ```
# benzene = pyscf.gto.Mole(atom = """
#     C  0.0000  1.3970  0.0000
#     C  1.2098  0.6985  0.0000
#     C  1.2098 -0.6985  0.0000
#     C  0.0000 -1.3970  0.0000
#     C -1.2098 -0.6985  0.0000
#     C -1.2098  0.6985  0.0000
#     H  0.0000  2.4810  0.0000
#     H  2.1486  1.2405  0.0000
#     H  2.1486 -1.2405  0.0000
#     H  0.0000 -2.4810  0.0000
#     H -2.1486 -1.2405  0.0000
#     H -2.1486  1.2405  0.0000
#     """,basis = "6-31G")
# benzene = benzene.build()
# ```
# 
# ```{margin}
# En este caso estamos llamando `benzene` a la molécula, usted puede usar cualquier otro nombre que prefiera.
# ```

# In[3]:


#Geometría


# In[4]:


benzene = pyscf.gto.Mole(atom = """
    C  0.0000  1.3970  0.0000
    C  1.2098  0.6985  0.0000
    C  1.2098 -0.6985  0.0000
    C  0.0000 -1.3970  0.0000
    C -1.2098 -0.6985  0.0000
    C -1.2098  0.6985  0.0000
    H  0.0000  2.4810  0.0000
    H  2.1486  1.2405  0.0000
    H  2.1486 -1.2405  0.0000
    H  0.0000 -2.4810  0.0000
    H -2.1486 -1.2405  0.0000
    H -2.1486  1.2405  0.0000
    """,basis = "6-31G")
benzene = benzene.build()


# En la mayoría de los software es común (pero no obligatorio) que antes de mandar el cálculo de una molécula se asigne una cantidad de memoria RAM, por ejemplo 2000 MB. 
# 
# En PySCF esto se hace mediante la instrucción
# ~~~python
# mol.max_memory = 2000
# ~~~
# 
# Cada vez que modifique algún parámetro de la molécula, debe volver a construirla antes de usarla.
# 
# **Asigne memoria a su cálculo**
# ```{margin}
# Note que `mol` es el nombre que le asignó a la molécula. En este ejemplo usamos el nombre de benzene.
# ```

# In[5]:


# Establezca memoria (y reconstruya la molécula)


# In[6]:


# Establezca memoria (en MB) (y reconstruya la molécula)
benzene.max_memory = 2000
benzene = benzene.build()


# ```{note}
# Los software usan la RAM asignada para guardar vectores y matrices como lo ha hecho en las prácticas anteriores. Si la memoria es suficiente, el programa guardará todo y el cálculo será más rápido, si no la hay, el cálculo será más lento. `A más memoria los cálculos tienden a ser igual o más rápidos`
# ```
# 
# ```{warning}
# La cantidad de memoria que puede asignar al cálculo depende de la cantidad de RAM que tenga su computadora. Recomendamos asignar menos memoria del total disponible ya que la memoria se reparte con los demás programas de su computadora. 
# ```

# Para realizar un cálculo de `energía` de una molécula `con la geometría y la base` especificada arriba, es necesario especificar un `método`. Por ejemplo, para hacer Hartree-Fock, escribiremos
# ~~~python
# from pyscf improt scf
# rhf = scf.RHF(mol)
# E = rhf.kernel()
# ~~~
# hemos de cambiar mol por el nombre que le asignamos a nuestra molécula.
# ```{margin}
# `SCF` significa Self-Consistent Field, `RHF` Significa Restricted Hartree-Fock. En moléculas de capa abierta se debe de usar `UHF` Unrestricted Hartree-Fock o `ROHF` Restricted-Open Hartree-Fock.
# ```
# 
# **Realice un cálculo con el método HF y la base 6-31G**
# ```{margin}
# En el caso de este ejemplo, la molécula se llama benzene, así que en vez de escribir mol, escribimos benzene.
# ```

# In[7]:


# Benceno HF


# In[8]:


# Benceno HF
from pyscf import scf
rhf = scf.RHF(benzene)
E = rhf.kernel()


# **Calcule la energía de benceno con MP2 y 6-31G**
# 
# ~~~python
# from pyscf import mp
# mp2 = mp.MP2(mol)
# Ecorr = mp2.kernel()
# ~~~
# ```{margin}
# `MP` significa Möller-Plesset y `MP2` indica que es de Segundo orden.
# ```
# ```{margin}
# Note que la última línea guarda en Ecorr solo la energía de correlación, no la energía de Hartree-Fock.
# ```

# In[9]:


# Benceno MP2/6-31G


# In[10]:


# Benceno MP2/6-31G
from pyscf import mp
mp2 = mp.MP2(benzene)
Ecorr = mp2.kernel()


# **Calcule la energía de benceno con el funcional M062X y la base 6-31G**
# 
# Use el código
# ~~~python
# from pyscf import dft
# rks = dft.RKS(mol)
# rks.xc = 'M062X'
# E = rks.kernel()
# ~~~
# ```{margin}
# `DFT` significa Density Functional Theory, y `RKS` significa Restricted Kohn-Sham. En moléculas de capa abierta se debe de usar `UKS` Unrestriced Kohn-Sham. 
# ```
# 
# ```{margin}
# La línea
# ~~~python
# rks.xc = 'M062X'
# ~~~
# define el funcional a usar.
# ```

# In[11]:


# Benceno M062X/6-31G


# In[12]:


# Benceno M062X/6-31G
from pyscf import dft
rks = dft.RKS(benzene)
rks.xc = 'M062X'
E = rks.kernel()


# Usualmente la geometría especificada no es necesariamente la geometría real. Es posible pedir al software que mueva los átomos hasta encontrar las coordenadas que representen un mínimo de energía con el método y base usados. Esto se llama `optimización de geometrías`, para lo cual hay que importar el optimizador y pasar la instancia del método. Por ejemplo, para optimizar con M062X, utilizaremos el rks que definimos antes
# ```
# from pyscf.geomopt.geometric_solver import optimize
# mol = geometric_opt(rks)
# ```
# 
# **Optimice la geometría de benceno con el método M062X y base 6-31G e imprima su energía**.
# 
# ```{warning}
# Este cálculo puede tardar entre segundos y unos minutos dependiendo del procesador de cada computadora.
# ```

# In[13]:


# Optimizacion de Geometría


# In[14]:


from pyscf.geomopt.geometric_solver import optimize
mol = optimize(rhf)


# ## Referencias

# - Q. Sun, T.C. Berkelbach, N.S. Blunt, G.H. Booth, S. Guo, Z. Li, J. Liu, J.D. McClain, E.R. Sayfutyarova, S. Sharma, S. Wouters, y G.K.-L. Chan, PySCF: the Python-based simulations of chemistry framework, Wiley Interdiscip. Rev. Comput. Mol. Sci. 8, e1340 (2018).
# - Q. Sun, X. Zhang, S. Banerjee, P. Bao, M. Barbry, N.S. Blunt, N.A. Bogdanov, G.H. Booth, J. Chen, Z.-H. Cui, J.J. Eriksen, Y. Gao, S. Guo, J. Hermann, M.R. Hermes, K. Koh, P. Koval, S. Lehtola, Z. Li, J. Liu, N. Mardirossian, J.D. McClain, M. Motta, B. Mussard, H.Q. Pham, A. Pulkin, W. Purwanto, P.J. Robinson, E. Ronca, E.R. Sayfutyarova, M. Scheurer, H.F. Schurkus, J.E.T. Smith, C. Sun, S.-N. Sun, S. Upadhyay, L.K. Wagner, X. Wang, A. White, J.D. Whitfield, M.J. Williamson, S. Wouters, J. Yang, J.M. Yu, T. Zhu, T.C. Berkelbach, S. Sharma, A.Y. Sokolov, y G.K.-L. Chan, Recent developments in the PySCF program package, J. Chem. Phys. 153, 024109 (2020).
# - Q. Sun, Libcint: An efficient general integral library for Gaussian basis functions, J. Comput. Chem. 36, 1664 (2015).
# - https://pyscf.org/quickstart.html
