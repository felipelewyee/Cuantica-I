#!/usr/bin/env python
# coding: utf-8

# # La Química Cuántica a tu alcance

# El presente libro digital tiene como objetivo ser un manual de prácticas para la asignatura de [`Química Cuántica I (clave 1404)`](https://quimica.unam.mx/wp-content/uploads/2017/03/1404qcuanticauno.pdf) de la Facultad de Química (UNAM). El libro digital es un esfuerzo compartido de Juan Felipe Huan Lew Yee, Eduardo Barrios Vargas y Jorge Martín del Campo Ramírez, en conjunto con los profesores del Departamento de Física y Química Teórica. Toda retroalimentación es bienvenida, y se le invita a compartirla a los correos:
# 
# M. en C. Juan Felipe Huan Lew Yee: [felipe.lew.yee@gmail.com](mailto:felipe.lew.yee@gmail.com)
# 
# Dr. Jorge Martín del Campo Ramirez: [jormacara@gmail.com](mailto:jormacara@gmail.com)
# 
# Dr. José Eduardo Barrios Vargas: [j.e.barrios@gmail.com](mailto:j.e.barrios@gmail.com)
# 
# ## Uso del manual
# 
# Se recomienda la lectura en el orden de los capítulos, el cual puede seguirse utilizando la barra de navegación que aparece del lado izquierdo en computadoras, o con el menú emergente en dispositivos móviles.
# 
# Se ha buscado reducir lo más posible los antecedentes matemáticos requeridos, por lo que en la mayoría de los ejercicios se plantea el problema y se pasa directamente a la solución sin realizar el desarrollo. Se asume que las ecuaciones han sido previamente revisadas en clase, sin embargo, en muchos notebooks pueden encontrarse partes clave del desarrollo matemático ocultas en una celda como la siguiente
# ```{admonition} Inserto matemático: Celda con contenido matemático oculto
# :class: dropdown
# 
# La lectura de este apartado matemático oculto es opcional si ya se conoce la deducción previa de las ecuaciones, pero aquí está por si alguien lo necesita. Aquí van algunas funciones trigonométricas
# 
# $$
# y = sin(x)
# $$
# 
# $$
# y = cos(x)
# $$
# 
# $$
# y = tan(x) = \frac{sin(x)}{cos(x)}
# $$
# ```

# A lo largo de los ejercicios encontrará celdas con el código correspondiente, y celdas vacías con una nota para indicar que debe ser llenada siguiendo alguna instrucción, **como la siguiente**

# In[1]:


# Llene la siguiente celda con su código


# La primera vez que se presente una instrucción se mostrará directamente para copiar y pegar en la celda, por ejemplo, el siguiente código muestra como realizar una impresión
# ```
# print("Bienvenidos al manual de Química Cuántica I")
# ```

# Posteriormente, cuando se solicite realizar algo usando código ya presentado, se mostrará la respuesta oculta para permitir que se intente resolver el ejercicio en la celda vacía por cuenta propia. Para revelar la respuesta hay que presionar en el botón con el símbolo + como en la siguiente celda

# In[2]:


print("Bienvenidos al manual de Química Cuántica I")


# ## Instalación y Ejecución
# 
# Este libro permite interactuar con el contenido. Seleccione la opción que le resulte más conveniente acorde a su sistema operativo.
# 
# ### En Computadora
# 
# Puede descargar los notebooks en formato *.ipynb con el botón con la flecha apuntando hacia abajo y dando clic en .ipynb. 
# 
# ![Descarga de ipynb](images/download_ipynb.png)
# 
# Para ejecutarlos localmente es su computadora se requiere una instalación previa de Python. En caso de no tenerla puede instalar [Anaconda](https://www.anaconda.com/products/individual). Una vez descargado puede interactuar con el notebook utilizando Jupyter Notebook.
# 
# ````{admonition} Windows
# :class: tip
# :class: dropdown
# 
# Abra la Powershell que encontrará en el menú de inicio y escriba
# ```
# Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
# ```
# Esto habilitará el `Windows Subsystem for Linux`. Reinicie su computadora una vez completada la instrucción. Nota: Puede que necesite ejecutar Powershell en modo administrador.
# 
# ![Powershell](images/powershell.png)
# 
# Abra la `Microsoft Store` y busque `Ubuntu`. De clic en obtener/instalar y espere a que termine la instalación. 
# ![ubuntu](images/ubuntu.png)
# 
# De clic en iniciar. Tardará unos minutos y le pedirá que estableza un usuario y contraseña.
# ![iniciar](images/iniciar.png)
# 
# Copie y pegue la siguiente línea para instalar `wget`. Le pedirá su contraseña.
# ```
# sudo apt install wget
# ```
# ![install-wget](images/install-wget.png)
# 
# Copie y pegue la siguiente línea para descargar `Anaconda`
# ```
# wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
# ```
# ![wget-anaconda](images/wget-anaconda.png)
# 
# Copie y pegue la siguiente línea para instalar Anaconda.
# ```
# bash Anaconda3-2020.07-Linux-x86_64.sh
# ```
# 
# Lea y realice las intrucciones que le indique el instalador. Espere unos minutos a que se realice la instalación. Asegúrese de dar `yes` al terminar.
# 
# ![bash-anaconda](images/bash-anaconda.png)
# 
# Copie y pegue las siguientes líneas para instalar los paquetes necesarios para el curso.
# ```
# conda create -n QCI
# conda activate QCI
# conda install psi4 psi4-rt python=3.7 -c psi4
# conda install jupyter notebook matplotlib numpy scipy sympy 
# ```
# 
# Copie y pegue la siguiente línea que hará que cada vez que abra Ubuntu se diriga a la carpeta de su usuario en windows. Reemplace `<Sustituya por su usuario de Windows>` por el nombre de su usuario.
# ```
# echo 'cd /mnt/c/Users/<Sustituya por su usuario de Windows>' >> .bashrc
# ```
# 
# Cierre la ventana, la instalación ha terminado.
# 
# **Cada vez que quiera ejecutar un notebook** deberá colocarlo en su carpeta QCI. Posteriormente abra Ubuntu, puede encontrarlo en el menú de inicio de Windows, deberá escribir las siguientes líneas
# 
# ```
# conda activate QCI
# jupyter notebook
# ```
# 
# Obtendrá una URL que puede copiar y pegar en su navegador
# ![url](images/url.png)
# 
# Se le abrirá una página con sus carpetas, notebooks existentes y la opción de crear nuevos notebooks.
# ![jupyter](images/jupyter.png)
# ````
# 
# ````{admonition} Mac y Linux
# :class: tip
# :class: dropdown
# 
# Vaya a https://www.anaconda.com/products/individual y descargue la versión de Anaconda correspondiente a su sistema operativo.
# 
# ![anaconda-webpage](images/anaconda-webpage.png)
# 
# Si descargó el instalador gráfico de Mac, puede instalarlo como cualquier otro paquete. Si descargó un instalador de línea de comandos, abra una terminal y copie la línea para ir a su carpeta de descargas. Nota: si su computadora está en inglés sustituya Descargas por Downloads
# ```
# cd Descargas
# ```
# 
# Copie y pegue la siguiente línea para instalar Anaconda. Sustituya `<Nombre del instalador>` por el nombre del archivo que descargó, por ejemplo Anaconda3-2020.07-MacOSX-x86_64.sh o Anaconda3-2020.07-Linux-x86_64.sh
# ```
# bash <Nombre del instalador>
# ```
# 
# Lea y realice las intrucciones que le indique el instalador. Espere unos minutos a que se realice la instalación. Asegúrese de dar `yes` al terminar.
# 
# ![bash-anaconda](images/bash-anaconda.png)
# 
# Copie y pegue las siguientes líneas para instalar los paquetes necesarios para el curso.
# ```
# conda create -n QCI
# conda activate QCI
# conda install psi4 psi4-rt python=3.7 -c psi4
# conda install jupyter notebook matplotlib numpy scipy sympy 
# ```
# Cierre la ventana, la instalación ha terminado.
# 
# **Cada vez que quiera ejecutar un notebook** deberá escribir las siguientes líneas
# 
# ```
# conda activate QCI
# jupyter notebook
# ```
# 
# Se le abrirá una página con sus carpetas, notebooks existentes y la opción de crear nuevos notebooks.
# ![jupyter](images/jupyter.png)
# ````
# 
# ### En línea
# 
# `````{admonition} Google Colab
# :class: tip
# :class: dropdown
# 
# Los notebooks pueden ejecutarse en línea desde `Google Colab` al presionar el ícono de la nave y dar clic sobre Colab, este ícono aparece en la parte superior de cada notebook. Los archivos generados en Google Colab se guardan automáticamente en su Google Drive. El único requisito es tener una cuenta de Gmail.
# 
# ![Google Colab example button](images/google_colab.png)
# 
# La ejecución en `Google Colab` de los notebooks que usen Psi4 requiere agregar manualmente lo siguiente como primera celda
# ```
# !wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh
# !bash Miniconda3-4.5.4-Linux-x86_64.sh -bfp /usr/local
# import sys
# sys.path.append('/usr/local/lib/python3.6/site-packages')
# !conda install -c psi4 psi4 -y
# ```
# `````
# 
# ````{admonition} En la misma página web
# :class: tip
# :class: dropdown
# Puede ejecutar los notebooks directamente en el libro haciendo clic sobre el botón.
# 
# ![Live Code](images/live_code.png)
# `````
# 
# Los ejercicios de este libro requieren un conocimiento básico de Python, mismo que puede adquirirse en el primer capítulo del libro. 

# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Python
# 
# Tutorial-Basico
# Tutorial-Avanzado
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Sistemas Sencillos
# 
# Caja
# Tunel_Region_Infinita
# Tunel_Region_Finita
# Oscilador
# Anillo
# Esfera
# Rotor
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Átomo de Hidrógeno
# 
# Radial
# Orbitales
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Métodos Aproximados
# 
# Variacional_Lineal
# PerturbacionesSinDegeneracion
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Sistemas de Muchos Electrones
# 
# Software
# Hartree-Fock
# Moller-Plesset
# DFT
# Dimerizacion
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Matemáticas
# 
# matematicas1
# matematicas2
# ```
# 
# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :numbered: 
# :caption: Anexos
# 
# CajaConTope
# Huckel
# ```
# 
