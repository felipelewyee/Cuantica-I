#!/usr/bin/env python
# coding: utf-8

# # La Química Cuántica a tu alcance

# El presente libro digital tiene como objetivo ser un manual de prácticas para la asignatura de [`Química Cuántica I (clave 1404)`](https://quimica.unam.mx/wp-content/uploads/2017/03/1404qcuanticauno.pdf) de la Facultad de Química (UNAM). El libro digital es un esfuerzo compartido de Juan Felipe Huan Lew Yee, José Eduardo Barrios Vargas y Jorge Martín del Campo Ramírez, en conjunto con los profesores del Departamento de Física y Química Teórica. Toda retroalimentación es bienvenida, y se le invita a compartirla a los correos:
# 
# M. en C. Juan Felipe Huan Lew Yee: [felipe.lew.yee@quimica.unam.mx](mailto:felipe.lew.yee@quimica.unam.mx)
# 
# Dr. Jorge Martín del Campo Ramírez: [jormacara@quimica.unam.mx](mailto:jormacara@quimica.unam.mx)
# 
# Dr. José Eduardo Barrios Vargas: [jebarrios@quimica.unam.mx](mailto:jebarrios@quimica.unam.mx)
# 
# ## Uso del libro digital
# 
# Se recomienda la lectura en el orden de los capítulos, el cual puede seguirse utilizando la barra de navegación que aparece del lado izquierdo en computadoras, o con el menú emergente en dispositivos móviles.
# 
# Se ha buscado reducir lo más posible los antecedentes matemáticos requeridos, por lo que en la mayoría de los ejercicios se plantea el problema y se pasa directamente al análisis de la solución. Aunque se asume que las ecuaciones han sido previamente revisadas en clase, el desarrollo matemático correspondiente se encuentra en celdas ocultas como la siguiente:
# ```{admonition} Inserto matemático: Celda con contenido matemático oculto
# :class: dropdown
# 
# La lectura de este apartado matemático oculto es opcional si ya se conoce la deducción previa de las ecuaciones, pero aquí está por si alguien lo necesita. Aquí van algunas funciones trigonométricas
# 
# $$
# y = \sin(x)
# $$
# 
# $$
# y = \cos(x)
# $$
# 
# $$
# y = \tan(x) = \frac{\sin(x)}{\cos(x)}
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
# `````{admonition} Google Colab
# :class: tip
# :class: dropdown
# 
# Los notebooks pueden ejecutarse en línea desde `Google Colab` al presionar el ícono de la nave y dar clic sobre Colab, este ícono aparece en la parte superior de cada notebook. Los archivos generados en Google Colab se guardan automáticamente en su Google Drive. El único requisito es tener una cuenta de Gmail.
# 
# ![Google Colab example button](images/google_colab.png)
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

# In[ ]:




