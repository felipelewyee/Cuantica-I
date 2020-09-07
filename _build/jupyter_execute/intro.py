# La Química Cuántica a tu alcance

El presente libro digital tiene como objetivo ser un manual de prácticas para la asignatura de [`Química Cuántica I (clave 1404)`](https://quimica.unam.mx/wp-content/uploads/2017/03/1404qcuanticauno.pdf) de la Facultad de Química (UNAM). La última sección puede considerarse una parte avanzada disponible de manera opcional para quien desee aprender más, y puede ser de utilidad para la asignatura de [`Química Cuántica II (clave 0038)`](https://quimica.unam.mx/wp-content/uploads/2017/06/0038qcuanticados_Q.pdf). El libro digital es un esfuerzo compartido de Juan Felipe Huan Lew Yee, Eduardo Barrios Vargas y Jorge Martín del Campo Ramírez, en conjunto con los profesores del Departamento de Física y Química Teórica. Toda retroalimentación es bienvenida, y se le invita a compartirla a los correos:

M. en C. Juan Felipe Huan Lew Yee: [felipe.lew.yee@gmail.com](mailto:felipe.lew.yee@gmail.com)

Dr. Jorge Martín del Campo Ramirez: [jormacara@gmail.com](mailto:jormacara@gmail.com)

Dr. José Eduardo Barrios Vargas: [j.e.barrios@gmail.com](mailto:j.e.barrios@gmail.com)

Se ha buscado reducir lo más posible los antecedentes matemáticos requeridos, por lo que en la mayoría de los ejercicios se plantea el problema y se pasa directamente a la solución sin realizar el desarrollo. Se asume que las ecuaciones han sido previamente revisadas en clase, sin embargo, en muchos notebooks pueden encontrarse partes clave del desarrollo matemático ocultas en una celda como la siguiente
```{admonition} Inserto matemático: Celda con contenido matemático oculto
:class: dropdown

La lectura de este apartado matemático oculto es opcional si ya se conoce la deducción previa de las ecuaciones, pero aquí está por si alguien lo necesita. Aquí van algunas funciones trigonométricas

$$
y = sin(x)
$$

$$
y = cos(x)
$$

$$
y = tan(x) = \frac{sin(x)}{cos(x)}
$$
```

Para los alumnos del curso de `Química Cuántca I` se recomienda la lectura en el orden de los capítulos, el cual puede seguirse utilizando la barra de navegación que aparece del lado izquierdo en computadoras, o con el menú emergente en dispositivos móviles. Este libro permite interactuar con el contenido, para lo cual se recomienda alguna de las dos opciones siguientes:

````{admonition} Ejecución de los notebooks (en línea)
:class: tip
Los notebooks pueden ejecutarse en línea desde **Google Colab** al presionar el ícono de la nave y dar clic sobre Colab, este ícono aparece en la parte superior de cada notebook. Los archivos generados en Google Colab se guardan automáticamente en su Google Drive. El único requisito este tener una cuenta asociada a **Google**.

![Google Colab example button](images/google_colab.png)
````

````{admonition} Ejecución de los notebooks (en computadora)
:class: tip
Puede descargar los *.ipynb y ejecutarlos directamente en su computadora dando click sobre el botón con la flecha apuntando hacia abajo. Para ejecutarlos localmente es su computadora se requiere una instalación previa de Python. En caso de no tenerla puede instalar [Anaconda](https://www.anaconda.com/products/individual). Una vez descargado puede interactuar con el notebook utilizando Jupyter Notebook.

En caso de tener una instalación local de Anaconda en Linux verifique que tenga instaladas las librerías utilizadas en los notebooks. Dichas librerías las puede instalar desde la terminal utilizando las siguientes líneas:
```
conda install matplotlib numpy scipy sympy 
```
Después, puede ejecutar un notebook con Jupyter
```
jupyter notebook NOMBRE.ipynb
```

Algunos ejercicios utilizan el software psi4, este se puede instalar desde Anaconda con
```
conda install psi4 psi4-rt python=3.7 -c psi4
```
````

Los ejercicios de este libro requieren un conocimiento básico de Python, mismo que puede adquirirse en el primer capítulo del libro. 

## Uso del manual

A lo largo de los ejercicios encontrará celdas con el código correspondiente, y celdas vacías con una nota para indicar que debe ser llenada siguiendo alguna instrucción, **como la siguiente**

# Llene la siguiente celda con su código

La primera vez que se presente una instrucción se mostrará directamente para copiar y pegar en la celda, por ejemplo, el siguiente código muestra como realizar una impresión
```
print("Bienvenidos al manual de Química Cuántica I")
```

Posteriormente, cuando se solicite realizar algo usando código ya presentado, se mostrará la respuesta oculta para permitir que se intente resolver el ejercicio en la celda vacía por cuenta propia. Para revelar la respuesta hay que presionar en el botón con el símbolo + como en la siguiente celda

print("Bienvenidos al manual de Química Cuántica I")


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Python

Tutorial-Basico
Tutorial-Avanzado
```


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Sistemas Sencillos

Caja
CajaConTope
Anillo_Esfera_Rotor
Oscilador
```


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Átomo de Hidrógeno

Radial
Orbitales
```


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Métodos Aproximados

Variacional_Lineal
PerturbacionesSinDegeneracion
```


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Sistemas de Muchos Electrones

Software
Hartree-Fock
Moller-Plesset
DFT
Dimerizacion
```


```{toctree}
:hidden:
:titlesonly:
:numbered: 
:caption: Matemáticas

matematicas1
matematicas2
```
