#!/usr/bin/python3
'''
IDENTIPY - 5.0
--------------
> Paquete para el análisis de ensamblados genómicos:
    - Alineamiento de las proteínas del ensamblado con una o varias querys
    - Alineamiento por muscle y creación del árbol filogenético
    - Búsqueda de dominios de la base de datos Prosite

> Uso recomendado: uso y control del paquete mediante el script main.py

> Módulos contenidos en este script:
    - iblast:
        + Generación de base de datos multifasta
        + Filtrado de la base de datos con querys
        + Graficado de resumen del blast

    - imuscle:
        + Generación de árbol filogenético
        + Visualización del árbol filogenético

    - iprosite:
        + Obtención de los dominios de las proteínas filtradas
        + Visualización de los dominios encontrados
'''

#Al importar el paquete se importan las funciones de los de módulos
from .iblast import *
from .imuscle import *
from .iprosite import *
