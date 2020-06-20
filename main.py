#!/usr/bin/python3
help_message = '''
main.py
-------
Script para el análisis de ensamblados genómicos en base al paquete  identipy:
    - Obtencion de proteínas de ensamblados genómicos que son homólogas a una o varias proteinas query
    - Alineamiento y obtención del árbol filogenético de dichas proteinas
    - Obtención de la lista de dominios presentes en cada proteina a partir de la base de datos de Prosite

Instrucciones de uso
---------------------
    1. Seleccionar/crear una carpeta en la que se va a realizar el análisis
    2. Crear en su interior una carpeta que contiene el grupo de genbanks de ensamblados genómicos que se van a usar (database)
    3. Crear también una carpeta que contiene el grupo de proteinas query en formato Fasta de a 1 proteína por archivo (query)
    4. Descargar y colocar también el archivo prosite.dat (ftp://ftp.expasy.org/databases/prosite/prosite.dat) en la carpeta
    5. Colocar finalmente el paquete de identipy y el script
    6. Ejecutar el script

Ejecución del script: \"main.py query_path database_path [evalue=num] [cov=num] [iden=num]\"
    - query_path: ruta del directorio donde están contenidas todas las proteinas query (una query por archivo en formato fasta)
    - database_path: ruta del directorio donde están contenidos todos los genbank referentes a los ensamblados genómicos

    Opcionales
    - evalue=num: num se sustituye por el número que se quiere aplicar como máximo de evalue a blastp (ej. evalue=0.01)
    - cov=num: num se sustituye por el % que se quiere aplicar como mínimo de cov a blastp (ej. cov=70)
    - iden=num: num se sustituye por el % que se quiere aplicar como mínimo de id a blastp (ej. id=70)

Archivos resultado
------------------
    - MultifaDB.fasta: multifasta con todas las secuencias del ensamblado
    - ID_Organism_table.csv: tabla con las equivalencias entre @id y organismo
    - Result_$query$: directorio que contiene el resto de archivos para la query $query$
        - Blast_result: resultado del blastp
        - MultifaFiltered: secuencias filtradas por el blastp
        - MuscleAlign: alineamiento de muscle de las secuencias filtradas
        - MuscleTree: árbol filogenético con las secuencias filtradas
        - Domains: tablas que contienen los datos de dominios encontrados

Gráficos durante la ejecución
-----------------------------
    - Gráfico scatter que resume de los resultados del blast
    - Gráfico de los árboles filogenéticos obtenidos
    - Gráfico de dominios encontrados y su posición en las proteínas
'''


#IMPORTACIÓN DE MODULOS
#######################
import os
import sys
import shutil
from subprocess import call
from Bio import SeqIO
import identipy as id


#Códigos ANSI para escritura en color utilizados en los print
#Para que una parte este en color se coloca entre el color y normal (color+string+normal)
header = '\033[95m'     #Morado
error = '\033[91m'      #Rojo
success = '\033[92m'    #Verde
question = '\033[93m'   #Amarillo
normal = '\033[0m'      #Texto normal


#CONTROL DE ARGUMENTOS
######################
if not len(sys.argv) >= 3: #Si no hay más 2 argumentos (3 contando el nombre) muestra la ayuda
    print(help_message)
    exit(1)
else: #Si hay el número correcto se extraen los path de querys y database y se comprueba que existen
    querys_path, database_path = sys.argv[1:3]
    if not os.path.isdir(querys_path):
        raise Exception('La ruta de los archivos query no ha sido encontrada')
    if not os.path.isdir(database_path):
        raise Exception('La ruta de los archivos database no ha sido encontrada')

    #El archivo prosite.dat DEBE estar en la carpeta donde se ejecuta el script. Se comprueba que esto es así
    if not os.path.isfile('prosite.dat'):
        raise Exception(error+'¡ERROR: falta el archivo prosite.dat en la carpeta en la que se hace el análisis!'+normal)

#El evalue, cov e id se inicializan
eval = 0.01
cov = 50
iden = 50

#Los argumentos opcionales son todos los que hay menos el nombre del script y las 2 rutas
extra_args = len(sys.argv)-3
for i in range(3,3+extra_args): #Para cada uno de estos argumentos adicionales
    #se separa por el = el nombre del argumento y su valor
    arg, value = sys.argv[i].split('=')
    #Se intenta convertir a numérico el valor, lo que confirma que es un número
    try:
        value = float(value)
    #Si la conversión numérica ha dado error se avisa pero se mantiene el valor por defecto
    except:
        print(error+'¡ERROR: has aportado un '+arg+' para el blast que no es de tipo numérico'+normal+ \
        ' por lo que ha sido omitido, en su lugar se usará el valor establecido por defecto\n')
    #Si el valor era numérico se establece su valor en la variable adecuada
    else:
        if arg == 'evalue':
            eval = value
        elif arg == 'cov':
            cov = value
        elif arg == 'iden':
            iden = value
        #Si el argumento indicado no existe se avisa pero se ignora
        else:
            print(error+'¡ERROR: has aportado el argumento opcional '+arg+' que no existe!' +normal+\
            ' por lo que ha sido ignorado\n')

#SCRIPT
#######
print('/****************************************************************************************************************************')
print('* '+header+'Ejecutando análisis de identipy en: '+os.getcwd()+normal)
print('\\****************************************************************************************************************************')

print('\n> Creando base de datos multifasta MultifaDB.fasta de los ensamblados genómicos...')
print('-----------------------------------------------------------------------------------')

#Creación del MultifaDB.fasta a partir de las secuencias de la database
id.make_multifa(database_path)
print('\t'+success+'{+}'+normal+' ¡La base de datos se ha creado exitosamente!')

for file in os.listdir(querys_path):
    #Para cada proteina query se hace un análisis de blast, muscle y prosite
    path = querys_path+"/"+file
    filename = file.split('.')[0] #Se obtiene el nombre separando por el punto y quedandose con la parte anterior

    print('\n> Ejecutando análisis para la query '+filename+'...')
    print('------------------------------------------------------------')

    #Control de argumentos para query:
    try: #Se intenta parsear por seqio
        SeqIO.read(path, 'fasta')
    except:
        print('\t'+error+'{x}'+normal+' El archivo query '+file+' no está en formato fasta y se ha omitido para el análisis')
    else:
        #Si no da error se procede a intentar crear los directorios de resultado
        #Si el directorio ya existe se pregunta si se desea sobreescribir el
        if os.path.isdir('Result_'+filename):
            print('\t'+error+'{x}'+normal+' ERROR: Ya existe una carpeta de resultados para la query '+filename)
            user_choice = input('\t'+question+'{?}'+normal+' Quieres sobreescribirla? [y/n]: ')
            if user_choice in ['Y','y','yes']:
                #Si se desea sobreescribir se borra el directorio y sus contenidos
                #y se crea de nuevo la carpeta
                shutil.rmtree('Result_'+filename)
                os.mkdir('Result_'+filename)
            else:
                #Pero si no se desea sobreescribir se pasa a la siguiente query
                continue
        #Si el directorio no existía simplemente se crea
        else:
            os.mkdir('Result_'+filename)

        #Se obtienen las IDs de las proteínas filtradas por blast
        #Durante esta llamada también se pregunta si se desea visualizar la gráfica de blast
        #(ver código de la función)
        hits = id.blastp_hits(path, 'Result_'+filename, eval, cov, iden)
        if len(hits) > 0: #Si se obtuvo algún hit se realiza el análisis

            #Si se ha filtrado alguna proteina se da la opción de ver el gráfico resumen
            #Durante el blast se ha generado el temporal 'blast_temporal' que se usa para el plot
            user_choice = input('\t\033[93m{?}\033[0m ¿Quieres observar un gráfico resumen del resultado de blast? [y/n]: ')
            if user_choice in ['Y','y','yes']:
                id.blast_plot('blast_temporal')
            #Una vez se pasa el plot se haya hecho o no se elimina el temporal
            os.remove('blast_temporal')

            #Se genera la base de datos filtrada (Result_$query$/MultifaFiltered)
            id.filter_database('Result_'+filename, path, hits)
            print('\t'+success+'{+}'+normal+' ¡La base de datos se ha filtrado exitosamente!')

            #Se crea el alineamiento de muscle. Una vez creado se pregunta si se desea visualizar la gráfica
            id.build_muscle('Result_'+filename)
            print('\t'+success+'{+}'+normal+' ¡El árbol filogenético se ha creado exitósamente!')
            user_choice = input('\t'+question+'{?}'+normal+' ¿Quieres observar un gráfico del árbol filogenético? [y/n]: ')
            if user_choice in ['Y','y','yes']:
                id.plot_muscle('Result_'+filename)

            #Se crea el archivo que contiene los dominios encontrados en cada proteina
            #Una vez creado se pregunta si se desea visualizar la gráfica
            id.make_domain('Result_'+filename)
            print('\t'+success+'{+}'+normal+' ¡El archivo de dominios ha sido creado exitosamente!')
            user_choice = input('\t'+question+'{?}'+normal+' ¿Quieres observar un gráfico de los dominios encontrados? [y/n]: ')
            if user_choice in ['Y','y','yes']:
                #Durante make_domain se crea el archivo output Result_$query$/Domains
                #Y un grupo de archivos en la carpeta Temporal que se usan para el plot
                id.plot_domains('Temporal')
            #Una vez pasa el plot se haya hecho o no se elimina el temporal
            shutil.rmtree('Temporal')

        else: #(viene de if len(hits)) : Si no se obtuvieron hits se omite el análisis para esta query
            print('\t'+error+'{x}'+normal+' No se encontró ningún hit. Análisis abortado.')
