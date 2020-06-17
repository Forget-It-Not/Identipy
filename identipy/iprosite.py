import os
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from Bio import SeqIO
from Bio.ExPASy import Prosite

def adapt_pattern (pattern):
    '''
    Adapta el patrón desde el formato prosite (input) al formato del módulo re (output)
    '''
    pattern = pattern.replace('x','.')
    pattern = pattern.replace('{','[^')
    pattern = pattern.replace('}',']')
    pattern = pattern.replace('-','')
    pattern = pattern.replace('(','{')
    pattern = pattern.replace(')','}')
    pattern = pattern.replace('<','^')
    pattern = pattern.replace('>','$')
    pattern = pattern[:-1]
    return pattern

#end adapt_pattern()


def make_domain (result_path):
    '''
    Crea un archivo que contiene los dominios encontrados en el multifasta de proteinas filtadas

    Input:
        - result_path: ruta donde esta el multifasta input y donde se guarda el resultado

    Output: archivo que contiene los dominios encontrados con un header que precede a cada proteina
    y los campos: "dominio, accesion, descripcion, patron encontrado". El return se silencia a 'None'
    '''
    #Para hacer el plot posterior de dominios de proteinas se hacen unos archivos temporales
    #que contienen el start y el end y la longitud total. Se almacenan en una carpeta temporal
    os.mkdir('Temporal')

    #Se abren el archivo output y el multifasta que se parsea
    with open(result_path+'/Domains.txt', 'w') as output:
        with open(result_path+'/MultifaFiltered.fasta','r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                #Si una proteína no se matchea con ningún patrón se va a eliminar su
                #archivo temporal, ya que por como se hace el plot añadiría un plot vacío
                #Inicialmente se supone que una proteína no matchea ningun patrón
                has_matched = False
                #Se escribe el header del archivo de dominios:
                output.write('>'+record.id+"\n-------------\n")
                output.write('Nombre dominio\tAccesión\tDescripción\tSecuencia\n')
                #Se abre el archivo temporal:
                with open('Temporal/'+record.id,'w') as temporal_output:
                    #Para cada dominio de prosite se comprueba si existe patrón,
                    #se adapta al formato de re y se intenta matchear
                    #El archivo prosite.dat se tiene que abrir para cada record
                    #ya que al escanear por el, como es un iterator se "consume" (
                    #si se abre al principio funciona para el primer record pero para los demás es como si estuviera vacío)
                    with open('prosite.dat', 'r') as prosite_file:
                        for domain in Prosite.parse(prosite_file):
                            if domain.pattern:
                                adapted_pattern = adapt_pattern(domain.pattern)
                                match = re.search(adapted_pattern, str(record.seq))
                                #Si ha matcheado se escriben los datos correspondientes al resultado
                                #y al output temporal y se cambia has_matched -> True
                                if match:
                                    output.write(domain.name+'\t'+domain.accession+'\t'+domain.description+'\t'+match.group()+'\n')
                                    temporal_output.write(domain.name+'\t'+str(match.start())+'\t'+str(match.end())+'\t'+str(len(record.seq))+'\n')
                                    has_matched = True

                        #Una vez se termina de escanear todos los dominios se añade
                        #un newline extra al output para separar de la proxima proteina
                        output.write('\n')

                #Finalmente una vez se ha acabado con la proteina y se ha cerrado el archivo
                #Si no se ha matcheado se elimina el temporal output creado
                if not has_matched:
                    os.remove('Temporal/'+record.id)

    return None

#end make_domain()


def color_dispenser():
    '''
    Generador para la obtención de un nuevo color (para el diccionario de colores de plot_domains())

    Use:
        generator = color_dispenser
        next(generator)
    '''
    #Obtenemos el color map y discretizamos mediante list_comprehension
    #Estimando hasta 15 colores distintos (puede que haya más dominios
    #y en ese caso se repetirán colores pero si se generan más las dif
    #entre colores son inapreciables)
    cm = plt.get_cmap('gist_rainbow')
    color_list = [cm(1.*i/15) for i in range(15)]
    idx = 0 #indice del color de la lista por el que se va

    #de forma infinita, cuando se llama al next se da el color[idx]
    #y se aumenta el idx para la siguiente llamada. Si va a salir del
    #rango de la lista se reinicia al 0
    while True:
        yield color_list[idx]

        if idx < 14:
            idx += 1
        elif idx == 14: #Pero cuando es 19 se reinicia para que no haya error
            idx = 0

#end color_dispenser()

def plot_domains (temporal_path):
    '''
    Crea los gráficos que muestran los dominios encontrados en las proteínas filtradas.

    Plot donde cada proteina se representa como una caja dentro de la cual hay cajas de
    colores que representan cada dominio encontrado.

    Input:
        - temporal_path: ruta de los archivos que contienen la información de los dominios
        (domain_id match_start match_end total_len_prot)
        siendo que cada archivo corresponde a UNA proteina.
    '''
    #Inicialización de variables y del plot
    height = 3      #Altura a la que se dibuja cada prot
    biggest = 0     #Maximo tamaño encontrado (para designar el xlim)
    prot_num = 0    #Número de proteinas encontradas (se plotea un máximo de 20)

    #Para hacer la leyenda de color de los dominios se mantiene un diccionario
    #de tipo dict[id] :: color, de forma que para cada id de dominio se asigna un color
    color_dict = {}
    color_generator = color_dispenser()    #Generador con el cual se genera un proximo color

    ax = plt.subplot()

    #Para cada proteina se abre el archivo y para cada linea se extraen los datos
    #Almacenados en los ficheros temporales
    for file in os.listdir(temporal_path):
        with open(temporal_path+'/'+file, 'r') as input:
            for line in input.readlines():
                id = line.split('\t')[0]
                start = int(line.split('\t')[1])
                end = int(line.split('\t')[2])
                total = int(line.split('\t')[3])
                #Si el tamaño total encontrado es mayor que el anterior mas grande
                #es el nuevo tamaño más grande
                if total > biggest:
                    biggest = total
                #Si se encuentra un nuevo dominio que no estaba en las keys se añade
                #como key y se le da como valor un nuevo color mediante el generador
                if not id in color_dict.keys():
                    color_dict[id] = next(color_generator)
                #Se dibuja el dominio como una caja con el color adecuado segun el id
                box = mpatch.FancyBboxPatch((start,height), end-start, 3, facecolor=color_dict[id])
                ax.add_patch(box)
            #Finalmente para cada proteina se dibuja también una caja que representa la proteína en total
            box = mpatch.FancyBboxPatch((0,height),total,3, facecolor='yellow', alpha=0.3)
            ax.add_patch(box)
            #Texto que va a la derecha de la caja de proteina y muestra su locus_tag
            ax.text(total+5, height+2.5, file, fontsize=8)
            #Finalmente antes de pasar a la siguiente proteina se aumenta la altura
            #a la que se va a dibujar y se incrementa el num prot
            height += 5
            prot_num += 1
            #Si num prot excede el valor maximo se paran las iteraciones
            if prot_num > 20:
                break

    #Creación de lista de Patches para la leyenda:
    #Para cada asociación dominio :: color encontrada se crea un Patch en una lista
    #y esta lista se puede usar directamente como leyenda
    patchlist = []
    for key in color_dict.keys():
        patchlist.append(mpatch.Patch(color=color_dict[key], label=key))
    plt.legend(handles=patchlist)
    #El limite de x se pone algo superior a biggest para el espacio del texto label
    plt.xlim(-2,biggest+60)
    plt.ylim(-2,102)
    plt.show()

    return None
