import os
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

def make_multifa (database_path):
    '''
    Crea el archivo multifasta con las secuencias de proteína de un grupo de ensamblados genómicos de tipo GenBank.
    También comprueba que las secuencias están en el formato GenBank.
    Si no lo están se omiten con un aviso al usuario en vez de cancelar todo el script.

    Input:
        - database_path: ruta de la carpeta que contiene los archivos

    Output: genera el archivo MultifaDB. El return se silencia a 'None'
    '''
    #Se abre el archivo donde se guarda el resultado
    #y el archivo donde se guarda la tabla de equivalencias id - organismo
    with open('MultifaDB.fasta', 'w') as output:
        with open('ID_Organism_table.csv', 'w') as tabla:
            #Cada archivo del directorio se abre y parsea
            for file in os.listdir(database_path):
                with open(database_path+'/'+file, 'r') as handle:
                    genbank = SeqIO.parse(handle, 'genbank')
                    #El control de argumentos es complicado porque SeqIO parsea el archivo
                    #sin producir error aunque no sea GenBank (no se puede usar try/except)
                    any = False #Se asume inicialmente que no tiene ningún record
                    for record in SeqIO.parse(handle, 'genbank'):
                        #Si tiene algún record se cambia any -> True y se asume que es GenBank
                        any = True
                        #Se añade la entrada correspondiente en la tabla ID-organismo
                        tabla.write(record.id+'\t'+record.annotations['organism']+'\n')
                        #Para cada feature se prueba a obtener el locus_tag del gen y la secuencia de proteinas
                        #Si se consigue se escribe en el archivo resultado y si no se pasa al siguiente archivo
                        for feature in record.features:
                            try:
                                feature.qualifiers['locus_tag'][0]
                                feature.qualifiers['translation'][0]
                            except:
                                pass
                            else:
                                output.write(">"+feature.qualifiers['locus_tag'][0]+'@'+record.id+"\n") #Se añade al final @id para identificar el organismo
                                output.write(str(feature.qualifiers['translation'][0]) +"\n\n")
                    #Si no se había encontrado ningún record se avisa al usuario de que el tipo de archivo es incorrecto
                    #pero se omite para continuar el script
                    if not any:
                         print('\t\033[91m{x}\033[0m El archivo database '+file+' no está en formato'+ \
                         ' genbank (o no contiene ninguna entrada) y se ha omitido para el análisis.')

    return None

#end make_multifa()


def blastp_hits (query_path, result_path, eval):
    '''
    Obtiene la lista de proteinas de un multifasta que han sido filtradas
    al hacer un blasteo frente a una query.

    Input:
        - query_path: ruta de la query para la cual se filtra
        - result_path: ruta del fichero de resultados para la query adecuada
        - eval: evalue usado en la llamada a blast

    Output: lista que contiene los nombres de las proteinas filtradas
    '''
    #Se abre el archivo donde guardar el resultado de blast
    #Y se realiza la llamada que se dirige a ese archivo
    with open('blast_temporal','w') as blast_file:
        call(['blastp','-query',query_path,'-subject', 'MultifaDB.fasta', '-evalue', eval, '-outfmt', '6 qseqid sseqid qcovs pident sseq'], stdout=blast_file)

    #Posteriormente se vuelve a abrir el archivo pero como lectura,
    #Se separan las líneas por tabulación y se recoge el 2do campo (nombres de genes)
    #Que se almacena en la lista de IDs
    with open('blast_temporal','r') as blast_file:
        hits = []
        lines = blast_file.readlines()
        for line in lines:
            hits.append(line.split('\t')[1])

    #Si se ha filtrado alguna proteina se da la opción de ver el gráfico resumen
    if len(hits) > 0:
        user_choice = input('\t\033[93m{?}\033[0m ¿Quieres observar un gráfico resumen del resultado de blast? [y/n]: ')
        if user_choice in ['Y','y','yes']:
            blast_plot()

    #Finalmente se crea el archivo final de resultado en la carpeta correspondiente
    #El archivo temporal carece de header porque facilita el tratamiento, y ahora se le añade
    with open(result_path+'/Blast_result', 'w') as blast_final:
        call(['echo', 'Query\tSubject\t% Coverage\t% Identity\t Sequence'], stdout=blast_final)
        call(['cat', 'blast_temporal'], stdout=blast_final)
    os.remove('blast_temporal')

    return hits

#end blastp_hits()


def filter_database (result_path, query_path, hits):
    '''
    Crea el archivo multifasta con las secuencias filtradas por blastp

    Input:
        - result_path: ruta de la carpeta donde se debe guardar el archivo output
        - query_path: ruta de la query usada para el filtrado
        - hits: lista de nombres de proteinas filtradas

    Output: genera el archivo MultifaFiltered. El return se silencia a 'None'
    '''
    #Se abre el archivo output y primero se introduce la query usada
    with open(result_path+'/MultifaFiltered.fasta', 'w') as output:
        with open(query_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                output.write(">"+record.id+"\n")
                output.write(str(record.seq)+"\n\n")

        #Posteriormente se abre y parsea el multifasta
        with open('MultifaDB.fasta','r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                #Si el nombre de la secuencia es un hit -> output
                if str(record.id) in hits:
                    output.write(">"+record.id+"\n")
                    output.write(str(record.seq)+"\n\n")
    return None

#end filter_database()


def blast_plot ():
    '''
    Dibuja un gráfico resumen del resultado de blast en forma de un scatterplot 2D
    con los valores de identidad y coverage obtenidos
    '''
    #Se inicializan variables cov e id donde se guardan los valores respectivos
    cov = []
    id = []

    #Se abre el archivo de blast y se sacan del 3er y 4to campo (ver llamada a blast)
    #los valores de cov e id y se van añadiendo a la lista
    with open('blast_temporal','r') as blast_file:
        lines = blast_file.readlines()
        for line in lines:
            str_cov = line.split('\t')[2]
            str_id = line.split('\t')[3]
            cov.append(float(str_cov))
            id.append(float(str_id))

    #Finalmente la lista se convierte en array para el plotting
    cov = np.array(cov)
    id = np.array(id)

    #Códigos del plot
    #################
    #Scatter plot con tamaño grande y transparencia
    #para que al acumularse varios puntos haya más color
    plt.scatter(cov,id,s=200, alpha=1/len(cov))
    plt.xlabel('Cobertura de alineamiento (%)')
    plt.ylabel('Identidad de alineamiento (%)')
    plt.title('Distribución de cobertura e identidad de los alineamientos')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(5,90,'Nº of hits: '+str(len(cov)), bbox=props) #Se dibuja una caja que indica el nº de hits
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.show()

    return None

#end blast_plot()
