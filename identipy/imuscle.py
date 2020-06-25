import os
from subprocess import call
from Bio import Phylo

def build_muscle (result_path):
    '''
    Crea el árbol filogenético de un grupo de proteinas filtradas

    Input:
        - result_path: ruta donde está el archivo input y
          donde dejar el resultado

    Output: genera el archivo newick del árbol filogenético.
            El return se silencia a 'None'.
    '''
    #Se abre devnull para silenciar el mensaje de stderr que genera muscle
    #(que no es un error en realidad)
    with open(os.devnull,'w') as devnull:
        call(['muscle', '-in',result_path+'/MultifaFiltered.fasta', '-out', \
              result_path+'/MuscleAlign'], stderr=devnull)

    return None

#end build_muscle()


def plot_muscle (result_path):
    '''
    Dibuja el árbol filogenético de un archivo .nw
    '''
    #De nuevo se abre devnull para silenciar muscle
    with open(os.devnull, 'w') as devnull:
        call(['muscle', '-maketree', '-in', result_path+'/MuscleAlign', '-out',\
              result_path+'/MuscleTree.nw'], stderr=devnull)
              
        #El Árbol se parsea y se gráfica con el Phylo de Biopython
        tree = Phylo.read(result_path+'/MuscleTree.nw', 'newick')
        Phylo.draw(tree)
    return None

#end plot_muscle()
