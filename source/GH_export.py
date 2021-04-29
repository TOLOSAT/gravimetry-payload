"""
@authors:
# =============================================================================
 Information:
    The functions in this script are used to export and create array files
todo:
    Something is wrong with: Store_temp_GLl
# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
#from numpy import pi, sin, cos
import matplotlib.pyplot as plt
#from time import gmtime, strftime

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
#import GH_terminal     as term
#import GH_harmonics    as harm
#import GH_geoMath      as gmath
#import GH_earthMap     as emap


# =============================================================================
# FUNCTIONS TO STORE FILES
# =============================================================================
def Store_Array(data, title, path="../Rendered"):
    """
    Stores an array into a text file that can later be imported again
    Input:
        data: the array in question
        title: a string, of the desired title for the file.
                Must incluse ".txt"
        path: path in which to store the array
    To import use:
        data = np.loadtxt(title)
    """
#    print(f"Writing \"{title}\" in \"{path}\"")

    file = open(f"{path}/{title}", "w+")
    for n in range (data.shape[0]):
#        print("Writing", n, "\tof", length-1)
        for m in range(data.shape[1]):
            file.write(str(data[n, m]))
            file.write("\t")
        file.write("\n")
    file.close()
#    print(f"\r\tDone writing {title}")


def Store_temp_GLl(G_Grid, G_Long, G_Lat, detail=""):
    """
    Stores arrays into text files for future import
    Should be used with imp.Load_GLl()
    If you want to keep the arrays, move them into the Randered/grid directory,
    or they might get written over
    """
    temp_GLl_path = "../Rendered/grid"
    Store_Array(G_Grid, f"{detail} G_Grid", temp_GLl_path)
    Store_Array(G_Long, f"{detail} G_Long", temp_GLl_path)
    Store_Array(G_Lat,  f"{detail} G_Lat",  temp_GLl_path)



# =============================================================================
# FUNCTIONS FOR FIGURES
# =============================================================================
def Store_Figure(fignum, title, time="", path="../Rendered/images", dpi=500):
    """
    Stores a figure into a .png format
    Input:
        fignum: matplotlib figure number
        title: title image name.
        path: image path location
        dpi: pixels per inch density
    """
    plt.figure(fignum)
#    mng = plt.get_current_fig_manager()
#    mng.window.showMaximized()
#    plt.show()
    file_name = f"{path}/{time} {title}.png"
    plt.savefig(file_name, dpi=dpi)

# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_store_temp():
    A = np.ones((1,5))
    B = np.ones((2,5))*2
    C = np.ones((5,2))*3
    Store_temp_GLl(A, B, C)



# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':

#    TEST_store_temp()



    '''
    Store_temp_GLl(G_Grid, G_Long, G_Lat, "TESTrr0")
    '''



    print("\nGH_export done")