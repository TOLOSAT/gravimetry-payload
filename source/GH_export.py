"""

@authors:

# =============================================================================
 Information:

    The functions in this script are used to export and create array files

# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt
from time import gmtime, strftime

#import GH_import       as imp
#import GH_convert      as conv
#import GH_generate     as gen
#import GH_solve        as solv
#import GH_displayGeoid as dgeo
#import GH_displaySat   as dsat
#import GH_export       as exp
#import GH_displayTopo  as dtopo
import GH_terminal     as term
#import GH_basemap      as bmp

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
    print(f"Writing \"{title}\" in \"{path}\"")
    length = len(data)
    file = open(f"{path}/{title}", "w+")
    for n in range (len(data)):
        print("Writing", n, "\tof", length-1)
        for m in range(len(data[0])):
            file.write(str(data[n, m]))
            file.write("\t")
        file.write("\n")
    file.close()
    print(f"\n\tDone writing {title}")


# =============================================================================
# TEST FUNCTIONS
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
# MAIN
# =============================================================================
if __name__ == '__main__':

    print("\nGH_export done")

