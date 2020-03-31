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

#import GH_import      as imp
#import GH_convert     as conv
#import GH_generate    as gen
#import GH_solve       as solv
#import GH_displayCoef as dcoef
#import GH_displaySat  as dsat

# =============================================================================
# FUNCTIONS TO STORE FILES
# =============================================================================
def Store_array(data, title, path=""):  
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
    length = len(data)
    file = open(title, "w+")
    for n in range (len(data)):
        print("Writing", n, "\tof", length-1)
        for m in range(len(data[0])):
            file.write(str(data[n, m]))
            file.write("\t")
        file.write("\n")
    file.close()
    print("\n\tWriting done\n")


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
  
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
 
    print("\nGH_export done")

