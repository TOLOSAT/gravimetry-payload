"""

@authors:

# =============================================================================
 Information: 

    The functions in this script are used to write text to the terminal 
    with a nice layout
        
# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import numpy as np
from numpy import pi, sin, cos
import time


#import GH_import      as imp
#import GH_convert     as conv
#import GH_generate    as gen
#import GH_solve       as solv
#import GH_displayCoef as dcoef
#import GH_displaySat  as dsat

# =============================================================================
# FUNCTIONS TO PRINT TEXT
# =============================================================================
def printProgressBar (iteration, total, decimals=1, length=15):
    """
    Prints a progress bar. Original code from: 
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
        
    Call in a loop to create a progress bar in the terminal
    Input:
        iteration: current iteration (Int)
        total: total iterations (Int)
        decimals: positive number of decimals in percent complete (Int)
        length: character length of bar 
        
    """    
    prefix="Progress: " #
    suffix="Complete" #
    fill="█"        # bar fill character
    printEnd="\r"   # end character (e.g. "\r", "\r\n")
    
    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_ProgressBar():
    # A List of Items
    items = list(range(0, 57))
    l = len(items)
    
    # Initial call to print 0% progress
    printProgressBar(0, l)
    for i, item in enumerate(items):
        # Do stuff...
        time.sleep(0.1)
        # Update Progress Bar
        printProgressBar(i + 1, l)
        
        
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    TEST_ProgressBar()
    
    print("\nGH_export done")

