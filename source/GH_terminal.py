"""
@authors:
# =============================================================================
 Information:
    The functions in this script are used to write text to the terminal
    with a nice printout
# =============================================================================
"""

# =============================================================================
# LIBRARIES
# =============================================================================
import time

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
# FUNCTIONS TO PRINT TEXT
# =============================================================================
def printProgressBar (iteration, total, length=20, decimals=1):
    """
    Prints a progress bar. Original code from:
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    Call in a loop to create a progress bar in the terminal
    Input:
        iteration: current iteration (Int)
        total: total iterations (Int)
        length: character length of bar
        decimals: positive number of decimals in percent complete
    """
    prefix="Progress:" #
    suffix="Complete" #

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = "█" * filledLength  +  '-' * (length-filledLength)

    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end = "\r")

    # Print new line when complete
    if iteration == total:
        print("\n")


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
def TEST_ProgressBar():
    # A List of Items
    items = list(range(0, 20))
    l = len(items)

    printProgressBar(0, l)
    for i, item in enumerate(items):
        # Do stuff...
        time.sleep(1)
        # Update Progress Bar
        printProgressBar(i + 1, l)


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    TEST_ProgressBar()

    print("\nGH_export done")