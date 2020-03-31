"""

@authors:

# =============================================================================
 Information: 
    
    The purpose of this script is to display various graphs, plots to show:
        "True" and "solved" Geoids
        Differences in coefficients
        
        
# =============================================================================
"""
# =============================================================================
# LIBRARIES
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np

#import GH_convert     as conv
import GH_import      as imp
#import GH_generate    as gen
#import GH_solve       as solv
#import GH_displaySat  as dsat


# =============================================================================
# DISPLAY FUNCTIONS  
# =============================================================================
def Plot_Array_Diff(HS_nm_slv, HC_nm_slv, fig_num = 6):
    print("plotting coeff difference")
    
    
    #resize the official coef
    HC, HS = imp.Fetch_Coef()
    HS_nm_sz = HS[:len(HS_nm_slv), :len(HS_nm_slv)]
    HC_nm_sz = HC[:len(HC_nm_slv), :len(HC_nm_slv)]

    #subtract calculated coeffs
    HS_nm_sz -= HS_nm_slv
    HC_nm_sz -= HC_nm_slv
    
    fig_HC = plt.figure(fig_num)
    plt.clf()
    plt.suptitle("Harmonic coeff difference between official and solved; degree: "+str(len(HS_nm_sz)-1))
    
    for n in range (0, len(HC_nm_sz)):
        Ms_n = np.arange(0, n+1)
        
        HC_ni = HC_nm_sz[n, :n+1]
        HS_ni = HS_nm_sz[n, :n+1]
        
        plt.subplot(211)
        plt.plot(Ms_n, HC_ni,'-*', label='n='+str(n))
        
        plt.subplot(212)
        plt.plot(Ms_n, HS_ni,'-*', label='n='+str(n))
    
    plt.subplot(211)
    plt.ylabel("COSINE coeff diff")
    plt.grid(True)
#    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HC_nm")
    plt.legend(loc = 'upper right', title = 'Degree n', fontsize = 5)
    
    plt.subplot(212)
    plt.ylabel("SINE coeff diff")
    plt.grid(True)
    plt.xlabel("order m of derivation (log)")
#    plt.ylabel("value of HS_nm")
#    plt.legend(loc = 'lower right', title = 'Degree n', fontsize = 5)
    
    plt.show()
    


# =============================================================================
# TEST FUNCTIONS
# =============================================================================
    
# =============================================================================
# MAIN 
# =============================================================================
if __name__ == '__main__':
    
    print("\nGH_displayCoef done")

