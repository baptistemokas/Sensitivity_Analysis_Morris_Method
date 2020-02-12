from SALib.sample import saltelli
from SALib.analyze import sobol
import sys

from SALib.analyze import morris
from SALib.sample.morris import sample
from SALib.test_functions import Sobol_G
from SALib.util import read_param_file
from SALib.plotting.morris import horizontal_bar_plot, covariance_plot, \
    sample_histograms
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



def ET(X):
    # formule d'origine 0.0031*X[:,0]*(X[:,1]+209)*(X[:,2]*(X[:,2]+15))**-1)
    # column 0 = C, column 1 = R, column 2 = t
    ############## CATEGORIE ORGANISATION INTERNE #################################
    x = np.linspace(0, 100, num=11, endpoint=True)
    y = np.cos(-x**2/9)
    f = interp1d(x, y, kind='cubic')
    
    int1 = (2*f(X[:,0])+1.9*f(X[:,1])+1.2*f(X[:,2])+1.5*f(X[:,3]))/6.7  # int1a - int1b -  int1c -  int1d
    
    int2a= (1.3*f(X[:,4])+1.1*f(X[:,5]))/2.4 # int2a1 - int2a2
    int2b= (2*f(X[:,6])+1*f(X[:,7])+3*f(X[:,8]))/6  #int2b1 - int2b2 - int2b3
    int2= (4*int2a+1*int2b)/5
    
    int3a= (1*f(X[:,9])+1.9*f(X[:,10]))/2.9 #int3a1 int3a2
    int3b= (1*f(X[:,11])+2*f(X[:,12])+4*f(X[:,13]))/7 #int3b1 int3b2 int3b3
    
    #### ....
         ....
         ....

         ....
         ....
         ....


    net1= net1a+net1b+net1c+net1d+net1e+net1f+net1g+net1h+net1i+net1j+net1k+net1l+net1m+net1n+net1o+net1p+net1q+net1r+net1s+net1t+net1u+net1v+net1w+net1x+net1y+net1z  
    net2a = f(X[:,522])+f(X[:,523])+f(X[:,524])+f(X[:,525])+f(X[:,526])+f(X[:,527])+f(X[:,528])
    net2= net2a+f(X[:,529])+f(X[:,530])
    catNET = (net1+net2)/2

    ######################################################
    ######################################################

#Int=(catINT+catNET+catSER+catMAJ+catACC+catMAN+catBCK+catMDP+catAPP+catEX+catPHY)

#return    int rssi jur op phy ex  mdp  bck man app acc maj ser net 

#######Créationsdesscénarios ##########
    AS1=(7*catINT+4*catRSSI+6*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/17
    AS2=(1*catINT+4*catRSSI+2*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/7
    AS3=(5*catINT+5*catRSSI+1*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/11
    AS4=(5*catINT+5*catRSSI+1*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/11
    
        #### 
         ....
         ....
         ....

         ....
         ....
         ....


    MO7=(8*catINT+3*catRSSI+1*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/12
    MO8=(1*catINT+2*catRSSI+4*catJUR+7*catOP+8*catPHY+3*catEX+5*catMDP+5*catBCK+2*catMAN+5*catAPP+5*catACC+3*catMAJ+2*catSER+3*catNET)/7
    
    s1  =AS1+FP1+MO1
    s2  =AS1+FP1+MO2
    s3  =AS1+FP1+MO3

        #### 
         ....
         ....
         ....

         ....
         ....
         ....

    s315=AS5+FP8+MO3
    s316=AS5+FP8+MO4
    s317=AS5+FP8+MO5
    s318=AS5+FP8+MO6
    s319=AS5+FP8+MO7
    s320=AS5+FP8+MO8
    
    D=(s1+s2+        #### ....   +s315+s316+s317+s318+s319+s320)/320
    I=(s1+s2+        #### ....   +s315+s316+s317+s318+s319+s320)/320
    C=(s1+s2+        #### ....   +s315+s316+s317+s318+s319+s320)/320
    P=(s1+s2+        #### ....   +s315+s316+s317+s318+s319+s320)/320
    
    
   
    return(8*D+7*I+0.7*C+2*P)
    
              


problem = {
    'num_vars': 533,
    'names': ['int1a',                   #[:,0]
    'int1b',                             #[:,1]
    'int1c',                             #[:,2]
    'int1d',                             #[:,3]
           
    'int2a1',                            #[:,4]
    'int2a2',                            #[:,5]
           
    'int2b1',                            #[:,6]
    'int2b2',                            #[:,7]
    'int2b3',                            #[:,8]
           
         #### 
         ....
         ....
         ....

         ....
         ....
         ....       

    'net1r1',                           #[:,506]
    'net1r2',                           #[:,507]

    'net1s1',                           #[:,508]
    'net1s2',                           #[:,509]

    'net1t1',                           #[:,510]
    'net1t2',                           #[:,511]

    'net1u1',                           #[:,512]
    'net1u2',                           #[:,513]

    'net1v1',                           #[:,514]
    'net1v2',                           #[:,515]

    'net1w1',                           #[:,516]
    'net1w2',                           #[:,517]

    'net1x1',                           #[:,518]
    'net1x2',                           #[:,519]

    'net1y1',                           #[:,520]
    'net1y2',                           #[:,521]

    'net1z1',                           #[:,522]
    'net1z2',                           #[:,523]

    'net2a1',                           #[:,524]
    'net2a2',                           #[:,525]
    'net2a3',                           #[:,526]
    'net2a4',                           #[:,527]
    'net2a5',                           #[:,528]
    'net2a6',                           #[:,529]
    'net2a7',                           #[:,530]

    'net2b',                            #[:,531]
    'net2c'],

    'groups': None,
    
    'bounds': [[90, 100], #int1a          #[:,0] ATTENTION PAS DE VIRGULE APRES LE DERNIER !!
        [8, 100], #int1b           #[:,1]
        [10, 100], #int1c           #[:,2]
        [60, 100], #int1d           #[:,3]
       
       #### 
         ....
         ....
         ....

         ....
         ....
         ....

        [57, 100], #int1c           #[:,524]
        [34, 100], #int1d           #[:,525]
        [96, 100], #int2a1          #[:,526]
        [59, 100], #int2a2          #[:,527]
        [34, 100], #int2b1          #[:,528]
        [63, 100], #int2b2          #[:,529]
        [81, 100], #int2b3          #[:,530]
        [57, 100], #int1a           #[:,531]
        [98, 100]]} #int1b           #[:,532]
                    




# Generate samples
param_values = sample(problem, N=1000, num_levels=4,
                      optimal_trajectories=None)

# To use optimized trajectories (brute force method),
# give an integer value for optimal_trajectories

# Run the "model" -- this will happen offline for external models
Y = ET(param_values)

# Perform the sensitivity analysis using the model output
# Specify which column of the output file to analyze (zero-indexed)
Si = morris.analyze(problem, param_values, Y, conf_level=0.95,
                    print_to_console=True,
                    num_levels=4, num_resamples=100)
# Returns a dictionary with keys 'mu', 'mu_star', 'sigma', and 'mu_star_conf'
# e.g. Si['mu_star'] contains the mu* value for each parameter, in the
# same order as the parameter file




def _sort_Si(Si, key, sortby='mu_star'):
    return np.array([Si[key][x] for x in np.argsort(Si[sortby])])


def _sort_Si_by_index(Si, key, index):
    return np.array([Si[key][x] for x in index])



def horizontal_bar_plot2(ax, Si, param_dict, sortby='mu_star', unit=''):
    '''Updates a matplotlib axes instance with a horizontal bar plot

    of mu_star, with error bars representing mu_star_conf
    '''

    assert sortby in ['mu_star', 'mu_star_conf', 'sigma', 'mu']

    # Sort all the plotted elements by mu_star (or optionally another
    # metric)
    names_sorted = _sort_Si(Si, 'names', sortby)
    mu_star_sorted = _sort_Si(Si, 'mu_star', sortby)
    mu_star_conf_sorted = _sort_Si(Si, 'mu_star_conf', sortby)

    # Plot horizontal barchart
    y_pos = np.arange(len(mu_star_sorted))
    plot_names = names_sorted

    out = ax.barh(y_pos,
                  mu_star_sorted,
                  xerr=mu_star_conf_sorted,
                  align='center',
                  ecolor='black',
                  **param_dict)

    ax.set_yticks(y_pos)
    ax.set_aspect(aspect=0.3)
    ax.set_yticklabels(plot_names)
    ax.set_xlabel(r'$\mu^\star$' + unit)

    ax.set_ylim(min(y_pos)-1, max(y_pos)+1)

    return out



plt.rc('xtick', labelsize=5) 
plt.rc('ytick', labelsize=0.3) 
plt
fig, ax1 = plt.subplots(1)
horizontal_bar_plot2(ax1, Si, {}, sortby='mu_star', unit="Coeficient de sensibilité")
#horizontal_bar_plot(ax1, Si, {}, unit=r"tCO$_2$/year")
#covariance_plot(ax2, Si, {}, unit=r"tCO$_2$/year")
#plt.tick_params(axis='x', labelsize=1)
plt.savefig("graphapp.pdf") 
plt.show()


















