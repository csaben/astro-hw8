"""
N2/Nbary = ( B/(1+B) ) * (1-X)

s.t.
X(T) = Np/Nbary
B(T) = N2/N1

given that the first eqn ratio is proportional to the relative strength of balmer lines

plot the rel Strength vs T
"""
import math
import numpy as np
import matplotlib.pyplot as plt

np.seterr('raise')

def main():
    #variables
    global Nbary; Nbary = 10e20 #m-3
    global chi; chi = -13.6 #eV
    # global k; k  = 1#1.38e-23 #m2 kg s-2 K-1  #boltzmann factor always causes an error
    global k; k  = 8.e-5#eV K-1
    global m; m = 9.109e-31 #kg
    global h; h = 6.626e-34 #m2 kg s-1
    global eV; eV = 1.602e-19 #convert to joules
    #when X is small and B is small str is 0, when X is big and B is big str is 0
    T = np.arange(1e4,1e6)
    strengths = strength(T)
    plot_fn(T, strengths)
    # max = np.max(strengths)

    #i expect Tmax ~ 10000K from class

def pseudo(T):
    # use algebra to exclude the B term
    return 1 / (Nbary*(1-X(T)))

def B(T):
    try:
        return 4*np.exp(-10.2/(k*T))
    except Exception as e:
        print(e, k*T)

def F(T):
    return Nbary * ( np.power( ( ( 2*np.pi*m*k*T ) / (np.square(h)) ) ,(3/2) ) ) * (np.exp(chi/(k*T)))
    # return Nbary * ( np.power(  ( 2*np.pi*m*k*T ) ,2) )
    # return 1

def X(T):
    return ( 1/(2*F(T))) * (-1+np.sqrt(1+4*F(T)))
    # return 1/2

def strength(T):
    return ( B(T) / (1+B(T)) ) * (1 - X(T))

def plot_fn(T, relStrength) :
    plt.plot(T, relStrength, linestyle= "dashed", color = "blue", label = "balmer lines") 
    plt.xlabel('T')
    plt.ylabel('relative strength')
    plt.legend()
    #plt.suptitle("problem 1: graph rel. strength of the balmer lines for nbary=10^20m^-3 as func of T and read off the temperature @ which you expect the strongest balmer lines")
    plt.suptitle("@ T=10000K")
    plt.title("hw 8 ; 9.1")
    plt.savefig(f"hw8p91", dpi=150) #saves the plot
    plt.show()


if __name__ == '__main__':
    main()
