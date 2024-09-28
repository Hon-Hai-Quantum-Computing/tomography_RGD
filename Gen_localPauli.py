
import numpy as np

def counting_elements(Nk, sym_Pauli):
    """to count numbers for Pauli operators of different weights of 0, 1, 2, etc.

    Args:
        Nk (int): number of quibits
        sym_Pauli (list): = ['X', 'Y', 'Z']

    Returns:
        list: a list of numbers for Pauli operators of weights 0, 1, 2, etc.
    """

    nSym = len(sym_Pauli)

    n0 = 1                                  #  all identity
    n1 = Nk    * nSym                       # num of single site replaced by Pauli
    n2 = int(Nk*(Nk-1)/2)        * nSym**2  # num of two sited replaced
    n3 = int(Nk*(Nk-1)*(Nk-2)/6) * nSym**3  # num of 3 sites replaced
    n4 = int(Nk*(Nk-1)*(Nk-2)*(Nk-3)/24) * nSym**4

    print(' (n0,n1,n2,n3,n4) = ({},{},{},{},{})'.format(n0,n1,n2,n3,n4))

    NumTerm_Pauli = {0: 1}   # identity (0 site replaced): only one term 
    Factorial = 1
    Pick2e    = 1
    sumALL    = 1
    for nR in np.arange(1, Nk+1):
        Pick2e *= (Nk-nR+1)
        Factorial *= nR
        power  = nSym**nR
        numR   = int(Pick2e / Factorial) * power
        sumALL = sumALL + numR

        NumTerm_Pauli[nR] = numR
        print('Rep {}: (Pick2e, Factorial, power, numR, sum) = ({}, {}, {}, {}, {})'.format(
               nR, Pick2e, Factorial, power, numR, sumALL))

    if sumALL != 4**Nk:
        print(' ***  ERROR in counting all elments')
        return
    else:
        print('  ***   num of ALL elements = {}'.format(4**Nk))

        return NumTerm_Pauli


# ------------------------------------------------- #
#       deal with single Pauli operator changed     #
# ------------------------------------------------- #        
def SingleSiteReplace(base, symbol, site):

    front = base[:site]
    after = base[site+1:]

    SingleSiteList = [front+ss+after for ss in symbol]

    #print(front)
    #print(after)
    #print(SingleSiteList)
    return SingleSiteList
        
def Gen_All_singlePauli(Nk, symbol):    
    """ Generate all single site Pauli matrices  {X,Y,Z}
    """

    base = 'I'*Nk

    ALLsingle = []
    for site in np.arange(Nk):
        Plist = SingleSiteReplace(base, symbol, site)
        ALLsingle += Plist

    if len(ALLsingle) !=  len(symbol)*Nk:
        print(ALLsingle)
        print(' ***  ERROR:  the length NOT consistent')
        return
    else:
        print(' ***   number of All replaced single Pauli = {}'.format(len(ALLsingle)))

    return ALLsingle

# ------------------------------------------------- #
#       deal with two Pauli operators changed       #
# ------------------------------------------------- #        

def Replace2site_s1s2(base, symb2s, s1, s2):
    First_Site  = [xx[0] for xx in symb2s]

    ChangeFirst = SingleSiteReplace(base, First_Site, s1) 

    New2s = [xx[:s2]+zz[1]+xx[s2+1:] for xx, zz in zip(ChangeFirst, symb2s)]

    #Replace2site = []
    #for xx, zz in zip(ChangeFirst, symb2s):
    #    front = xx[:s2]
    #    backW = xx[s2+1:]
    #    tt = front + zz[1] + backW
    #    Replace2site.append(tt)

        #print(xx, zz, zz[1])
        #print(front, backW, tt)

    return New2s

def ReplaceIden2site(Nk, symb_P9):

    print(' -----------------------------------------------')
    base = 'I'*Nk
    ALLnew2S = []

    num2s = 0
    for s1 in np.arange(Nk-1):
        for s2 in np.arange(s1+1, Nk):
            NewS1S2 = Replace2site_s1s2(base, symb_P9, s1, s2)
            ALLnew2S += NewS1S2

            num2s += 1
            print('s1 = {}, s2= {}, New2s = {}'.format(s1, s2, NewS1S2))
    print(' num2s = {}'.format(num2s))

    if num2s*len(symb_P9) !=  len(ALLnew2S):
        print(' ###   number of replaced 2 indices = {}'.format(len(ALLnew2S)))
        print(' *** WRONG:  some repeated indices happen  ***')
        return
    else:
        print(' ###   number of ALLnew2S = {}'.format(len(ALLnew2S)))

    return ALLnew2S

# ------------------------------------------- #
#   Generate single site &  two site Pauli    #
# ------------------------------------------- #

def MultiplyList(sym1, sym2):

    ListMultiply = []
    for xx in sym1:
        yy = [xx + zz for zz in sym2]
        ListMultiply += yy
        #print(yy)
        #print(ListMultiply)

    return ListMultiply



def Gen_Site_2s_Pauli(Nk):

    symbol_P4  = ['I', 'X', 'Y', 'Z']
    sym_Pauli  = ['X', 'Y', 'Z']
    
    symb_16 = MultiplyList(symbol_P4, symbol_P4)
    
    symb_P9 = MultiplyList(sym_Pauli, sym_Pauli)

    ALL_2sPauli = ReplaceIden2site(Nk, symb_P9)

    ALL_sitePauli = Gen_All_singlePauli(Nk, sym_Pauli)

    counting_elements(Nk, sym_Pauli)

    numB2s = len(ALL_sitePauli) + len(ALL_2sPauli)
    print('  ***   num of single site &  two site Pauli = {}'.format(numB2s))

    return ALL_sitePauli, ALL_2sPauli, symb_P9, numB2s

