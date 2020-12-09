# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:02:12 2020

@author: Pran Regu
"""

# Ring numbering

from rdkit import Chem
from rdkit.Chem import Draw

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

def rnumber(Compd):
    index = 1
    Comp = Compd
    C2 = ''
    for i in Comp:
        j = i
        if i == '(':
            index+=1
        elif i == ')':
            index-=1
        if is_number(i) == True:
            ind_str = str(index)
            j = ind_str
            Comp = Comp.replace(i, ind_str, 1)
        C2+= j    
    return C2   


if __name__ == "__main__":
    #Compound = 'C4=CC=C(O)C=C4(C5=CC(C6=C(O)C=CC=C6)=CC(O)=C5)'
    Compound = 'C4=CC=CC=C4(C5=CC(C6=CC=C(C6=CC=C(C6=CC=C(C2=CC(C#N)=C(C3=CC=CC=C3)C=C2)C=C6)C=C6)C=C6)=CC=C5)'
    Mod_Comp = rnumber(Compound)
    print('\nOriginal SMILES: ', Compound)
    print('\nModified SMILES: ', Mod_Comp)

    Draw.MolToMPL(Chem.MolFromSmiles(Compound), size = (200, 200))
    Draw.MolToMPL(Chem.MolFromSmiles(Mod_Comp), size = (200, 200))
   