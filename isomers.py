from rdkit import Chem
from rdkit.Chem.inchi import *
import pandas as pd
import pubchempy as pcp


def get_formula_from_inchi(inchi):
    formula = inchi.split('/')[1]
    return formula

def find_salts(smiles):
    return smiles.find('.')==-1   
    
def get_isomers_for_formula(formula,output_path = 'temp.csv'):
    #get all isomers from PubChem for formula
    cmpds=pcp.get_cids(formula, 'formula')
    pcp.download('CSV', output_path, cmpds, operation='property/InChi', overwrite=True)
    data=pd.read_csv(output_path, sep=',')
    l1 = len(data)
    
    #Remove bad mols
    data['mol']=data['InChI'].apply(MolFromInchi)
    data = data.dropna(subset=['mol'])
    #Remove Hydrogens
    data['mol']=data['mol'].apply(Chem.RemoveHs)
    #Remove isomers and salts
    data['smiles']=data['mol'].apply(Chem.MolToSmiles, isomericSmiles=False)
    data['OK']=data['smiles'].apply(find_salts)
    data=data[data['OK']==True]
    #Drop duplicates
    data.drop_duplicates(subset ="InChI", inplace = True)
    data.drop_duplicates(subset ="smiles", inplace = True)
    data =data.reset_index()
    l2 = len(data)
    
    with open ('./results/'+formula+'.txt', 'w') as f:
        f.write('Formula: ' + formula+'\n')
        f.write('Total number of retrieved molecules: '+str(l1)+'\n')
        f.write('Number of candidates after cleaning: '+ str( l2)+'\n')


    return data[['InChI']]
    
    
