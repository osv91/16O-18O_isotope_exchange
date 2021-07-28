from rdkit import Chem
from rdkit.Chem.inchi import *
import itertools
def identify_O(mol, substructs):
    exch_O = set()
    
    if len(substructs)>0:
        for i in substructs:
            if len(i)>0:
                
                for j in i:
                    
                    if mol.GetAtomWithIdx(j).GetSymbol()=='O':
                        exch_O.add(j)
                        
    
    return exch_O

def calculate_N_O18(inchi):
    mol = MolFromInchi(inchi)
    exchangeable_O = set()
    
    
    core_keto ='[#6]([#6])(=O)[#6]'
    core_keto_mol = Chem.MolFromSmarts(core_keto)
    exchangeable_O = exchangeable_O.union(identify_O(mol, mol.GetSubstructMatches(core_keto_mol)))
    
    core_ald = '[CX3H1](=O)[#6]'
    core_ald_mol  = Chem.MolFromSmarts(core_ald)
    exchangeable_O =exchangeable_O.union(identify_O(mol, mol.GetSubstructMatches(core_ald_mol)))
    
    core_carboxyl= '[CX3](=O)[OX2H1]'
    core_carboxyl_mol  = Chem.MolFromSmarts(core_carboxyl)
    exchangeable_O = exchangeable_O.union(identify_O(mol, mol.GetSubstructMatches(core_carboxyl_mol)))
    
    core_enol_oh = '[#6]=[#6][#6]([!R;OH])[#6]'
    core_enol_oh_mol = Chem.MolFromSmarts(core_enol_oh)
    exchangeable_O =exchangeable_O.union(identify_O(mol, mol.GetSubstructMatches(core_enol_oh_mol)))
    
    core_aryl_oh = '[Rc]C([!O])[OH]'
    core_aryl_oh_mol = Chem.MolFromSmarts(core_aryl_oh)
    exchangeable_O =exchangeable_O.union(identify_O(mol, mol.GetSubstructMatches(core_aryl_oh_mol)))
       
    return len(exchangeable_O)


def identify_exchangeable_O(inchi):
    l = set()
    mol = MolFromInchi(inchi)
    
    core_keto ='[#6]([#6])(=O)[#6]'
    core_keto_mol = Chem.MolFromSmarts(core_keto)
    
    l.add(mol.GetSubstructMatches(core_keto_mol))
    
    core_ald = '[CX3H1](=O)[#6]'
    core_ald_mol  = Chem.MolFromSmarts(core_ald)
    l.add(mol.GetSubstructMatches(core_ald_mol))
    
    core_carboxyl= '[CX3](=O)[OX2H1]'
    core_carboxyl_mol  = Chem.MolFromSmarts(core_carboxyl)
    l.add(mol.GetSubstructMatches(core_carboxyl_mol))
    
    core_enol_oh = '[#6]=[#6][#6]([OX2H])[#6]'
    core_enol_oh_mol = Chem.MolFromSmarts(core_enol_oh)
    l.add(mol.GetSubstructMatches(core_enol_oh_mol))
    
    core_aryl_oh = '[Rc]C([!O])[OH]'
    core_aryl_oh_mol = Chem.MolFromSmarts(core_aryl_oh)
    l.add(mol.GetSubstructMatches(core_aryl_oh_mol))
    s=set()
    
    for i in l:
        for j in i:
            for z in j:
             
                if mol.GetAtomWithIdx(z).GetSymbol()=='O':
                    s.add(z)
      
       
    return s


def run_16O18O_exchange (inchi, N_label):
    mol = MolFromInchi(inchi)
    labeled_mols = []
    s= identify_exchangeable_O(inchi)
    for i in itertools.combinations(s, N_label):
        molO18 = MolFromInchi(inchi)
        for j in range(N_label):
            
            molO18.GetAtomWithIdx(i[j]).SetIsotope(18)
        labeled_mols.append(molO18)
    return (labeled_mols)    
    
def clean_O18_MS1(InChIs, N, formula):
    InChIs['N'] = InChIs['InChI'].apply(calculate_N_O18)
    InChIs = InChIs[InChIs.N>=N]
    with open ('./results/'+formula+'.txt', 'a') as f:
        f.write('Number of candidates match O18: ' + str(len(InChIs))+'\n')
        
    return InChIs
    
    
