#Import libraries

from rdkit import Chem
import re
from rdkit.Chem.rdchem import *
from rdkit.Chem.rdmolops import *

from collections import Counter
from pyparsing import *
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import copy
from time import sleep
from itertools import chain
from pyparsing import *
import pandas as pd

import networkx as nx
import matplotlib.pyplot as plt
from rdkit.Chem import Descriptors

from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM

from subprocess import check_call
import os
import copy 



# Function for producing unique values for the list.
# Fastest way to uniqify a list in Python

def f5(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result
   
   
   
   
   
   
   
# Calculate molecular mass from formula

def get_mass(formula):

    parts = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    mass = 0

    for index in range(len(parts)):
        if parts[index].isnumeric():
            continue

        atom = Chem.Atom(parts[index])
        multiplier = int(parts[index + 1]) if len(parts) > index + 1 and parts[index + 1].isnumeric() else 1
        mass += atom.GetMass() * multiplier

    return mass

# Function to parse molecular formula

def parse_formula(formula):
    elements = [(1, "H"),(2, "He"),(3, "Li"),(4, "Be"),(5, "B"),(6, "C"),(7, "N"),(8, "O"),(9, "F"),(10, "Ne"),(11, "Na"),(12, "Mg"),(13, "Al"),(14, "Si"),(15, "P"),(16, "S"),(17, "Cl"),(18, "Ar"),(19, "K"),(20, "Ca"),(21, "Sc"),(22, "Ti"),(23, "V"),(24, "Cr"),(25, "Mn"),(26, "Fe"),(27, "Co"),(28, "Ni"),(29, "Cu"),(30, "Zn"),(31, "Ga"),(32, "Ge"),(33, "As"),(34, "Se"),(35, "Br"),(36, "Kr"),(37, "Rb"),(38, "Sr"),(39, "Y"),(40, "Zr"),(41, "Nb"),(42, "Mo"),(43, "Tc"),(44, "Ru"),(45, "Rh"),(46, "Pd"),(47, "Ag"),(48, "Cd"),(49, "In"),(50, "Sn"),(51, "Sb"),(52, "Te"),(53, "I"),(54, "Xe"),(55, "Cs"),(56, "Ba"),(57, "La"),(58, "Ce"),(59, "Pr"),(60, "Nd"),(61, "Pm"),(62, "Sm"),(63, "Eu"),(64, "Gd"),(65, "Tb"),(66, "Dy"),(67, "Ho"),(68, "Er"),(69, "Tm"),(70, "Yb"),(71, "Lu"),(72, "Hf"),(73, "Ta"),(74, "W"),(75, "Re"),(76, "Os"),(77, "Ir"),(78, "Pt"),(79, "Au"),(80, "Hg"),(81, "Tl"),(82, "Pb"),(83, "Bi"),(84, "Po"),(85, "At"),(86, "Rn"),(87, "Fr"),(88, "Ra"),(89, "Ac"),(90, "Th"),(91, "Pa"),(92, "U"),(93, "Np"),(94, "Pu"),(95, "Am"),(96, "Cm"),(97, "Bk"),(98, "Cf"),(99, "Es"),(100, "Fm"),(101, "Md"),(102, "No"),(103, "Lr"),(104, "Rf"),(105, "Db"),(106, "Sg"),(107, "Bh"),(108, "Hs"),(109, "Mt"),(110, "Ds"),(111, "Rg"),(112, "Cn"),(113, "Uut"),(114, "Uuq"),(115, "Uup"),(116, "Uuh"),(118, "Uuo")]
    el_list = sorted(list(zip(*elements))[1],key=lambda x:len(x),reverse=True)
    el = Or([Word(x) for x in el_list])
    isotop = el ^ Group(Suppress("[") + OneOrMore(Word(nums)) + el + Suppress("]")).setParseAction(lambda str, location, tokens: "".join(tokens[0]))
    parser = OneOrMore(Group(isotop + ZeroOrMore(Word(nums))))
    parsed = parser.parseString(formula)
    res = set()
    for p in parsed:
        if p[0] =='H': continue
        if len(p) == 2:
            res.add((p[0],int(p[1])))
        else:
            res.add((p[0],1))
    return res
    
    
# Calculate fragments
# Main function. Enumerates all paths on molecular graph and select only those with correct brutto-formula

def calc_fragments(mol,target_with_d,check_one_molecule_rule=False, parent_indexes_condition=False):
    rdDepictor.Compute2DCoords(mol)
    target_mol = copy.deepcopy(mol)
    
    # Remove deuterium
    target=set()
    for idx in target_with_d:
        if idx[0]!='2H':
            target.add(idx)
    #????????????????
    target=target_with_d
    
    N = sum([v for (_,v) in target])
    M = N -1
    
    # It works for bonds! Not for atoms!!
    frags = rdkit.Chem.rdmolops.FindAllSubgraphsOfLengthMToN(target_mol,N-2,N+3, True) #it was M, N+1
    frags = chain.from_iterable(frags)
    
    # Finding atoms connected to found bonds
    indexes = set()
    indexes_list=[]
    All_Atoms_list=set()
    All_Atoms_list_list=[]
    for fragment in frags:
        bonds = [mol.GetBondWithIdx(idx) for idx in fragment]
        idxs = set()
        atoms = []
        Atoms_list=set()

        for bond in bonds:
            begin = bond.GetBeginAtom()
            end   = bond.GetEndAtom()

            if begin.GetIdx() not in idxs:
                idxs.add(begin.GetIdx())
                Atoms_list.add(begin)
                if begin.GetIsotope() != 0:
                    atoms.append(str(begin.GetIsotope()) + begin.GetSymbol())
                else:    
                    atoms.append(begin.GetSymbol())                    

            if end.GetIdx() not in idxs:
                idxs.add(end.GetIdx())
                Atoms_list.add(end)
                if end.GetIsotope() != 0:
                    atoms.append(str(end.GetIsotope()) + end.GetSymbol())     
                else:    
                    atoms.append(end.GetSymbol())
          
        atom_set = frozenset(Counter(atoms).items()) #sorted(set(Counter(atoms).items()),key=lambda x: x[0])        
       
        if (parent_indexes_condition==True):
            #parent_indexes=set(int(idx.GetProp('pi')) for idx in Atoms_list)
            parent_indexes=set(int(target_mol.GetAtomWithIdx(idx).GetProp("pi")) for idx in idxs)
            
        
        if (atom_set == target):

            indexes.add(frozenset(idxs))
            indexes_list.append(frozenset(idxs))
            if (parent_indexes_condition==True):
                All_Atoms_list.add(frozenset(parent_indexes))
                All_Atoms_list_list.append(frozenset(parent_indexes))
                #print(idxs)
                #print(parent_indexes)
            
    if (parent_indexes_condition==True):        
        return f5(All_Atoms_list_list)
    else:
        return f5(indexes_list)
       
       
       
       
# Calculate the number of fragments remaining in the molecule after removing selected atoms

def validate_connectivity(mol,idxs, parent_indexes_condition=False):
    edit_mol = Chem.RWMol(mol)
    if (parent_indexes_condition==False):
        for idx in sorted(list(idxs), reverse=True):
            edit_mol.RemoveAtom(idx)    
        return len(GetMolFrags(edit_mol)) #GetMolFrags(edit_mol) #
    if (parent_indexes_condition==True):
        
        mol_pi_list=[int(edit_mol.GetAtomWithIdx(idx).GetProp("pi")) for idx in range(len(edit_mol.GetAtoms()))]
        fragments_local_i=[mol_pi_list.index(idx) for idx in frag]
        for idx1 in sorted(list(fragments_local_i),reverse=True):
            edit_mol.RemoveAtom(idx1)
        return len(GetMolFrags(edit_mol))
        
        
        
        
        
        
# Function to get molecular weigth of the selected fragment of the molecule. M/z can be calculated using RDKit functions
# or taken from input Fragment_data

def get_mw_elem_of_fr(Parent_mol,highlight, TakeMzFromInput=False,Fragment_data=pd.DataFrame({"formula":["C9H9O4"], "mz":[181.0495],"formulaNoH":["C9O4"]})):
    edit_mol = Chem.RWMol(Parent_mol)
    Remaining_mol=set([i for i in range(len(edit_mol.GetAtoms()))])-highlight
    for idx in sorted(list(Remaining_mol), reverse=True):
        edit_mol.RemoveAtom(idx) 
    elem_formula=rdMolDescriptors.CalcMolFormula(edit_mol,True)
    mol_weight=round(Descriptors.ExactMolWt(edit_mol),5)
    
    if TakeMzFromInput==True:
        elem_formula_NoH=re.sub("H\d*", "", elem_formula)
        mol_weight_1=Fragment_data[Fragment_data["formulaNoH"]==elem_formula_NoH]["mz"]
        if (mol_weight_1.empty==False):
            mol_weight=float(mol_weight_1.iat[0])
        else:
            mol_weight=0
        
    return ([mol_weight,elem_formula]) #GetMolFrags(edit_mol) #
    
    
    
    
    
    
# Creates an SVG file to be drawn. 

def Gererate_SVG (Parent_mol,legend,highlight,Fragment_data):
    size_x=450 
    size_y=150
    
    elem_formula=get_mw_elem_of_fr(Parent_mol,highlight)[1]
    mol_weight=get_mw_elem_of_fr(Parent_mol,highlight,True,Fragment_data)[0]
    
    # Draw
    try:
        mc_mol = rdMolDraw2D.PrepareMolForDrawing(Parent_mol, kekulize=True)
    except ValueError: # <- can happen on a kekulization failure
        mc_mol = rdMolDraw2D.PrepareMolForDrawing(Parent_mol, kekulize=False) #edit_mol

    drawer = rdMolDraw2D.MolDraw2DSVG(size_x, size_y)     
    drawer.DrawMolecule(mc_mol,legend=legend+"; m="+str(mol_weight)+"; "+elem_formula,highlightAtoms=highlight)


    #drawer.DrawMolecule(mc_mol,legend=str(fragment_index))
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return(svg)
    
    
    
    
    
###########################################################################
# Check the correctness of the input data
#
###########################################################################
def Check_input_data (Fragment_data):

    for index in range(len(Fragment_data)):
        mass=round(get_mass(Fragment_data.iloc[index]['formula']),3)
        mass_ini=round(Fragment_data.iloc[index]['mz'],3)
        Check_fragment=""
        if abs(mass_ini-mass)>1:
            Check_fragment="!!!!!!!!!!!CHECK FRAGMENT!!!!!!!!!!!!"

        print(Fragment_data.iloc[index]['formula'],mass_ini,mass,mass_ini-mass,Check_fragment)




# Generate all possible fragments of molecule for given input and sorting according to the filtration rule.

def Generate_all_fragments (mol,Fragment_data):


    All_Fragments=[]
    For_order=pd.DataFrame()
    fragment_formula_index=0
    fragment_index=0

    for fragment_formula in Fragment_data["formula"]:
        target = parse_formula(fragment_formula)#{('C',11),('O',1)}
        fragments = calc_fragments(mol,target,False)
        All_Fragments=All_Fragments + (list(fragments))
        Mz=Mz_list[fragment_formula_index]
        fragment_formula_index=fragment_formula_index+1

        for fragment in fragments:

            Nummber_of_remaining_fragments = validate_connectivity(mol,fragment)     

            Bonds_between_fr_and_restofmol=0
            Remaining_mol=set([i for i in range(len(mol.GetAtoms()))])-fragment
            for fr_atom in fragment:
                for mol_atom in Remaining_mol:
                    bond=mol.GetBondBetweenAtoms(fr_atom,mol_atom)
                    if bond !=None:
                        Bonds_between_fr_and_restofmol=Bonds_between_fr_and_restofmol+1
                    #print(fr_atom,mol_atom,bond)

            For_order=For_order.append(pd.DataFrame({'nr':[Nummber_of_remaining_fragments],'fi':[fragment_index],'elem':[fragment_formula],'mz':[Mz],'bb':Bonds_between_fr_and_restofmol}),ignore_index=True)
            fragment_index=fragment_index+1

    For_order=For_order.sort_values(by=['nr','bb'])

    # Restoring row indexing to normal. Without it it will remain the same even after sorting!!!!
    For_order = For_order.reset_index(drop=True)
    Ordered_fragments=list( All_Fragments[i] for i in For_order['fi'] )
    
    return ([For_order,Ordered_fragments])



# Draw fragments

def Draw_fragments(mol, For_order, Ordered_fragments):

    fragment_index=0
    for fragment in Ordered_fragments:

        AllChem.Compute2DCoords(mol)
        try:
            mol.GetAtomWithIdx(0).GetExplicitValence()
        except RuntimeError:
            mol.UpdatePropertyCache(False)
        try:
            mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
        except ValueError: # <- can happen on a kekulization failure
            mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        drawer = rdMolDraw2D.MolDraw2DSVG(size_x, size_y)

        opts = drawer.drawOptions()

        for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetIsotope()!=0:
                opts.atomLabels[i] = str(mol.GetAtomWithIdx(i).GetIsotope())+mol.GetAtomWithIdx(i).GetSymbol()+str(i)
            else:
                opts.atomLabels[i] = mol.GetAtomWithIdx(i).GetSymbol()+str(i)



        drawer.DrawMolecule(mc_mol,highlightAtoms=fragment,legend=str(fragment_index)+'; m/z='+str(For_order['mz'][fragment_index])+"; nr="+str(For_order['nr'][fragment_index]) + "; b="+str(For_order['bb'][fragment_index])+ "                                          ")

        #drawer.DrawMolecule(mc_mol,legend=str(fragment_index))
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        display(SVG(svg))
        fragment_index=fragment_index+1
    
    
    
def Calculate_fragmentation_tree (Parent_mol,Number_of_levels,Fragment_data, Breaking_bonds_max, Remaining_fragments_max):

    # Define fragmentation graph
    Fragmentation_graph={Parent_mol:[]}

    # Store atom indexes in fragmentation graph. Both keys and values
    Parent_mol_set=frozenset(int(i.GetProp("pi")) for i in Parent_mol.GetAtoms())

    ###############################################################
    # Comment here to run loops!!!!!!!

    Fragmentation_graph_sets={Parent_mol_set: []}

    Fragmentation_graph_sets_edge={Parent_mol_set: []}

    ###############################################################

    level=0
    while level <= Number_of_levels:
        level=level+1
        # Get all parent fragments
        Parent_fragments=list(Fragmentation_graph.keys())+list(Fragmentation_graph.values())[0]
        #Parent_fragments_set=list(Fragmentation_graph_sets.values())[0]
        #Parent_fragments_set.add(list(Fragmentation_graph_sets.keys())[0])

        Parent_fragments_set={(list(Fragmentation_graph_sets.keys())[0])}
        Parent_fragments_set=Parent_fragments_set.union(set(chain.from_iterable(Fragmentation_graph_sets.values())))

        #for mol in Parent_fragments:

        for mol_pi in Parent_fragments_set:
            ###############
            Remaining_mol=set([i for i in range(len(Parent_mol.GetAtoms()))])-mol_pi
            ################
            mol = Chem.RWMol(Parent_mol)
            for idx in sorted(list(Remaining_mol), reverse=True):
                mol.RemoveAtom(idx)

            Fragment_structure_list=[]
            Fragment_structure_list_pi=[]
            Fragment_structure_list_pi_edge=[]

            fragment_formula_index=0

            for fragment_formula in Fragment_data["formula"]:
                target = parse_formula(fragment_formula)#{('C',11),('O',1)}
                fragments = calc_fragments(mol,target,False)
                fragments_pi = calc_fragments(mol,target,False,True)

                if fragments !=[]: 
                    Mz=Mz_list[fragment_formula_index]


                    For_order=pd.DataFrame()
                    fragment_index=0


                    for fragment in fragments:

                        Nummber_of_remaining_fragments = validate_connectivity(mol,fragment)     

                        Bonds_between_fr_and_restofmol=0
                        Remaining_mol=set([i for i in range(len(mol.GetAtoms()))])-fragment
                        for fr_atom in fragment:
                            for mol_atom in Remaining_mol:
                                bond=mol.GetBondBetweenAtoms(fr_atom,mol_atom)
                                if bond !=None:
                                    Bonds_between_fr_and_restofmol=Bonds_between_fr_and_restofmol+1
                            #print(fr_atom,mol_atom,bond)

                        For_order=For_order.append(pd.DataFrame({'nr':[Nummber_of_remaining_fragments],'fi':[fragment_index],'elem':[fragment_formula],'mz':[Mz],'bb':Bonds_between_fr_and_restofmol}),ignore_index=True)
                        fragment_index=fragment_index+1

                    For_order=For_order.sort_values(by=['nr','bb'])

                    # Restoring row indexing to normal. Without it it will remain the same even after sorting!!!!
                    For_order = For_order.reset_index(drop=True)
                    Ordered_fragments=list( list(fragments)[i] for i in For_order['fi'] )

                    # Condition to filter real fragments
                    # If working with deuterium this should be incresased by to 4 and 3 respectively!!!!
                    # Without deuterium I used 3 and 2
                    #
                    #Possible_Candidates=For_order[(For_order['bb']<3) & (For_order['nr']<2)]   #4 and 3            
                    Possible_Candidates=For_order[(For_order['nr']==For_order.iloc[0]['nr']) & (For_order['bb']==For_order.iloc[0]['bb'])]
                    Possible_Candidates=Possible_Candidates[(Possible_Candidates['bb']<Breaking_bonds_max) & (Possible_Candidates['nr']<Remaining_fragments_max)]
                    #
                    ##################################

                    Possible_fragments=list( list(fragments)[i] for i in Possible_Candidates['fi'] )
                    Possible_fragments_pi=list( list(fragments_pi)[i] for i in Possible_Candidates['fi'] )
                    Possible_fragments_pi_edge=[[int(For_order[For_order["fi"]==i]["nr"]),int(For_order[For_order["fi"]==i]["bb"])] for i in For_order["fi"]]


                    if (Possible_fragments!=[]):
                        print("Parent fragment "+str(mol_pi))
                        print("Generated fragment. Current indexes" + str(Possible_fragments))

                        print("Generated fragment. Parrent indexes" + str(Possible_fragments_pi))

                    for selected_fragment_index in range(len(Possible_fragments)):
                        selected_fragment = Possible_fragments[selected_fragment_index] 

                        Remaining_mol=set([i for i in range(len(mol.GetAtoms()))])-selected_fragment
                        edit_mol = Chem.RWMol(mol)
                        for idx in sorted(list(Remaining_mol), reverse=True):
                            edit_mol.RemoveAtom(idx) 

                        Fragment_structure_list.append(edit_mol)
                        Fragment_structure_list_pi.append(Possible_fragments_pi[selected_fragment_index])
                        Fragment_structure_list_pi_edge.append(Possible_fragments_pi_edge[selected_fragment_index])

                        # Draw
                        try:
                            mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
                        except ValueError: # <- can happen on a kekulization failure
                            mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False) #edit_mol

                        drawer = rdMolDraw2D.MolDraw2DSVG(size_x, size_y)     
                        drawer.DrawMolecule(mc_mol,legend=str(Mz),highlightAtoms=selected_fragment)


                        #drawer.DrawMolecule(mc_mol,legend=str(fragment_index))
                        drawer.FinishDrawing()
                        svg = drawer.GetDrawingText()
                        display(SVG(svg))

                    Fragmentation_graph[mol]=Fragment_structure_list                  #
                    Fragmentation_graph_sets[mol_pi]=Fragment_structure_list_pi       #
                    Fragmentation_graph_sets_edge[mol_pi]=Fragment_structure_list_pi_edge

                    print(Fragmentation_graph_sets)
                fragment_formula_index=fragment_formula_index+1
    
    for parent in Fragmentation_graph_sets.copy():
        for frag in Fragmentation_graph_sets[parent]:
            print(parent,frag)
            if (frag==parent):
                Fragmentation_graph_sets[parent].remove(frag)
                
    return(Fragmentation_graph_sets,Fragmentation_graph_sets_edge)


################################################################################################
# Draw generated fragmentatoin graph
#
#
#
################################################################################################

def Draw_tree (Parent_mol,Fragment_data,Fragmentation_graph_sets,Fragmentation_graph_sets_edge,path_to_pictures):

    parent_index=0
    fragment_index=0

    Unique_Fragments=set()
    All_Nodes_description=[]
    All_Edges_description=[]

    Fragments_Id=pd.DataFrame()
    fragment_id=0

    All_Nodes_description_Id=[]
    All_Edges_description_Id=[]

    for parent in Fragmentation_graph_sets:
        print("Displaying parent structure " + str(parent))
        parent_Colors={a: (0,1,0) for a in parent}

        Node_Name_P=list(parent)
        Node_Name_P=' '.join(map(str,Node_Name_P))
        Node_Name_P = Node_Name_P.replace(" ","a")

        if parent not in Unique_Fragments:
            print("dddddddddddddddd")
            Unique_Fragments.add(parent)
            elem_formula=get_mw_elem_of_fr(Parent_mol,parent)[1]
            mol_weight=get_mw_elem_of_fr(Parent_mol,parent,True,Fragment_data)[0]
            if (mol_weight==0):
                print("!!!!!! mol_weight=0")


            Fragments_Id=Fragments_Id.append(pd.DataFrame({'name':[Node_Name_P],'fi':[str(fragment_id)],"formula":[elem_formula],"mz":[mol_weight],"formulaNoH":[re.sub("H\d*", "", elem_formula)]}),ignore_index=True)
            fragment_id=fragment_id+1

            Node_description=Node_Name_P+" [label=\"\",image=\""+path_to_pictures+"Fig\\"+Node_Name_P+".png"+"\",labelloc=b];"
            All_Nodes_description.append(Node_description)

            node_id_P=Fragments_Id[Fragments_Id['name']==Node_Name_P].iat[0,1]
            Node_description_id=node_id_P+" [label=\"\",image=\""+path_to_pictures+"Fig\\"+Node_Name_P+".png"+"\",labelloc=b];"

            All_Nodes_description_Id.append(Node_description_id)

            svg=Gererate_SVG(Parent_mol,str(node_id_P),parent,Fragment_data)
            #write svg picture to file
            file_name=path_to_pictures+"Fig\\"+Node_Name_P
            with open (file_name+".svg",'w') as f:
                f.write(svg)

            drawing = svg2rlg(file_name +".svg")
            renderPM.drawToFile(drawing, file_name+".png",fmt="PNG")

        else:
            node_id_P=Fragments_Id[Fragments_Id['name']==Node_Name_P].iat[0,1]

        fragment_index=0

        svg=Gererate_SVG(Parent_mol,str(node_id_P),parent,Fragment_data)
        display(SVG(svg))

        for frag in Fragmentation_graph_sets[parent]:
            print("Displaying fragments " + str(frag))

            Node_Name_F=list(frag)
            Node_Name_F=' '.join(map(str,Node_Name_F))
            Node_Name_F = Node_Name_F.replace(" ","a")

            if frag not in Unique_Fragments:
                print("dddddddddddddddd")
                Unique_Fragments.add(frag)
                elem_formula=get_mw_elem_of_fr(Parent_mol,frag)[1]
                mol_weight=get_mw_elem_of_fr(Parent_mol,frag,True,Fragment_data)[0]


                Fragments_Id=Fragments_Id.append(pd.DataFrame({'name':[Node_Name_F],'fi':[str(fragment_id)],"formula":[elem_formula],"mz":[mol_weight],"formulaNoH":[re.sub("H\d*", "", elem_formula)]}),ignore_index=True)
                fragment_id=fragment_id+1

                Node_description=Node_Name_F+" [label=\"\",image=\""+path_to_pictures+"Fig\\"+Node_Name_F+".png"+"\",labelloc=b];"
                All_Nodes_description.append(Node_description)

                node_id_F=Fragments_Id[Fragments_Id['name']==Node_Name_F].iat[0,1]
                Node_description_id=node_id_F+" [label=\"\",image=\""+path_to_pictures+"Fig\\"+Node_Name_F+".png"+"\",labelloc=b];"

                All_Nodes_description_Id.append(Node_description_id)

                svg=Gererate_SVG(Parent_mol,str(node_id_F),frag,Fragment_data)
                #write svg picture to file
                file_name=path_to_pictures+"Fig\\"+Node_Name_F
                with open (file_name+".svg",'w') as f:
                    f.write(svg)

                drawing = svg2rlg(file_name +".svg")
                renderPM.drawToFile(drawing, file_name+".png",fmt="PNG")

            else:
                node_id_F=Fragments_Id[Fragments_Id['name']==Node_Name_F].iat[0,1]

            # Draw
            svg=Gererate_SVG(Parent_mol,str(node_id_F),frag,Fragment_data)
            display(SVG(svg))

            Edge_description=Node_Name_P+" -- "+Node_Name_F+";"
            All_Edges_description.append(Edge_description)

            Edge_description_id=node_id_P+" -- "+node_id_F+";"
            All_Edges_description_Id.append(Edge_description_id)

            fragment_index=fragment_index+1
        parent_index=parent_index+1

    #####
    edges_pairs=[]
    for i in list(Fragmentation_graph_sets_edge.values()):
        edges_pairs=edges_pairs+i

    Ind1=0
    Ind2=0
    k=0
    for edge in edges_pairs:
        if (edge!=[0,0]):
            k=k+1
            Ind1=Ind1+edge[0]
            Ind2=Ind2+edge[1]

    print (Ind1/k,Ind2/k)

    UsedFragments=len(f5(Fragments_Id['formulaNoH']))/(len(Fragment_data['formulaNoH']))
    print(UsedFragments)
    # Write grap description to file

    f=open(path_to_pictures+'For_Graphviz.txt','w')
    f.write("graph G {"+'\n')

    f.write("label = <"+"Remaining fragments Index = "+str(Ind1/k)+ "; Breaking bond index = " +str(round(Ind2/k,3))+ "; Used fragments (%) index = " +str(round(UsedFragments,3))+ ">;"+'\n')

    f.write("labelloc = \"t\";"+'\n')

    ####################################
    # Add table with ions
    f.write("aHtmlTable [shape=plaintext;  label=< <table border='1' cellborder='0'> \n")
    f.write("<tr>")
    for formulaNoH in Fragment_data['formulaNoH']:
        f.write("<td>"+formulaNoH+"</td>\n")
    f.write("<td>"+str(len(Fragment_data['formulaNoH']))+"</td>\n")
    f.write("</tr>\n")
    f.write("<tr>")
    for formulaNoH in Fragment_data['formulaNoH']:
        if formulaNoH in list(Fragments_Id['formulaNoH']):
            f.write("<td>"+"+"+"</td>\n")
        else:
            f.write("<td>"+"-"+"</td>\n")
    f.write("<td>"+str(len(f5(Fragments_Id['formulaNoH'])))+"</td>\n")

    f.write("</tr>\n")

    f.write("</table> >]; \n")

    for ele in All_Nodes_description_Id:
        f.write(ele+'\n')
    f.write('\n')
    for ele in All_Edges_description_Id:
        f.write(ele+'\n')

    f.write("}"+'\n')
    f.close()

    with open(path_to_pictures+'Fragmentation_graph_sets.txt', 'w') as f:
        print(Fragmentation_graph_sets, file=f)

    with open(path_to_pictures+'Fragmentation_graph_sets_edge.txt', 'w') as f:
        print(Fragmentation_graph_sets_edge, file=f)

    Fragments_Id.to_csv(path_to_pictures+'Fragments_Id.txt', sep='\t')
    
    
    
