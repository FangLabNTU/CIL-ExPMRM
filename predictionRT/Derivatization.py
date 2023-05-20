import sys, os

import kora as kora
import pandas as pd


# install rdkit
#import pip
#pip install kora

# install molmass for molecular mass calculation
#!pip install molmass

import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from molmass import Formula
from rdkit.Chem import MolStandardize

def canonicalize_smiles(smiles):
  '''canonicalize compound smiles'''
  mol = Chem.MolFromSmiles(smiles)
  if mol is not None:
    return Chem.MolToSmiles(mol, isomericSmiles=True)
  else:
    return ''

def whether_COOH(smiles):
  '''determine whether compound contains carboxyl group'''
  mol = Chem.MolFromSmiles(smiles)
  COOH_patt = Chem.MolFromSmarts('[#6]C(=O)[O;H,-1]')
  COOH_info = mol.HasSubstructMatch(COOH_patt)
  return COOH_info

def whether_OH(smiles):
  '''determine whether compound contains hydroxyl group'''
  mol = Chem.MolFromSmiles(smiles)
  ArOH_patt = Chem.MolFromSmarts('c[OH]')
  ArOH_info = mol.HasSubstructMatch(ArOH_patt)
  AlOH_patt = Chem.MolFromSmarts('[C!$(C=O)]-[OH]')
  AlOH_info = mol.HasSubstructMatch(AlOH_patt)
  if ArOH_info or AlOH_info:
    return True
  else:
    return False

def MEPA_derivatization(smiles):
  '''generate MPEA product'''
  mol = Chem.MolFromSmiles(smiles)
  repl = Chem.MolFromSmiles('CN(C)CCC1=CC=CC=C1')
  patt = Chem.MolFromSmarts('[C$(C=O)]-[OH]')
  rms = AllChem.ReplaceSubstructs(mol, patt, repl)
  try:
    return Chem.MolToSmiles(rms[0])
  except:
    return ''

def DnsCl_derivatization(smiles):
  '''generate DnsCl product'''
  mol = Chem.MolFromSmiles(smiles)
  repl = Chem.MolFromSmiles('CN(C)c(ccc1)c2c1c(S(=O)(OC)=O)ccc2')
  if mol.HasSubstructMatch(Chem.MolFromSmarts('c[OH]')):
    rms = AllChem.ReplaceSubstructs(mol,Chem.MolFromSmarts('c[OH]'),repl,False,13)
    return Chem.MolToSmiles(rms[0])
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[C!$(C=O)]-[OH]')):
    rms = AllChem.ReplaceSubstructs(mol,Chem.MolFromSmarts('[C!$(C=O)]-[OH]'),repl,False,13)
    return Chem.MolToSmiles(rms[0])
  else:
    return ''

def calculate_mass(smiles):
  '''calculate molecular mass'''
  mol = Chem.MolFromSmiles(smiles)
  f = Formula(rdMolDescriptors.CalcMolFormula(mol))
  return f.isotope.mass

def derivatization1(compound_smiles):
  '''get derivatization product and its molecular mass'''
  if whether_OH(compound_smiles):
    product = DnsCl_derivatization(canonicalize_smiles(compound_smiles))
    mass = calculate_mass(product)
    return product, mass
  else:
    return '', ''

def derivatization2(compound_smiles):
  '''get derivatization product and its molecular mass'''
  if whether_COOH(compound_smiles):
    product = MEPA_derivatization(canonicalize_smiles(compound_smiles))
    mass = calculate_mass(product)
    return product, mass
  else:
    return '', ''

#derivatization("CCC(C)(C1=CC=C(O)C=C1)C1=CC=C(O)C=C1", "DnsCl")
#define function for derivatization and output derivatized SMILES and mass
#输入数据，表格，表格列参数[CARSN,name,SMILES]
def derivatization_operator(inputfile,outfilename):
  ###对每行的smiles 用转化函数derivatization1()
  #保存所得的两个函数输出结果
  print(inputfile)
  inputdata = pd.read_csv(inputfile)
  smilestr = inputdata['SMILES']
  DnsCl_smi_list = []
  DnsCl_mass_list = []
  MPEA_smi_list = []
  MPEA_mass_list = []
  for index in smilestr.index:
    x = smilestr.loc[index]
    x = standardize_smiles(x)
    DnsCl_smi, DnsCl_mass = derivatization1(x)
    MPEA_smi, MPEA_mass = derivatization2(x)
    DnsCl_smi_list.append(DnsCl_smi)
    DnsCl_mass_list.append(DnsCl_mass)
    MPEA_smi_list.append(MPEA_smi)
    MPEA_mass_list.append(MPEA_mass)

  print(len(DnsCl_smi_list))
  inputdata['DnsCl_SMILES'] = DnsCl_smi_list
  inputdata['DnsCl_mass'] = DnsCl_mass_list
  inputdata['MPEA_SMILES'] = MPEA_smi_list
  inputdata['MPEA_mass'] = MPEA_mass_list

  inputdata.to_csv(outfilename,index =False)
  return inputdata

def standardize_smiles(smiles):
    try:
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)

        # Standardize the molecule
        lfc = MolStandardize.fragment.LargestFragmentChooser()
        standard_mol = lfc.choose(mol)

        # Remove any remaining charges
        udc = MolStandardize.charge.Uncharger()
        standard_mol = udc.uncharge(standard_mol)

        # Convert standardized molecule back to SMILES
        standard_smiles = Chem.MolToSmiles(standard_mol, isomericSmiles=True)

        return standard_smiles
    except Exception as e:
        print(f"Error standardizing SMILES: {smiles}")
        print(e)
        return None


if __name__ == '__main__':
  derivatization_operator("MSdatabase_query_data.csv","test.csv")
