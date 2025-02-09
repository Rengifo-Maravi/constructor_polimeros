from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def build_polymer(monomer_smiles, num_units, reaction_type):
    """
    Construye un polímero según el tipo de polimerización.

    reaction_type:
        - "heterolítica" -> Une los monómeros sin modificar su estructura.
        - "radicalaria" -> Rompe enlaces dobles antes de unir los monómeros.
        - "policondensación" -> Elimina un grupo funcional de un extremo antes de unir los monómeros.
    """
    monomer = Chem.MolFromSmiles(monomer_smiles)  # Convertir SMILES a Mol
    monomer = Chem.RWMol(monomer)  # Convertimos a Mol editable
    
    if reaction_type == "radicalaria":
        # Busca y rompe los enlaces dobles para permitir la polimerización
        for bond in monomer.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:  
                bond.SetBondType(Chem.rdchem.BondType.SINGLE)  

    elif reaction_type == "policondensación":
        # Buscar grupos funcionales reactivos en los extremos
        reactive_groups = ["O", "N", "Cl"]
        atoms_to_remove = []
        
        for atom in monomer.GetAtoms():
            if atom.GetSymbol() in reactive_groups and len(atom.GetBonds()) == 1:
                atoms_to_remove.append(atom.GetIdx())
        
        if atoms_to_remove:
            monomer.RemoveAtom(atoms_to_remove[0])  # Eliminar solo un extremo

    # Convertimos nuevamente a SMILES para repetir el monómero y formar el polímero
    monomer_smiles_modified = Chem.MolToSmiles(monomer)
    polymer_smiles = monomer_smiles_modified * num_units  
    
    # Convertimos el polímero a RDKit Mol y optimizamos su geometría 3D
    polymer_mol = Chem.MolFromSmiles(polymer_smiles)
    polymer_mol = Chem.AddHs(polymer_mol)
    AllChem.EmbedMolecule(polymer_mol, AllChem.ETKDG())  
    AllChem.UFFOptimizeMolecule(polymer_mol)  
    
    return polymer_mol  

def visualize_molecule(mol):
    """
    Visualiza una molécula en 3D usando py3Dmol.
    """
    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer
