import streamlit as st
from polymer_utils import build_polymer, visualize_molecule
from rdkit import Chem

st.title("Generador de Polímeros 🧪")

# Entrada del usuario
monomer_smiles = st.text_input("Introduce el SMILES del monómero:", "CC(=O)O")  # Por defecto, ácido láctico
num_units = st.number_input("Número de repeticiones:", min_value=1, max_value=100, value=10)
reaction_type = st.selectbox("Tipo de reacción:", ["heterolítica", "radicalaria", "policondensación"])

if st.button("Generar Polímero"):
    try:
        polymer_mol = build_polymer(monomer_smiles, num_units, reaction_type)
        st.success("¡Polímero generado con éxito!")

        # Mostrar el SMILES del polímero
        polymer_smiles = Chem.MolToSmiles(polymer_mol)
        st.text(f"SMILES del polímero: {polymer_smiles}")

        # Visualización en 3D
        st.subheader("Visualización en 3D:")
        viewer = visualize_molecule(polymer_mol)
        st.components.v1.html(viewer.show(), height=450)

    except Exception as e:
        st.error(f"Error en la generación del polímero: {e}")
