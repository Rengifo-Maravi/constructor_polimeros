import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from pysoftk.linear_polymer.linear_polymer import Lp
from pysoftk.topologies.ranpol import Rnp
from pysoftk.topologies.diblock import Db
from pysoftk.format_printers.format_mol import Fmt


# Función para obtener la molécula desde un SMILES
def obtener_molecula(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"SMILES no válido: {smiles}")
    return mol

# Función para construir polímeros lineales
def generar_polimero_lineal(smiles, group, n_units, shift, output_filename):
    try:
        monomero = obtener_molecula(smiles)
        AllChem.EmbedMolecule(monomero)
        polimero = Lp(monomero, group, n_units, shift).linear_polymer("MMFF", 350)
        Fmt(polimero).xyz_print(output_filename)
        return f"Polímero lineal guardado en: {output_filename}"
    except Exception as e:
        return f"Error: {str(e)}"

# Función para construir polímeros aleatorios
def generar_polimero_random(smiles_a, smiles_b, group, ratio, n_units, output_filename):
    try:
        mol_a = obtener_molecula(smiles_a)
        mol_b = obtener_molecula(smiles_b)
        random_polymer = Rnp(mol_a, mol_b, group).random_ab_copolymer(n_units, ratio, 350)
        Fmt(random_polymer).xyz_print(output_filename)
        return f"Polímero aleatorio guardado en: {output_filename}"
    except Exception as e:
        return f"Error: {str(e)}"

# Función para construir polímeros dibloques
def generar_polimero_dibloque(smiles_a, smiles_b, group, n_a, n_b, output_filename):
    try:
        mol_a = obtener_molecula(smiles_a)
        mol_b = obtener_molecula(smiles_b)
        diblock_polymer = Db(mol_a, mol_b, group).diblock_copolymer(n_a, n_b, "MMFF", 350)
        Fmt(diblock_polymer).xyz_print(output_filename)
        return f"Polímero dibloque guardado en: {output_filename}"
    except Exception as e:
        return f"Error: {str(e)}"

# Función para construir polímeros tribloques
def generar_polimero_tribloque(smiles_a, smiles_b, smiles_c, group, n_units, ratio_ab, ratio_c, output_filename):
    try:
        mol_a = obtener_molecula(smiles_a)
        mol_b = obtener_molecula(smiles_b)
        mol_c = obtener_molecula(smiles_c)
        triblock_polymer = Rnp(mol_a, mol_b, group).random_abc_copolymer(mol_c, n_units, ratio_ab, ratio_c, 350)
        Fmt(triblock_polymer).xyz_print(output_filename)
        return f"Polímero tribloque guardado en: {output_filename}"
    except Exception as e:
        return f"Error: {str(e)}"

# Interfaz de Streamlit
st.title("Constructor de Polímeros")
st.write("Construye polímeros lineales, dibloques, tribloques o aleatorios usando RDKit y PySoftK.")

tipo_polimero = st.selectbox("Selecciona el tipo de polímero", ["Lineal", "Aleatorio", "Dibloque", "Tribloque"])

# Formulario para polímero lineal
if tipo_polimero == "Lineal":
    smiles = st.text_input("SMILES del Monómero", value="c1cc(sc1Br)Br")
    group = st.text_input("Grupo de Conexión", value="Br")
    n_units = st.number_input("Número de Unidades", min_value=1, max_value=100, value=5)
    shift = st.number_input("Desplazamiento entre Unidades", min_value=0.5, max_value=5.0, value=1.25)
    output_filename = st.text_input("Nombre del Archivo de Salida (.xyz)", value="lineal.xyz")
    if st.button("Generar Polímero Lineal"):
        st.write(generar_polimero_lineal(smiles, group, n_units, shift, output_filename))

# Formulario para polímero aleatorio
elif tipo_polimero == "Aleatorio":
    smiles_a = st.text_input("SMILES del Monómero A", value="BrCOCBr")
    smiles_b = st.text_input("SMILES del Monómero B", value="c1(ccc(cc1)Br)Br")
    group = st.text_input("Grupo de Conexión", value="Br")
    ratio = st.slider("Proporción A/B", min_value=0.1, max_value=1.0, value=0.5, step=0.1)
    n_units = st.number_input("Número de Unidades", min_value=1, max_value=100, value=10)
    output_filename = st.text_input("Nombre del Archivo de Salida (.xyz)", value="random.xyz")
    if st.button("Generar Polímero Aleatorio"):
        st.write(generar_polimero_random(smiles_a, smiles_b, group, ratio, n_units, output_filename))

# Formulario para polímero dibloque
elif tipo_polimero == "Dibloque":
    smiles_a = st.text_input("SMILES del Monómero A", value="BrCOCBr")
    smiles_b = st.text_input("SMILES del Monómero B", value="c1cc(sc1Br)Br")
    group = st.text_input("Grupo de Conexión", value="Br")
    n_a = st.number_input("Número de Unidades en Bloque A", min_value=1, max_value=100, value=5)
    n_b = st.number_input("Número de Unidades en Bloque B", min_value=1, max_value=100, value=5)
    output_filename = st.text_input("Nombre del Archivo de Salida (.xyz)", value="diblock.xyz")
    if st.button("Generar Polímero Dibloque"):
        st.write(generar_polimero_dibloque(smiles_a, smiles_b, group, n_a, n_b, output_filename))

# Formulario para polímero tribloque
elif tipo_polimero == "Tribloque":
    smiles_a = st.text_input("SMILES del Monómero A", value="BrCOCBr")
    smiles_b = st.text_input("SMILES del Monómero B", value="c1cc(sc1Br)Br")
    smiles_c = st.text_input("SMILES del Monómero C", value="c1cc(oc1Br)Br")
    group = st.text_input("Grupo de Conexión", value="Br")
    n_units = st.number_input("Número de Unidades Totales", min_value=1, max_value=100, value=10)
    ratio_ab = st.slider("Proporción A/B", min_value=0.1, max_value=1.0, value=0.4, step=0.1)
    ratio_c = st.slider("Proporción C", min_value=0.1, max_value=1.0, value=0.2, step=0.1)
    output_filename = st.text_input("Nombre del Archivo de Salida (.xyz)", value="triblock.xyz")
    if st.button("Generar Polímero Tribloque"):
        st.write(generar_polimero_tribloque(smiles_a, smiles_b, smiles_c, group, n_units, ratio_ab, ratio_c, output_filename))
