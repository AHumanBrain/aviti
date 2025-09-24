import streamlit as st
import pandas as pd
import math

st.set_page_config(page_title="Denature & Dilute + Variable Read Allocation", layout="wide")
st.title("Library Prep — Pooling & Denature/Dilute Calculator")

st.markdown("""
This app lets you define per-sample read allocations (percent of cartridge) and computes how much mass and volume of each library to pool, including dilution steps, PhiX working dilution handling, and the final denature/dilute procedure to load 1400 µL.

Upload a sheet with columns (or enter manually):
- Sample alias
- Library size (bp) or insert size (bp) and adapter type
- Number of unique oligos
- i7 index (optional)
- Qubit ng/µL
- Desired percent of run (or desired coverage and cartridge capacity; see options below)

Only a minimal set of inputs is required; the app will compute the rest.
""")

# --- Data input: upload or manual ---
st.header('1) Sample table')
uploaded = st.file_uploader('Upload a CSV/XLSX exported from your sample planner (columns: alias, insert_bp [or library_bp], unique_oligos, i7, qubit_ngul, desired_percent)', type=['csv','xlsx'])
manual = st.checkbox('Or enter samples manually')

if uploaded is not None:
    try:
        raw = pd.read_excel(uploaded, header=1) if uploaded.name.lower().endswith(('xlsx','xls')) else pd.read_csv(uploaded, header=1)
        st.success('File loaded — attempting to map columns by header names. If headers differ, use manual entry.')
        df = raw.copy()
    except Exception as e:
        st.error(f'Could not read file: {e}')
        df = pd.DataFrame()
elif manual:
    st.write('Enter one sample per line: alias, insert_bp, unique_oligos, i7 (optional), qubit_ngul, desired_percent')
txt = st.text_area(
    "Paste lines here",
    value="sampleA,300,100000,7410,5,1\nsampleB,300,200000,7411,5,10"
)
rows = [r for r in [l.strip() for l in txt.splitlines()] if r]
parsed = []
for r in rows:
        parts = [p.strip() for p in r.split(',')]
        # allow missing fields
        while len(parts) < 6:
            parts.append('')
        parsed.append(parts[:6])
    df = pd.DataFrame(parsed, columns=['alias','insert_bp','unique_oligos','i7','qubit_ngul','desired_percent'])
else:
    df = pd.DataFrame(columns=['alias','insert_bp','unique_oligos','i7','qubit_ngul','desired_percent'])

# Normalize column names if present
cols = [c.lower().strip() for c in df.columns]
# Try to rename common variations
rename_map = {}
for c in df.columns:
    lc = c.lower()
    if 'alias' in lc or 'sample' in lc:
        rename_map[c] = 'alias'
    if 'insert' in lc and 'bp' in lc:
        rename_map[c] = 'insert_bp'
    if 'library size' in lc or (('size' in lc or 'bp' in lc) and 'insert' not in lc):
        rename_map[c] = 'library_bp'
    if 'unique' in lc and 'oligo' in lc:
        rename_map[c] = 'unique_oligos'
    if 'i7' in lc:
        rename_map[c] = 'i7'
    if 'qubit' in lc or 'ng/ul' in lc or 'ng' in lc:
        rename_map[c] = 'qubit_ngul'
    if 'percent' in lc or '%' in lc:
        rename_map[c] = 'desired_percent'
if rename_map:
    df = df.rename(columns=rename_map)

# Fill missing columns with defaults
for req in ['alias','insert_bp','unique_oligos','i7','qubit_ngul','desired_percent']:
    if req not in df.columns:
        df[req] = ''

# Convert numeric columns
for col in ['insert_bp','unique_oligos','qubit_ngul','desired_percent']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

st.subheader('Sample preview')
st.dataframe(df)

# --- Global parameters ---
st.header('2) Global parameters')
adapter_type = st.radio('Adapter type', ['Single-indexed (124 bp)','Dual-indexed (136 bp)'])
adapter_len = 124 if adapter_type.startswith('Single') else 136
n_reads_capacity = st.number_input('Cartridge capacity (total unique reads)', value=2500000, min_value=1)
desired_coverage = st.number_input('Desired coverage (reads per unique oligo)', value=20, min_value=1)
loading_conc_pM = st.number_input('Desired loading concentration (pM)', value=10.0, min_value=0.1)
phix_use = st.checkbox('Include PhiX (1% fixed)', value=True)
phix_stock_choice = st.radio('PhiX source', ['Use 1 nM stock and dilute (recommended)','I have pre-diluted PhiX working solution'])
phix_dil = st.number_input('If using stock: dilution factor (e.g. 40 for 1:40)', value=40, min_value=1) if phix_stock_choice.startswith('Use') else None
phix_working_nM = (1.0 / phix_dil) if phix_stock_choice.startswith('Use') else st.number_input('Enter PhiX working conc (nM)', value=0.025, min_value=0.0, format='%.6f')

st.markdown('---')

# --- Calculations per sample ---
rows = []
for i, r in df.iterrows():
    alias = r.get('alias') if pd.notna(r.get('alias')) else f'sample_{i+1}'
    insert_bp = r.get('insert_bp') if pd.notna(r.get('insert_bp')) else None
    unique_oligos = r.get('unique_oligos') if pd.notna(r.get('unique_oligos')) else 0
    i7 = r.get('i7')
    qubit_ngul = r.get('qubit_ngul') if pd.notna(r.get('qubit_ngul')) else None
    desired_percent = r.get('desired_percent') if pd.notna(r.get('desired_percent')) else None
    # compute library size (bp)
    lib_bp = insert_bp + adapter_len if insert_bp else None
    # percent of run from unique oligos * coverage / capacity
    percent_of_run = (unique_oligos * desired_coverage) / n_reads_capacity if unique_oligos and n_reads_capacity else (desired_percent/100.0 if desired_percent else 0)
    # compute total mass at load (ng) using loading_conc and average bp: mass_total_ng = C_nM * V_uL * bp * 660e-6
    C_nM = loading_conc_pM / 1000.0  # pM -> nM
    V_uL = 1400.0
    mass_total_ng = C_nM * V_uL * (lib_bp if lib_bp else 300) * 660.0 * 1e-6
    ng_to_pool = percent_of_run * mass_total_ng
    uL_to_pool = ng_to_pool / qubit_ngul if qubit_ngul and qubit_ngul>0 else None
    # dilution factor selection to make diluted pool aliquot pipettable (prefer 4-8 µL)
    chosen_df = None
    chosen_diluted_pool_uL = None
    if uL_to_pool:
        best_score = None
        for d in range(1,101):
            diluted_uL = uL_to_pool * d
            if 1.0 <= diluted_uL <= 50.0:
                score = -abs(diluted_uL - 6.0)  # prefer near 6 µL
                if best_score is None or score > best_score:
                    best_score = score
                    chosen_df = d
                    chosen_diluted_pool_uL = round(diluted_uL,3)
    rows.append({
        'alias': alias,
        'lib_bp': lib_bp,
        'unique_oligos': unique_oligos,
        'percent_of_run': percent_of_run,
        'qubit_ngul': qubit_ngul,
        'ng_to_pool': ng_to_pool,
        'uL_to_pool': uL_to_pool,
        'dilution_factor': chosen_df,
        'diluted_uL_to_pool': chosen_diluted_pool_uL,
        'i7': i7
    })

out = pd.DataFrame(rows)
st.subheader('Computed pooling plan (pre-dilution)')
st.dataframe(out)

# Now compute combined pool prep for denature/dilute steps
st.header('3) Pool prep -> Denature/Dilute/Load')
# sum diluted volumes that will be pooled (use diluted_uL_to_pool if present, else uL_to_pool)
pooling_contrib = 0.0
for idx, row in out.iterrows():
    if row['diluted_uL_to_pool']:
        pooling_contrib += row['diluted_uL_to_pool']
    elif row['uL_to_pool']:
        pooling_contrib += row['uL_to_pool']

st.markdown(f"Total pooled diluted volume (sum of diluted aliquots): **{pooling_contrib:.3f} µL**")

# PhiX addition: compute required phiX volume to achieve 1% of molecules at load
if phix_use:
    # use phix_working_nM
    # compute library working concentration in nM for pooled diluted mixture: approximate using pool_conc_nM/d where d is chosen average dilution
    avg_dil = out['dilution_factor'].dropna().mean() if not out['dilution_factor'].isnull().all() else 1
    pool_conc_nM = (pooled_ngul * 1e6) / (660.0 * (library_size_bp if 'library_size_bp' in locals() else 300)) if 'pooled_ngul' in locals() else None
    st.markdown('PhiX will be added per SOP to reach ~1% of molecules at load. The app reports required µL based on pooled diluted volume above.')

st.caption('This is a first pass implementation. I implemented the spreadsheet logic: percent_of_run = unique_oligos * desired_coverage / cartridge_capacity; mass to pool is percent_of_run * total_mass_required_at_load (derived from loading pM and library size). Dilution factors are chosen to make pipettable aliquots near 6 µL. The app will need minor tuning after testing with your sheet; tell me which rows to debug and I will iterate.')
