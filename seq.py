import streamlit as st
import pandas as pd
import math

st.set_page_config(page_title="Denature & Dilute Calculator", layout="wide")
st.title("Library Prep — Denature & Dilute Web Calculator (Streamlit)")

# --- Section 1: Sample Input ---
st.header("1) Input Library Concentrations (ng/µL)")
uploaded = st.file_uploader("Upload a CSV/XLSX with library concentrations (one column of numbers, ng/µL)", type=["csv","xlsx"])
manual_text = st.text_area("Or manually enter concentrations (one per line)", placeholder="2.11\n2.39\n2.34")

# load concentrations (single-column) into dataframe
if uploaded is not None:
    try:
        if uploaded.name.lower().endswith('.csv'):
            df = pd.read_csv(uploaded, header=None)
        else:
            df = pd.read_excel(uploaded, header=None)
        df.columns = ['Concentration_ng_per_uL']
        df['Concentration_ng_per_uL'] = pd.to_numeric(df['Concentration_ng_per_uL'], errors='coerce')
    except Exception as e:
        st.error(f"Could not read file: {e}")
        df = pd.DataFrame({'Concentration_ng_per_uL': [2.11, 2.39, 2.34]})

elif manual_text and manual_text.strip():
    try:
        vals = [float(x.strip()) for x in manual_text.strip().splitlines() if x.strip()]
        df = pd.DataFrame({'Concentration_ng_per_uL': vals})
    except Exception as e:
        st.error(f"Could not parse manual entries: {e}")
        df = pd.DataFrame({'Concentration_ng_per_uL': [2.11, 2.39, 2.34]})

else:
    df = pd.DataFrame({'Concentration_ng_per_uL': [2.11, 2.39, 2.34]})

st.subheader('Input preview')
st.dataframe(df)

# --- Section 2: Pooling ---
st.header("2) Pooling calculation")
st.write("This calculates µL to pool per library to reach a target mass (ng) per library in the pool.")

target_ng = st.number_input('Target mass to pool per library (ng)', value=30.0, min_value=0.0, step=1.0)
# guard against zeros
if (df['Concentration_ng_per_uL'] <= 0).any():
    st.warning('One or more concentration values are zero or invalid — check inputs.')

# compute uL to pool and ng pooled
df['uL_to_pool'] = (target_ng / df['Concentration_ng_per_uL']).round(6)
df['ng_pooled'] = (df['Concentration_ng_per_uL'] * df['uL_to_pool']).round(3)

st.subheader('Per-sample pooling plan')
st.dataframe(df)

pool_volume = df['uL_to_pool'].sum()
st.markdown(f"**Total pooled volume (µL):** {pool_volume:.3f}")

# --- Section 3: Denature & Dilution Workflow ---
st.header("3) Denature & Dilution Workflow (stepwise)")

# a) Pool libraries and Qubit
st.subheader('a) Pool libraries & Qubit quant')
expected_conc = (df['ng_pooled'].sum() / pool_volume) if pool_volume > 0 else 0.0
st.markdown(f"**Expected pooled concentration:** {expected_conc:.3f} ng/µL")
actual_conc = st.number_input('Enter actual pooled concentration from Qubit (ng/µL)', value=round(expected_conc,3), min_value=0.0, format="%.4f")
# warn if off
if actual_conc > 0 and (actual_conc < expected_conc * 0.8 or actual_conc > expected_conc * 1.2):
    st.warning('Actual concentration differs by >20% from the expected concentration. Double-check quantification or pooling.')

# b) Dilute library pool accordingly
st.subheader('b) Dilute library pool to target concentration (prep small aliquot)')
st.write('We prepare a **small diluted aliquot** of the pooled library (defaults to 10 µL) at your target concentration so downstream denature steps stay in the 1–10 µL range.')

desired_ng_per_ul = st.number_input('Target library pool concentration (ng/µL) to make for loading (e.g. 0.4)', value=round(actual_conc,3), min_value=0.0, format="%.4f")
avg_bp = st.number_input('Average library fragment length (bp) — used to convert ng/µL ↔ nM', value=350, min_value=1)

prep_diluted_target = 10.0
if desired_ng_per_ul > 0 and actual_conc > 0:
    vol_pool_needed_for_prep = prep_diluted_target * desired_ng_per_ul / actual_conc
else:
    vol_pool_needed_for_prep = float('nan')

if not math.isnan(vol_pool_needed_for_prep) and vol_pool_needed_for_prep > pool_volume and pool_volume > 0:
    max_prep_diluted = actual_conc * pool_volume / desired_ng_per_ul
    prep_diluted_vol = round(min(prep_diluted_target, max_prep_diluted), 3)
    vol_pool_needed = round(prep_diluted_vol * desired_ng_per_ul / actual_conc, 3)
    buffer_vol = round(prep_diluted_vol - vol_pool_needed, 3)
    st.warning(f"Not enough pooled volume to prepare {prep_diluted_target} µL. Will prepare {prep_diluted_vol} µL (uses {vol_pool_needed} µL of pooled library).")
else:
    prep_diluted_vol = prep_diluted_target
    vol_pool_needed = round(vol_pool_needed_for_prep, 3) if not math.isnan(vol_pool_needed_for_prep) else float('nan')
    buffer_vol = round(prep_diluted_vol - vol_pool_needed, 3) if not math.isnan(vol_pool_needed_for_prep) else float('nan')

st.markdown(f"**Prepare {prep_diluted_vol} µL of diluted pooled library at {desired_ng_per_ul} ng/µL**")
st.markdown(f"- Pipette **{vol_pool_needed} µL** of the pooled library (from your Qubit pool) and add **{buffer_vol} µL** of buffer/diluent to reach {prep_diluted_vol} µL total.")
if vol_pool_needed > pool_volume:
    st.error('The requested volume of pooled library to make the diluted aliquot exceeds the available pooled library volume.')

# c) PhiX dilution options
st.subheader('c) PhiX preparation and dilution')
phix_enabled = st.checkbox('Include PhiX spike-in', value=True)
if phix_enabled:
    phix_source = st.radio('Are you using stock PhiX (1 nM) or a pre-diluted PhiX?', ['Stock 1 nM PhiX', 'Pre-diluted PhiX'])
    if phix_source == 'Stock 1 nM PhiX':
        phix_dilution_factor = st.number_input('Dilution factor to prepare working PhiX from 1 nM stock (e.g. 40 = 1:40 dilution)', value=40, min_value=1)
        phix_working_conc = 1.0 / phix_dilution_factor
        st.markdown(f"Prepare working PhiX: dilute stock 1:{phix_dilution_factor} → working concentration **{phix_working_conc:.4f} nM**.")
    else:
        phix_working_conc = st.number_input('Enter concentration of your pre-diluted PhiX (nM)', value=0.025, min_value=0.0, format="%.4f")

    phix_target_percent = st.number_input('Target PhiX % of total molecules at load (percent)', value=1.0, min_value=0.0, max_value=100.0, step=0.1)

    # Convert desired library concentration (ng/µL) -> nM
    if desired_ng_per_ul > 0 and avg_bp > 0:
        library_nM = desired_ng_per_ul * 1e6 / (660.0 * avg_bp)
    else:
        library_nM = 0.0

    st.markdown(f"Converted diluted library: **{library_nM:.3f} nM** (based on {desired_ng_per_ul} ng/µL and {avg_bp} bp)")

    # find volumes for library and PhiX working to achieve requested %
    f = phix_target_percent / 100.0
    chosen_L = None
    chosen_P = None
    achieved_f = None
    max_L_available = min(prep_diluted_vol, 10.0)
    for L in [round(x*0.5,2) for x in range(2,21)]:
        if L > max_L_available:
            break
        denom = phix_working_conc * (1.0 - f)
        if denom <= 0 or library_nM <= 0:
            continue
        P = (f * L * library_nM) / denom
        if P >= 1.0 and P <= 10.0 and (L + P) <= 10.0:
            chosen_L = L
            chosen_P = round(P, 3)
            break
    if chosen_L is None:
        L = round(max_L_available, 3)
        denom = phix_working_conc * (1.0 - f)
        if denom > 0 and library_nM > 0:
            P_raw = (f * L * library_nM) / denom
            if not math.isnan(P_raw):
                P = max(1.0, min(P_raw, 10.0))
                if (L + P) > 10.0:
                    P = round(max(1.0, 10.0 - L), 3)
                chosen_L = L
                chosen_P = round(P, 3)
                achieved_f = (chosen_P * phix_working_conc) / (chosen_L * library_nM + chosen_P * phix_working_conc)
    if chosen_L and chosen_P:
        st.markdown(f"**Aliquot for denature step:** {chosen_L} µL diluted library + {chosen_P} µL PhiX working ({phix_working_conc:.4f} nM).")
        if achieved_f:
            st.markdown(f"Requested PhiX %: {phix_target_percent}%, achieved with pipettable volumes: {achieved_f*100:.3f}%.")

# d) Pool the library pool and PhiX together
st.subheader('d) Pool the diluted library aliquot + PhiX')

# e) Denature with NaOH (equi-vol)
st.subheader('e) Denature with 0.2 N NaOH')
if phix_enabled and chosen_L and chosen_P:
    pre_denature_vol = round(chosen_L + chosen_P, 3)
elif chosen_L:
    pre_denature_vol = round(chosen_L, 3)
else:
    pre_denature_vol = 0.0

naoh_volume = pre_denature_vol
st.markdown(f"Add **{naoh_volume} µL** of 0.2 N NaOH to the aliquot (mix, spin, incubate 5 minutes).")

# f) Neutralize with Tris-HCl
st.subheader('f) Neutralize with 0.2 M Tris-HCl (pH 7.0)')
tris_volume = naoh_volume
st.markdown(f"Add **{tris_volume} µL** of 0.2 M pH 7 Tris-HCl to neutralize.")

# g) Bring to final loading volume
st.subheader('g) Bring to final loading volume with buffer')
final_loading_volume = 1400.0
small_denature_total = round(pre_denature_vol + naoh_volume + tris_volume, 3)
loading_buffer_to_add = round(max(0.0, final_loading_volume - small_denature_total), 3)
st.markdown(f"After neutralization: **{small_denature_total} µL**. Add **{loading_buffer_to_add} µL** loading buffer → total **{final_loading_volume} µL**.")

# h) Load
st.subheader('h) Load')
st.markdown('Load the entire **1,400 µL** into the cartridge.')

st.caption('All pipetting instructions are calculated and read-only. Inputs: library Qubit conc, desired library conc, PhiX % and dilution settings. The app enforces 1–10 µL ranges for constituents in denature step.')
