import streamlit as st
import pandas as pd
import math

st.set_page_config(page_title="Denature & Dilute Calculator", layout="wide")
st.title("Library Prep — Denature & Dilute Web Calculator (Streamlit)")

# --- Section 1: Sample Input ---
st.header("1) Input Library Concentrations (ng/µL)")
uploaded = st.file_uploader("Upload a CSV/XLSX with library concentrations (one column of numbers, ng/µL)", type=["csv","xlsx"])
manual_text = st.text_area("Or manually enter concentrations (one per line)", placeholder="2.11
2.39
2.34")

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

# compute uL to pool and ng pooled
df['uL_to_pool'] = (target_ng / df['Concentration_ng_per_uL']).round(6)
df['ng_pooled'] = (df['Concentration_ng_per_uL'] * df['uL_to_pool']).round(3)

st.subheader('Per-sample pooling plan')
st.dataframe(df)

pool_volume = df['uL_to_pool'].sum()
st.markdown(f"**Total pooled volume (µL):** {pool_volume:.3f}")

# --- Section 3: Denature & Dilution Workflow ---
st.header("3) Denature & Dilute — automated volumes")

# Ensure variables
chosen_L = None
chosen_P = None
achieved_f = None

# User inputs minimal: insert size and desired loading conc
st.subheader('User inputs (minimal)')
insert_length_bp = st.number_input('Library insert size (bp)', value=400, min_value=1)
loading_conc_pM = st.number_input('Desired loading concentration (pM)', value=10.0, min_value=0.1)
# effective loading (99%)
effective_loading_pM = 0.99 * loading_conc_pM
st.markdown(f"Effective loading concentration used for calculations: **{effective_loading_pM:.3f} pM**")

# Compute expected pooled concentration from pool plan (ng/µL) and convert to nM
expected_pooled_ngul = (df['ng_pooled'].sum() / pool_volume) if pool_volume > 0 else 0.0
st.markdown(f"Expected pooled concentration (from pooling plan): **{expected_pooled_ngul:.4f} ng/µL**")
actual_qubit = st.number_input('(Optional) Enter actual pooled concentration from Qubit (ng/µL) — leave as expected to use calculated value', value=round(expected_pooled_ngul,6), format="%.6f")
pooled_ngul = actual_qubit if actual_qubit and actual_qubit > 0 else expected_pooled_ngul
# convert to nM: nM = (ng/µL * 1e6) / (660 * bp)
pool_conc_nM = pooled_ngul * 1e6 / (660.0 * insert_length_bp)
st.markdown(f"Pooled library concentration: **{pool_conc_nM:.4f} nM** (based on {pooled_ngul:.6f} ng/µL and {insert_length_bp} bp)")

# PhiX working concentration options
st.subheader('PhiX working solution')
phix_default_from_stock = True
phix_choice = st.radio('PhiX source', ['Use 1 nM stock and dilute (recommended)','I have a pre-diluted PhiX working solution'])
if phix_choice.startswith('Use 1 nM'):
    phix_dil_factor = st.number_input('Dilution factor for stock -> working (e.g. 40 for 1:40)', value=40, min_value=1)
    phix_working_nM = 1.0 / phix_dil_factor
else:
    phix_working_nM = st.number_input('Enter PhiX working concentration (nM)', value=0.025, min_value=0.0, format="%.6f")
st.markdown(f"PhiX working concentration used: **{phix_working_nM:.6f} nM**")

# Target PhiX fraction (fixed at 1% per SOP)
phix_fraction = 0.01
st.markdown(f"PhiX target fraction of molecules at load: **{phix_fraction*100:.2f}% (fixed)**")

# Now find a dilution factor for the pooled library (integer 1..50) that yields pipettable aliquot L (µL)
# Equation: final_pM = L * C_L * 0.7142857 ; where C_L = pool_conc_nM / d ; target final_pM = effective_loading_pM
# => L = (effective_loading_pM / 0.7142857) / C_L
K = effective_loading_pM * 1.4  # since 1/0.7142857 ≈ 1.4
best = None
candidates = []
for d in range(1,51):
    C_L = pool_conc_nM / d if d>0 else 0
    if C_L <= 0:
        continue
    L = K / C_L
    if L <= 0:
        continue
    # compute required P for PhiX: P = (f * L * C_L) / ((1-f) * C_P)
    denom = (1.0 - phix_fraction) * phix_working_nM
    if denom <= 0:
        P = float('nan')
    else:
        P = (phix_fraction * L * C_L) / denom
    pre_total = L + P
    # check constraints: individual volumes between 1 and 10, prefer L in 4-8
    ok = (1.0 <= L <= 10.0) and (1.0 <= P <= 10.0) and (pre_total <= 30.0)
    score = 0
    # scoring: prefer ok, then prefer L near 6, smaller pre_total
    score += (100 if ok else 0)
    score -= abs(6.0 - L)
    score -= pre_total * 0.01
    candidates.append((score, d, round(L,3), round(P,3), round(pre_total,3)))

if candidates:
    candidates.sort(reverse=True)
    best = candidates[0]

if best:
    _, best_d, best_L, best_P, best_pre = best
    st.markdown(f"Selected pooled-library dilution: **1:{best_d}** → working conc {pool_conc_nM/best_d:.4f} nM")
    st.markdown(f"Pipette **{best_L} µL** of diluted pooled library and **{best_P} µL** of PhiX working solution into the denature mix. (Pre-denature volume: {best_pre} µL)")
    chosen_L = best_L
    chosen_P = best_P
else:
    st.error('Could not find pipettable volumes that satisfy constraints with current inputs. Try adjusting loading concentration or insert size.')

# Compute NaOH and Tris volumes (equal to pre-denature total) and final buffer top-up
pre_denature_vol = round((chosen_L or 0) + (chosen_P or 0),3)
naoh_vol = pre_denature_vol
tris_vol = naoh_vol
small_denature_total = round(pre_denature_vol + naoh_vol + tris_vol,3)
loading_buffer_to_add = round(max(0.0, 1400.0 - small_denature_total),3)

st.subheader('Denature & neutralize steps (read-only)')
st.markdown(f"Pre-denature volume (library + PhiX): **{pre_denature_vol} µL**")
st.markdown(f"Add **{naoh_vol} µL** of 0.2 N NaOH (mix, spin, incubate 5 min).")
st.markdown(f"Add **{tris_vol} µL** of 0.2 M pH 7 Tris-HCl to neutralize.")
st.markdown(f"Total after neutralization: **{small_denature_total} µL**")
st.markdown(f"Add **{loading_buffer_to_add} µL** of loading buffer to bring to **1400 µL** and load entire volume into cartridge.")

st.caption('Only visible user inputs: library insert size and desired loading concentration. All other volumes are calculated to meet pipetting constraints (1–10 µL per constituent) and SOP rules (PhiX 1%). Verify before proceeding.')
