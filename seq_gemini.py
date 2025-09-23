import streamlit as st
import pandas as pd
import math

# --- TROUBLESHOOTING GUIDE ---
#
# The error "TypeError: Failed to fetch dynamically imported module..." is typically
# a front-end or deployment issue, not a bug in your Python code logic below.
# It means the user's browser could not download a necessary JavaScript file from
# Streamlit's server.
#
# Here are the most common causes and solutions:
#
# 1. Browser Cache: The browser might have a stale or corrupted version of the app's
#    files.
#    >> SOLUTION: Do a "hard refresh" (Ctrl+Shift+R or Cmd+Shift+R) and clear
#       your browser cache.
#
# 2. Streamlit Cloud Environment: Sometimes the Streamlit Cloud environment can get
#    into a bad state.
#    >> SOLUTION: In your Streamlit Cloud workspace, find the app and use the menu
#       to "Reboot" it. This often resolves transient issues.
#
# 3. Dependency Mismatch: A mismatch between the `streamlit` library version (or other
#    libraries) in your `requirements.txt` and what's expected can cause errors.
#    >> SOLUTION: Ensure your `requirements.txt` file pins specific versions, for example:
#       streamlit==1.32.2
#       pandas==2.2.1
#       openpyxl==3.1.2  # Required for .xlsx file support in pandas
#
# 4. Network Issues or Browser Extensions: A firewall, ad-blocker, or other browser
#    extension could be blocking the script from loading.
#    >> SOLUTION: Try accessing the app in an incognito/private browser window or
#       disabling extensions to see if the problem persists.
#
# Your Python code logic below is robust and unlikely to be the direct cause.
# --- END TROUBLESHOOTING GUIDE ---


st.set_page_config(page_title="Denature & Dilute Calculator", layout="wide")
st.title("Library Prep — Denature & Dilute Web Calculator (Streamlit)")

# --- Section 1: Sample Input ---
st.header("1) Input Library Concentrations (ng/µL)")
uploaded = st.file_uploader("Upload a CSV/XLSX with library concentrations (one column of numbers, ng/µL)", type=["csv","xlsx"])
manual_text = st.text_area(
    "Or manually enter concentrations (one per line)",
    placeholder="2.11\n2.39\n2.34"
)

# --- Code Refinement ---
# Defined a default dataframe once to avoid repeating code.
default_df = pd.DataFrame({'Concentration_ng_per_uL': [2.11, 2.39, 2.34]})
df = None # Initialize df as None

# load concentrations (single-column) into dataframe
if uploaded is not None:
    try:
        if uploaded.name.lower().endswith('.csv'):
            df = pd.read_csv(uploaded, header=None)
        else:
            # Note: pandas needs the 'openpyxl' library to read .xlsx files.
            # Make sure it's in your requirements.txt file.
            df = pd.read_excel(uploaded, header=None)
        df.columns = ['Concentration_ng_per_uL']
        df['Concentration_ng_per_uL'] = pd.to_numeric(df['Concentration_ng_per_uL'], errors='coerce')
        df.dropna(inplace=True) # Remove rows that couldn't be converted to numbers
    except Exception as e:
        st.error(f"Could not read file: {e}")
        df = default_df

elif manual_text and manual_text.strip():
    try:
        vals = [float(x.strip()) for x in manual_text.strip().splitlines() if x.strip()]
        df = pd.DataFrame({'Concentration_ng_per_uL': vals})
    except Exception as e:
        st.error(f"Could not parse manual entries: {e}")
        df = default_df

else:
    df = default_df

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

# Ensure variables are initialized
chosen_L = None
chosen_P = None

# User inputs minimal: library size and desired loading conc
st.subheader('User inputs (minimal)')
library_size_bp = st.number_input('Average library size (bp)', value=400, min_value=1)
loading_conc_pM = st.number_input('Desired loading concentration (pM)', value=10.0, min_value=0.1)
# effective loading (99%)
effective_loading_pM = 0.99 * loading_conc_pM
st.markdown(f"Effective loading concentration used for calculations: **{effective_loading_pM:.3f} pM**")

# Compute expected pooled concentration from pool plan (ng/µL) and convert to nM
expected_pooled_ngul = (df['ng_pooled'].sum() / pool_volume) if pool_volume > 0 else 0.0
st.markdown(f"Expected pooled concentration (from pooling plan): **{expected_pooled_ngul:.4f} ng/µL**")
actual_qubit = st.number_input('(Optional) Enter actual pooled concentration from Qubit (ng/µL) — leave as expected to use calculated value', value=round(expected_pooled_ngul,6), format="%.6f")
pooled_ngul = actual_qubit if actual_qubit and actual_qubit > 0 else expected_pooled_ngul

# convert ng/µL to nM using library size
pool_conc_nM = pooled_ngul * 0.8 * 1e6 / 660.0 / library_size_bp
st.markdown(f"Pooled library concentration: **{pool_conc_nM:.4f} nM** (based on {pooled_ngul:.6f} ng/µL and {library_size_bp} bp)")

# PhiX working concentration options
st.subheader('PhiX working solution')
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

# --- Formulae from SOP / Spreadsheet ---
# For each dilution factor d, calculate required library and PhiX volumes
best = None
candidates = []
for d in range(1,51):
    C_L = pool_conc_nM / d if d > 0 else 0  # diluted library concentration (nM)
    if C_L <= 0:
        continue
    # Using sheet formula: P = (f * L * C_L) / ((1-f) * C_P)
    for L in [x/10.0 for x in range(10,101)]:  # 1.0 to 10.0 µL in 0.1 increments
        denom = (1.0 - phix_fraction) * phix_working_nM
        if denom <= 0:
            continue
        P = (phix_fraction * L * C_L) / denom
        if P < 1.0 or P > 10.0:
            continue
        pre_total = L + P
        if pre_total <= 0 or pre_total > 30:
            continue
        # effective concentration achieved
        achieved_pM = (L * C_L + P * phix_working_nM) / pre_total * 1000.0  # nM->pM
        diff = abs(achieved_pM - effective_loading_pM)
        # Score penalizes deviation from target and favors L near 6.0 uL
        score = -diff - abs(6.0 - L)
        candidates.append((score, d, round(L,2), round(P,2), round(pre_total,2), achieved_pM))

if candidates:
    candidates.sort(reverse=True)
    best = candidates[0]

st.subheader('Automated Denature & Dilute Plan')
if best:
    _, best_d, best_L, best_P, best_pre, achieved_pM = best
    st.success(
        f"**Optimal plan found:**\n\n"
        f"1. **Dilute pooled library 1:{best_d}** (working conc. will be {pool_conc_nM/best_d:.4f} nM).\n\n"
        f"2. Pipette **{best_L} µL** of diluted library.\n\n"
        f"3. Pipette **{best_P} µL** of PhiX working solution.\n\n"
        f"This gives a pre-denature total of **{best_pre} µL** and achieves ~**{achieved_pM:.2f} pM** (target: {effective_loading_pM:.2f} pM)."
    )
    chosen_L = best_L
    chosen_P = best_P
else:
    st.error('Could not find pipettable volumes that satisfy constraints with current inputs. Try adjusting loading concentration or library size.')

# Compute NaOH and Tris volumes (equal to pre-denature total) and final buffer top-up
pre_denature_vol = round((chosen_L or 0) + (chosen_P or 0),3)
naoh_vol = pre_denature_vol
tris_vol = naoh_vol
small_denature_total = round(pre_denature_vol + naoh_vol + tris_vol,3)
loading_buffer_to_add = round(max(0.0, 1400.0 - small_denature_total),3)

st.subheader('Final Denature & Neutralize Steps')
denature_plan = {
    "Pre-denature volume (library + PhiX)": f"{pre_denature_vol} µL",
    "Add 0.2 N NaOH (mix, spin, incubate 5 min)": f"{naoh_vol} µL",
    "Add 0.2 M pH 7 Tris-HCl to neutralize": f"{tris_vol} µL",
    "Total after neutralization": f"{small_denature_total} µL",
    "Add loading buffer to bring to 1400 µL": f"{loading_buffer_to_add} µL",
    "Final volume to load into cartridge": "1400 µL"
}
st.table(pd.DataFrame(denature_plan.items(), columns=['Step', 'Volume']))


st.caption('This calculator provides a suggested plan based on your inputs. All calculations should be verified in the lab before proceeding.')
