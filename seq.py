import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="Denature & Dilute Calculator", layout="wide")
st.title("Library Prep — Denature & Dilute Web Calculator (Streamlit)")

# --- Section 1: Sample Input ---
st.header("1) Input Library Concentrations")
uploaded = st.file_uploader("Upload a CSV/XLSX with library concentrations (ng/µL)", type=["csv","xlsx"])

if uploaded is not None:
    try:
        if uploaded.name.endswith('.csv'):
            df = pd.read_csv(uploaded, header=None)
        else:
            df = pd.read_excel(uploaded, header=None)
        df.columns = ['Concentration_ng_per_uL']
    except Exception as e:
        st.error(f"Could not read file: {e}")
        df = pd.DataFrame({'Concentration_ng_per_uL':[2.11,2.39,2.34]})
else:
    df = pd.DataFrame({'Concentration_ng_per_uL':[2.11,2.39,2.34]})

st.dataframe(df)

# --- Section 2: Pooling ---
st.header("2) Pooling Calculation")
target_ng = st.number_input('Target mass to pool per library (ng)', value=30.0, min_value=0.0, step=1.0)
df['uL_to_pool'] = (target_ng / df['Concentration_ng_per_uL']).round(6)
df['ng_pooled'] = (df['Concentration_ng_per_uL'] * df['uL_to_pool']).round(3)

st.dataframe(df)

pool_volume = df['uL_to_pool'].sum()
st.markdown(f"**Total pooled volume (µL):** {pool_volume:.3f}")

# --- Section 3: Denature & Dilution ---
st.header("3) Denature & Dilution Workflow")

st.subheader("a) Pool Libraries & Qubit Quant")
expected_conc = (df['ng_pooled'].sum() / pool_volume) if pool_volume > 0 else 0
st.markdown(f"**Expected pooled concentration:** {expected_conc:.2f} ng/µL")
actual_conc = st.number_input('Enter actual pooled concentration from Qubit (ng/µL)', value=expected_conc, min_value=0.0)
if actual_conc < expected_conc * 0.8 or actual_conc > expected_conc * 1.2:
    st.warning('Actual concentration differs >20% from expected — double-check quantification!')

st.subheader("b) Dilute Library Pool")
desired_ng_per_ul = st.number_input('Target library pool concentration (ng/µL) for loading', value=actual_conc, min_value=0.0)
library_dilution_factor = actual_conc / desired_ng_per_ul if desired_ng_per_ul > 0 else 1
st.markdown(f"Dilution factor for library pool: **{library_dilution_factor:.2f}x**")

st.subheader("c) Dilute PhiX")
use_phix = st.checkbox("Use PhiX (spike-in)", value=True)
if use_phix:
    phix_conc = st.number_input('PhiX stock concentration (nM or ng/µL depending on SOP)', value=10.0)
    phix_target = st.number_input('Target PhiX fraction (%)', value=1.0)
    st.markdown(f"Dilute PhiX to match desired spike-in concentration based on overall loading conc ({desired_ng_per_ul} ng/µL).")

st.subheader("d) Pool Library + PhiX")
st.markdown("Combine diluted library pool and diluted PhiX in desired ratio before denaturation.")

st.subheader("e) Denature with NaOH")
naoh_volume = st.number_input('Volume of 0.2 N NaOH to add (µL)', value=pool_volume)
st.markdown(f"Add {naoh_volume:.2f} µL of 0.2 N NaOH to the combined library+PhiX pool, mix, spin, incubate 5 min.")

st.subheader("f) Neutralize with Tris-HCl")
st.markdown(f"Add {naoh_volume:.2f} µL of 0.2 M pH 7 Tris-HCl to neutralize.")

st.subheader("g) Bring to Final Volume")
final_volume = 1400.0
loading_buffer = st.number_input('Volume of loading buffer to add (µL)', value=max(0.0, final_volume - (pool_volume + naoh_volume*2)))
st.markdown(f"Bring total volume up to **{final_volume} µL** with loading buffer.")

st.subheader("h) Load")
st.markdown("Load all 1400 µL into cartridge.")

# Export results
st.header('Export')
buf = io.BytesIO()
out_df = df.copy()
out_df['Total_pooled_volume_uL'] = pool_volume
out_df.to_csv(buf, index=False)
buf.seek(0)
st.download_button('Download pool plan as CSV', data=buf, file_name='pool_plan.csv')

st.caption('Verify all calculations against your SOP and instrument documentation before proceeding.')
