import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="Denature & Dilute Calculator", layout="wide")
st.title("Library Prep — Denature & Dilute Web Calculator (Streamlit)")
st.markdown(
    """
    Upload a table of sample concentrations (ng/µL) or enter them manually. The app computes the uL to pool for a chosen target mass per library (ng) or lets you set a target pooled concentration and do stepwise denature + dilution calculations.

    **Notes:** Default protocol parameters are *editable* — change them to match your instrument (Aviti/MiSeq) or SOP.
    """
)

# --- Input: samples ---
st.header("1) Samples / Input")
col1, col2 = st.columns([2,1])
with col1:
    uploaded = st.file_uploader("Upload a CSV/XLSX with at least columns: Sample, Concentration_ng_per_uL", type=["csv","xlsx"])
    st.write("Or paste a simple two-column table (Sample, concentration)")
    txt = st.text_area("Paste sample rows (one per line: name,conc)", height=80)

with col2:
    st.write("Example input format:")
    st.code("Sample1,2.34\nSample2,10.1\nSample3,12.1")

# load dataframe
df = None
if uploaded is not None:
    try:
        if uploaded.name.endswith('.csv'):
            df = pd.read_csv(uploaded)
        else:
            df = pd.read_excel(uploaded)
    except Exception as e:
        st.error(f"Could not read file: {e}")

if txt.strip() and df is None:
    try:
        df = pd.read_csv(io.StringIO(txt), header=None)
        if df.shape[1] >= 2:
            df = df.iloc[:,0:2]
            df.columns = ['Sample','Concentration_ng_per_uL']
    except Exception as e:
        st.error(f"Could not parse pasted text: {e}")

if df is None:
    # blank template
    df = pd.DataFrame({
        'Sample': ['Sample_1','Sample_2','Sample_3'],
        'Concentration_ng_per_uL': [2.11, 2.39, 2.34]
    })

# clean
if 'Concentration_ng_per_uL' in df.columns:
    df['Concentration_ng_per_uL'] = pd.to_numeric(df['Concentration_ng_per_uL'], errors='coerce')
else:
    # try to guess a second column
    cols = df.columns.tolist()
    if len(cols) >= 2:
        df.columns = ['Sample','Concentration_ng_per_uL']
        df['Concentration_ng_per_uL'] = pd.to_numeric(df['Concentration_ng_per_uL'], errors='coerce')

st.dataframe(df)

# --- Pooling calculation ---
st.header("2) Pooling calculation")
mode = st.radio("Choose pooling target mode", ['Target mass per library (ng)','Target pooled concentration (nM) — optional'])

if mode.startswith('Target mass'):
    target_ng = st.number_input('Target mass to pool per library (ng)', value=30.0, min_value=0.0, step=1.0)
    df['uL_to_pool'] = (target_ng / df['Concentration_ng_per_uL']).round(6)
    df['ng_pooled'] = (df['Concentration_ng_per_uL'] * df['uL_to_pool']).round(3)
    st.write('uL to pool (calculated from target mass)')
else:
    st.info('Target pooled concentration mode uses molecular conversions and requires average fragment length in bp and library molar mass assumptions.')
    target_nM = st.number_input('Target pooled concentration (nM)', value=4.0, min_value=0.0)
    avg_bp = st.number_input('Average library insert length (bp)', value=350)
    # Convert ng/ul -> nM: nM = (ng/ul * 1e6) / (660 * bp)
    df['nM_per_uL'] = (df['Concentration_ng_per_uL'] * 1e6) / (660.0 * avg_bp)
    # uL needed to reach target nM in final volume: This is an approximation; we present uL to contribute to 1 uL-equivalent
    # We will compute uL to pool for a hypothetical 1 uL final volume of library material; better to use the mass-based method.
    df['uL_to_pool'] = (target_nM / df['nM_per_uL']).replace([pd.NA, pd.Inf], None).round(6)
    df['nM_pooled'] = (df['nM_per_uL'] * df['uL_to_pool']).round(3)

st.dataframe(df[['Sample','Concentration_ng_per_uL','uL_to_pool','ng_pooled' if 'ng_pooled' in df.columns else 'nM_pooled']])

total_uL = df['uL_to_pool'].sum()
st.markdown(f"**Total pooled volume (sum of uL to pool):** {total_uL:.3f} µL")

# Option to add pooling dead volume / make-up
dead_vol = st.number_input('Pool dead volume / extra µL (to account for pipetting)', value=5.0)
pool_volume = total_uL + dead_vol
st.markdown(f"**Planned pooled volume including dead volume:** {pool_volume:.2f} µL")

# --- Denature & Dilute ---
st.header("3) Denature & Dilution steps")
st.markdown('Editable protocol parameters — change to match your SOP or instrument.')
colA, colB, colC = st.columns(3)
with colA:
    vol_library_to_use = st.number_input('Volume of pooled library used for denaturation (µL)', value=min(50.0, pool_volume))
    naoh_conc_N = st.number_input('NaOH concentration (N) to add (e.g. 0.2)', value=0.2)
    naoh_volume = st.number_input('Volume of NaOH to add (µL)', value=vol_library_to_use)  # often equal volume
with colB:
    incubation_min = st.number_input('Incubation time (min)', value=5)
    neutralize_with = st.selectbox('Neutralize with', ['Tris-HCl pH 7.0 (0.2 M)', 'Dilution Buffer (HT1)', 'Custom'])
    tris_volume = st.number_input('Volume of Tris-HCl (µL) to add after neutralization', value=1200.0)
with colC:
    final_loading_volume = st.number_input('Final volume for loading (µL)', value=1400.0)
    final_loading_conc_nM = st.number_input('Target final loading concentration (pM or nM as desired)', value=12.5)

# Basic denature calculation (conservative)
st.subheader('Calculated denature mix')
st.markdown(f"- Library aliquot: **{vol_library_to_use:.2f} µL**")
st.markdown(f"- Add NaOH {naoh_conc_N} N — volume: **{naoh_volume:.2f} µL**, incubate **{incubation_min} min** at RT")
st.markdown(f"- Then add {tris_volume:.2f} µL {neutralize_with} to neutralize/dilute to a final volume of approx **{vol_library_to_use + naoh_volume + tris_volume:.2f} µL**")

# Final computed available volume for loading and whether enough
available_after_dilution = vol_library_to_use + naoh_volume + tris_volume
st.markdown(f"**Available volume after dilution:** {available_after_dilution:.2f} µL. Planned final loading volume: {final_loading_volume:.2f} µL.")
if available_after_dilution < final_loading_volume:
    st.warning('Available diluted volume is smaller than your planned loading volume. Increase pool aliquot or Tris/diluent volume.')
else:
    st.success('Diluted volume meets planned loading volume.')

# Export results
st.header('Export')
buf = io.BytesIO()
out_df = df.copy()
out_df['Total_pooled_volume_uL'] = pool_volume
out_df.to_csv(buf, index=False)
buf.seek(0)
st.download_button('Download pool plan as CSV', data=buf, file_name='pool_plan.csv')

st.markdown('---')
st.write('If you want, download this app script and run with:')
st.code('pip install streamlit pandas\nstreamlit run denature_dilute_app.py')

st.caption('This calculator provides **estimates** and flexible parameters. Verify volumes and concentrations against your laboratory SOP and instrument documentation before proceeding with experiments.')
