import streamlit as st
import pandas as pd
import io

st.title("Library Pooling, Dilution, and Loading Calculator")

st.markdown("""
Paste **tab-separated values (TSV)** below with the following columns only:  

**Alias | Library Size (bp) | Unique Oligos | Qubit Quant (ng/µL)**
""")

# Example TSV input
placeholder_text = (
    "sampleA\t300\t100000\t2.1\n"
    "sampleB\t350\t200000\t5.0\n"
)

txt = st.text_area(
    "Paste TSV here",
    value=placeholder_text,
    height=200,
    help="Copy from Google Sheets and paste here (tab-separated)."
)

# Global inputs
cartridge_capacity = st.selectbox(
    "Cartridge Capacity (reads)",
    options=[100_000_000, 500_000_000, 1_000_000_000],
    index=1,
    format_func=lambda x: f"{x:,}"
)

# Global desired coverage
desired_coverage = st.number_input(
    "Desired Coverage (applied to all libraries)",
    min_value=1,
    max_value=1000,
    value=40,
    step=1
)

loading_conc = st.number_input("Desired loading concentration (pM)", value=10.0, step=0.5)
include_phix = st.checkbox("Include PhiX spike-in", value=True)
phix_input_type = st.radio("PhiX stock type", ["1 nM stock", "Dilution factor"], horizontal=True)
if phix_input_type == "Dilution factor":
    phix_dilution = st.number_input("Enter PhiX dilution factor (e.g., 40 for 1:40)", value=40, step=1)
else:
    phix_dilution = 1

if txt.strip():
    try:
        df = pd.read_csv(io.StringIO(txt), sep="\t", header=None)
        df.columns = ["Alias", "Library Size", "Unique Oligos", "Qubit Quant (ng/µL)"]

        # Fraction of reads for each library (based on desired coverage)
        df["Frac of Cart (%)"] = ((df["Unique Oligos"] * desired_coverage) / cartridge_capacity * 100).round(3)
        df["Read Fraction"] = (df["Frac of Cart (%)"] / df["Frac of Cart (%)"].sum())

        # Weighted average library size
        weighted_avg_size = (df["Library Size"] * df["Read Fraction"]).sum()

        # Updated Mass Needed formula (line 71)
        df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)
        df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

        # Cartridge Utilization Percentage
        total_reads_required = (df["Unique Oligos"].astype(float) * desired_coverage).sum()
        utilization_pct = (total_reads_required / cartridge_capacity) * 100
        st.markdown(f"**Cartridge Utilization:** {utilization_pct:.2f}% of {cartridge_capacity:,} reads")

        # Per-library dilution factor calculation
        raw_vols = df["Volume Needed (µL)"].fillna(0).astype(float)
        dilution_factors = []
        diluted_vols = []

        for raw in raw_vols:
            if raw <= 0:
                d = 1.00
                diluted = 0.00
            else:
                d_min = 1.0 / raw
                d_max = 10.0 / raw
                d = max(d_min, 1.0)
                if d > d_max:
                    d = d_max
                d = round(d, 2)
                diluted = round(raw * d, 2)
            dilution_factors.append(d)
            diluted_vols.append(diluted)

        df["Dilution Factor"] = dilution_factors
        df["Diluted Vol (µL)"] = diluted_vols

        st.subheader("📊 Input and Calculations")
        st.dataframe(df)

        # --- Corrected Pool Concentration ---
        st.subheader("📌 Pool Concentration")

        # Calculated pool concentration (ng/µL) = sum of masses / sum of volumes
        calculated_pool_conc = df["Mass Needed (ng)"].sum() / df["Volume Needed (µL)"].sum()
        st.write(f"**Calculated pool concentration (ng/µL):** {calculated_pool_conc:.3f} ng/µL")

        # Measured pool concentration (ng/µL) input
        measured_pool_conc = st.number_input(
            "Measured pool concentration (ng/µL)",
            min_value=0.0,
            value=float(calculated_pool_conc),
            step=0.01
        )

        # Pooled concentration in nM using measured value
        pool_conc_nM = measured_pool_conc * 0.8 * 1e6 / (660 * weighted_avg_size)
        st.write(f"**Pooled library concentration (nM) based on measured value:** {pool_conc_nM:.2f} nM")

        # Determine if dilution is needed
        if pool_conc_nM > loading_conc:
            dilution_factor = pool_conc_nM / loading_conc
            diluted_pool_vol = round(df["Volume Needed (µL)"].sum() * dilution_factor, 1)
            instructions_pool = f"Dilute the pool: **{diluted_pool_vol} µL** pooled libraries."
        else:
            instructions_pool = f"No dilution needed; combine **{df['Volume Needed (µL)'].sum():.1f} µL** pool directly."

        # PhiX addition
        if include_phix:
            phix_vol = round(df["Volume Needed (µL)"].sum() * (phix_dilution / 40), 1)
        else:
            phix_vol = 0.0

        total_mix = df["Volume Needed (µL)"].sum() + phix_vol

        st.subheader("🧪 Step-by-step Instructions")
        st.markdown(f"""
1. Pool the libraries according to calculated volumes above.  
2. {instructions_pool}  
3. Add PhiX: **{phix_vol} µL** ({'diluted PhiX' if phix_dilution > 1 else '1 nM PhiX stock'}).  
4. Mix for a total of **{total_mix:.1f} µL**.  
5. Load into the cartridge.
""")

    except Exception as e:
        st.error(f"Error parsing TSV: {e}")
