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

        # Convert ng/µL to nM for each library
        df["Qubit Conc (nM)"] = (df["Qubit Quant (ng/µL)"] * 1e6) / (660 * df["Library Size"])

        # Fraction of cartridge (%)
        df["Frac of Cart (%)"] = (
            (df["Unique Oligos"] * desired_coverage) / cartridge_capacity * 100
        ).round(3)

        # Mass and volume (ng, µL) - updated formula
        df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * ((df["Frac of Cart (%)"]) / 100)
        df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

        # --- Cartridge Utilization Percentage ---
        total_reads_required = (df["Unique Oligos"].astype(float) * desired_coverage).sum()
        utilization_pct = (total_reads_required / cartridge_capacity) * 100
        st.markdown(f"**Cartridge Utilization:** {utilization_pct:.2f}% of {cartridge_capacity:,} reads")

        # --- Per-library dilution factor calculation ---
        raw_vols = df["Volume Needed (µL)"].fillna(0).astype(float)
        dilution_factors = []
        diluted_vols = []

        for raw in raw_vols:
            if raw <= 0:
                d = 1.00
                diluted = 0.00
            else:
                # calculate dilution factor ensuring diluted volume between 1–10 µL
                d_min = 1.0 / raw
                d_max = 10.0 / raw
                d = max(d_min, 1.0)  # must be at least 1
                if d > d_max:
                    d = d_max
                d = round(d, 2)
                diluted = round(raw * d, 2)
            dilution_factors.append(d)
            diluted_vols.append(diluted)

        df["Dilution Factor"] = dilution_factors
        df["Diluted Vol (µL)"] = diluted_vols

        # --- Pool concentration calculation (fixed) ---
        total_mass_ng = df["Mass Needed (ng)"].sum()
        total_volume_uL = df["Volume Needed (µL)"].sum()   # use actual pipetted volumes

        pool_conc_ng_uL = total_mass_ng / total_volume_uL

        # Weighted average library size
        weighted_avg_size = (
            (df["Library Size"] * (df["Unique Oligos"] * desired_coverage)).sum()
            / (df["Unique Oligos"] * desired_coverage).sum()
        )

        # Convert to nM
        pool_conc_nM = pool_conc_ng_uL * 0.8 * 1e6 / (660 * weighted_avg_size)

        st.subheader("📊 Input and Calculations")
        st.dataframe(df)

        # --- Pool concentration display ---
        st.subheader("📌 Pool Concentration")
        st.write(f"**Total pooled volume:** {total_volume_uL:.2f} µL")
        st.write(f"**Calculated pool concentration (ng/µL):** {pool_conc_ng_uL:.3f} ng/µL")

        measured_pool_conc = st.number_input(
            "Measured pool concentration (ng/µL)",
            min_value=0.0,
            value=float(pool_conc_ng_uL),
            step=0.01
        )

        pool_conc_nM = measured_pool_conc * 0.8 * 1e6 / (660 * weighted_avg_size)
        st.write(f"**Pooled library concentration (nM) based on measured value:** {pool_conc_nM:.2f} nM")

        # --- Step-by-step Instructions ---
        st.subheader("🧪 Step-by-step Instructions")

        if include_phix:
            # Simple pooling: combine pooled library with PhiX directly
            # Example: 7.3 µL pool + 1.4 µL PhiX
            pool_vol_for_loading = 7.3  # you could compute this dynamically
            phix_vol = 1.4
            st.markdown(f"""
            1. Pool the libraries according to calculated volumes above.  
            2. Take **{pool_vol_for_loading:.1f} µL** pooled libraries (≈ {pool_conc_nM:.1f} nM).  
            3. Add **{phix_vol:.1f} µL** 1 nM PhiX.  
            4. Mix gently.  
            5. Load directly at **{loading_conc:.1f} pM** without further dilution.  
            """)
        else:
            st.markdown(f"""
            1. Pool the libraries according to calculated volumes above.  
            2. Use **{total_volume_uL:.1f} µL** pooled libraries at {pool_conc_nM:.2f} nM.  
            3. Dilute/adjust if necessary for {loading_conc:.1f} pM loading.  
            """)

    except Exception as e:
        st.error(f"Error parsing TSV: {e}")
