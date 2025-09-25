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
    "1\t285\t116560\t11.00\n"
    "2\t285\t116560\t7.72\n"
    "3\t285\t116560\t6.36\n"
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
        df["Read Fraction"] = df["Frac of Cart (%)"] / df["Frac of Cart (%)"].sum()

        # Weighted average library size
        weighted_avg_size = (df["Library Size"] * df["Read Fraction"]).sum()

        # Mass Needed formula
        df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)

        # Volume Needed for pipetting (for user guidance)
        df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

        # Optional per-library dilution for pipette convenience
        raw_vols = df["Volume Needed (µL)"].fillna(0).astype(float)
        dilution_factors = []
        diluted_vols = []

        for raw in raw_vols:
            if raw <= 0:
                d = 1.0
                diluted = 0.0
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
        df["Diluted Vol (µL)"] = diluted_vols  # for user pipetting guidance

        # Cartridge Utilization Percentage
        total_reads_required = (df["Unique Oligos"].astype(float) * desired_coverage).sum()
        utilization_pct = (total_reads_required / cartridge_capacity) * 100
        st.markdown(f"**Cartridge Utilization:** {utilization_pct:.2f}% of {cartridge_capacity:,} reads")

        st.subheader("📊 Input and Calculations")
        st.dataframe(df)

        # --- Corrected Pool Concentration ---
        st.subheader("📌 Pool Concentration")

        # Use Mass Needed / sum of actual volumes to be combined
        total_mass_ng = df["Mass Needed (ng)"].sum()
        total_volume_uL = df["Volume Needed (µL)"].sum()  # actual volumes added to pool

        pool_conc_ng_uL = total_mass_ng / total_volume_uL
        pool_conc_nM = pool_conc_ng_uL * 0.8 * 1e6 / (660 * weighted_avg_size)

        st.write(f"**Total mass pooled:** {total_mass_ng:.2f} ng")
        st.write(f"**Total volume pooled:** {total_volume_uL:.2f} µL")
        st.write(f"**Calculated pool concentration (ng/µL):** {pool_conc_ng_uL:.3f}")
        st.write(f"**Pooled library concentration (nM):** {pool_conc_nM:.2f} nM")

        # Determine if denature/dilute steps are required
        if pool_conc_nM > loading_conc:
            dilution_factor = pool_conc_nM / loading_conc
            denature_vol = total_volume_uL / dilution_factor
            instructions_pool = f"Dilute the pool: **{denature_vol:.1f} µL** pooled libraries."
        else:
            instructions_pool = f"No dilution needed; combine **{total_volume_uL:.1f} µL** pool directly."

        # PhiX addition
        if include_phix:
            phix_vol = round(total_volume_uL * (phix_dilution / 40), 1)
        else:
            phix_vol = 0.0

        total_mix = total_volume_uL + phix_vol

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
