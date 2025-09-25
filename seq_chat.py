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

coverage = st.number_input(
    "Desired coverage (applies to all libraries)",
    value=40,
    step=5
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

        # Convert ng/µL to nM: (ng/µL * 10^6) / (660 g/mol/bp * bp length)
        df["Qubit Conc (nM)"] = (df["Qubit Quant (ng/µL)"] * 1e6) / (660 * df["Library Size"])

        # Fraction of cartridge (%)
        df["Frac of Cart (%)"] = (
            (df["Unique Oligos"] * coverage) / cartridge_capacity * 100
        ).round(3)

        # Mass and volume (ng, µL)
        df["Mass Needed (ng)"] = df["Unique Oligos"] * coverage / cartridge_capacity * df["Qubit Quant (ng/µL)"]
        df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

        # Dilution factor
        df["Dilution Factor"] = (df["Qubit Conc (nM)"] / loading_conc).apply(lambda x: max(1, round(x)))
        df["Diluted Vol (µL)"] = (df["Volume Needed (µL)"] * df["Dilution Factor"]).apply(lambda x: max(3, round(x, 1)))

        st.subheader("📊 Input and Calculations")
        st.dataframe(df)

        # Pool concentration (before PhiX)
        pool_conc = df["Qubit Conc (nM)"].mean()
        st.write(f"**Pooled library concentration:** {pool_conc:.2f} nM")

        # Pool dilution
        dilution_factor = pool_conc / (loading_conc * 0.99)
        diluted_pool_vol = max(4.5, round(30 / (1 + dilution_factor), 1))  # keep pipetteable

        # PhiX addition
        if include_phix:
            phix_vol = round(diluted_pool_vol * (phix_dilution / 40), 1)
        else:
            phix_vol = 0.0

        total_mix = diluted_pool_vol + phix_vol
        naoh_vol = total_mix
        tris_vol = total_mix
        final_vol = 1400
        buffer_vol = final_vol - (total_mix + naoh_vol + tris_vol)

        st.subheader("🧪 Step-by-step Instructions")
        st.markdown(f"""
        1. Pool the libraries according to calculated volumes above.  
        2. Dilute the pool: **{diluted_pool_vol} µL** pooled libraries.  
        3. Add PhiX: **{phix_vol} µL** ({'diluted PhiX' if phix_dilution > 1 else '1 nM PhiX stock'}).  
        4. Mix for a total of **{total_mix:.1f} µL**.  
        5. Add **{naoh_vol:.1f} µL** 0.2 N NaOH, mix, spin, and incubate 5 minutes at room temp.  
        6. Add **{tris_vol:.1f} µL** 0.2 M pH 7 Tris-HCl to neutralize.  
        7. Add **{buffer_vol:.1f} µL** loading buffer to bring total volume to {final_vol} µL.  
        8. Load all **{final_vol} µL** into the cartridge.  
        """)

    except Exception as e:
        st.error(f"Error parsing TSV: {e}")
