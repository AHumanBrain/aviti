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
    "4\t285\t116560\t6.25\n"
    "5\t285\t116560\t4.56\n"
    "6\t285\t116560\t7.00\n"
    "7\t245\t19125\t6.13\n"
    "8\t285\t18097\t6.41\n"
    "9\t285\t3857\t11.50\n"
    "10\t285\t79937\t5.86\n"
    "11\t223\t5460\t20.80\n"
    "12\t223\t17370\t19.40\n"
    "13\t223\t673\t15.60\n"
    "14\t324\t100000\t30.00\n"
)

txt = st.text_area(
    "Paste TSV here",
    value=placeholder_text,
    height=200,
    help="Copy from Google Sheets and paste here (tab-separated)."
)

cartridge_capacity = st.selectbox(
    "Cartridge Capacity (reads)",
    options=[100_000_000, 500_000_000, 1_000_000_000],
    index=1,
    format_func=lambda x: f"{x:,}"
)

desired_coverage = st.number_input(
    "Desired Coverage (applied to all libraries)",
    min_value=1,
    max_value=1000,
    value=40,
    step=1
)

loading_conc = st.number_input("Desired loading concentration (pM)", value=10.0, step=0.5)
include_phix = st.checkbox("Include PhiX spike-in", value=True)
phiX_pct = 10.0  # default 10% spike-in

if include_phix:
    phix_input_type = st.radio("PhiX stock type", ["1 nM stock", "Dilution factor"], horizontal=True)
    if phix_input_type == "Dilution factor":
        phix_dilution = st.number_input("Enter PhiX dilution factor (e.g., 40 for 1:40)", value=40, step=1)
    else:
        phix_dilution = 1
else:
    phix_input_type = "1 nM stock"
    phix_dilution = 1

final_volume_uL = 1400.0

if txt.strip():
    try:
        df = pd.read_csv(io.StringIO(txt), sep="\t", header=None)
        df.columns = ["Alias", "Library Size", "Unique Oligos", "Qubit Quant (ng/µL)"]

        # Qubit ng/µL → nM
        df["Qubit Conc (nM)"] = (df["Qubit Quant (ng/µL)"] * 1e6) / (660 * df["Library Size"])

        # Fraction of cartridge
        df["Frac of Cart (%)"] = ((df["Unique Oligos"] * desired_coverage) / cartridge_capacity * 100).round(3)

        # Mass needed
        df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)
        df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

        # Cartridge utilization
        total_reads_required = (df["Unique Oligos"].astype(float) * desired_coverage).sum()
        utilization_pct = (total_reads_required / cartridge_capacity) * 100
        st.markdown(f"**Cartridge Utilization:** {utilization_pct:.2f}% of {cartridge_capacity:,} reads")

        # Per-library dilution factor
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

        # Compute pool concentration internally (preserves original float precision)
        total_mass_ng = df["Mass Needed (ng)"].sum()
        total_pooled_volume_uL = df["Diluted Vol (µL)"].sum()
        calculated_pool_conc_ng_uL = total_mass_ng / total_pooled_volume_uL if total_pooled_volume_uL > 0 else 0.0

        measured_pool_conc_ng_uL_temp = calculated_pool_conc_ng_uL  # internal precise value

        # Weighted average library size
        weighted_lib_size = (df["Library Size"] * df["Mass Needed (ng)"]).sum() / df["Mass Needed (ng)"].sum()

        # Pool concentration in nM
        pool_conc_nM_measured = measured_pool_conc_ng_uL_temp * 0.8 * 1e6 / (660 * weighted_lib_size)

        st.write(f"**Calculated pool concentration (ng/µL):** {calculated_pool_conc_ng_uL:.3f}")
        st.write(f"**Measured pool concentration (nM):** {pool_conc_nM_measured:.3f}")

        # Show number_input for user only (does not affect calculation unless user changes)
        measured_pool_conc_ng_uL = st.number_input(
            "Measured pooled library concentration (ng/µL)",
            value=calculated_pool_conc_ng_uL,
            step=0.01
        )

        # Pool + PhiX volumes (exact formulas)
        V_pool_uL = loading_conc * (100 - phiX_pct) / 100 * final_volume_uL / (pool_conc_nM_measured * 1000)
        if phix_input_type == "1 nM stock":
            phix_stock_pM = 1000 / phix_dilution
        else:
            phix_stock_pM = 1 / phix_dilution
        V_phix_uL = loading_conc * phiX_pct / 100 * final_volume_uL / phix_stock_pM
        total_mix_uL = V_pool_uL + V_phix_uL

        # Step-by-step instructions
        try:
            st.subheader("🧪 Step-by-step (high-level / follow your lab SOP)")

            pool_phix_mix_uL = total_mix_uL
            naoh_vol_uL = pool_phix_mix_uL
            neutralize_vol_uL = pool_phix_mix_uL
            buffer_vol_uL = final_volume_uL - (pool_phix_mix_uL + naoh_vol_uL + neutralize_vol_uL)
            if buffer_vol_uL < 0:
                buffer_vol_uL = 0.0
                st.warning("Computed loading buffer volume < 0 µL. Check PhiX fraction, final volume, or measured pool concentration.")

            instructions_md = f"""
1. **Prepare individual libraries**: pipette each library at the **Diluted Vol (µL)** listed above.  
   - Total pooled volume: **{total_pooled_volume_uL:.2f} µL**  
   - Total mass pooled: **{total_mass_ng:.6f} ng**

2. **Combine pooled libraries** into a single tube. Mix gently.

3. **Measure pooled concentration** (Qubit or equivalent) and update the "Measured pool concentration" if different.

4. **Mix pool + PhiX**: transfer **{V_pool_uL:.2f} µL** of the pooled library and **{V_phix_uL:.2f} µL** PhiX (if using) into a clean tube.  
   - Total pool+PhiX mixture: **{pool_phix_mix_uL:.2f} µL**

5. **Denature with NaOH**:  
   - Add **{naoh_vol_uL:.2f} µL** of 0.2 N NaOH (equal to pool+PhiX volume), mix gently, spin briefly, and incubate at room temperature for 5 minutes.  
   - Add **{neutralize_vol_uL:.2f} µL** of pool+PhiX mixture to neutralize.

6. **Add loading buffer** to bring total volume to **{final_volume_uL:.0f} µL**:  
   - Volume of buffer needed: **{buffer_vol_uL:.2f} µL**  

7. **Load** all **{final_volume_uL:.0f} µL** into the cartridge according to your sequencer’s SOP.
"""
            st.markdown(instructions_md)

        except Exception as e:
            st.error(f"Error generating step-by-step instructions: {e}")

    except Exception as e:
        st.error(f"Error parsing TSV or computing values: {e}")
