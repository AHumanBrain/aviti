import streamlit as st
import pandas as pd
import io

st.set_page_config(layout="wide")
st.title("Library Pooling, Dilution, and Loading Calculator")

st.markdown("""
Paste **tab-separated values (TSV)** below with the following columns only:

**Alias | Library Size (bp) | Unique Oligos | Qubit Quant (ng/µL)**

# --------------------
# Example TSV
# --------------------
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
    height=300,
    help="Copy from Google Sheets and paste here (tab-separated)."
)

# --------------------
# Global controls
# --------------------
cartridge_capacity = st.selectbox(
    "Cartridge Capacity (reads)",
    options=[100_000_000, 500_000_000, 1_000_000_000],
    index=0,
    format_func=lambda x: f"{x:,}"
)

desired_coverage = st.number_input(
    "Desired Coverage (applied to all libraries)",
    min_value=1,
    max_value=1000,
    value=40,
    step=1
)

loading_conc_pM = st.number_input("Target loading concentration (pM)", value=10.0, step=0.5)

# PhiX options
include_phix = st.checkbox("Include PhiX spike-in", value=True)
phix_input_type = st.radio("PhiX stock type", ["1 nM stock", "Dilution factor"], horizontal=True)
if phix_input_type == "Dilution factor":
    phix_dilution = st.number_input("Enter PhiX dilution factor (e.g., 40 for 1:40)", value=40, step=1)
else:
    phix_dilution = 1

phiX_pct = st.number_input("Desired PhiX (molar %) in pre-denature mix", min_value=0.0, max_value=50.0, value=10.0, step=0.5)
final_volume_uL = st.number_input("Final volume after neutralization/dilution (µL)", value=1400.0, step=10.0)

# --------------------
# Parse TSV and calculations
# --------------------
if not txt.strip():
    st.info("Paste TSV values to get started.")
else:
    try:
        df = pd.read_csv(io.StringIO(txt.strip()), sep="\t", header=None)
        if df.shape[1] != 4:
            st.error("TSV must have exactly 4 columns: Alias, Library Size (bp), Unique Oligos, Qubit Quant (ng/µL)")
        else:
            df.columns = ["Alias", "Library Size", "Unique Oligos", "Qubit Quant (ng/µL)"]

            # Fraction of reads for weighted average
            df["Read Weight"] = df["Unique Oligos"] * desired_coverage
            df["Read Fraction"] = df["Read Weight"] / df["Read Weight"].sum()

            # Weighted average library size (bp)
            weighted_avg_size = (df["Library Size"] * df["Read Fraction"]).sum()

            # Frac of cartridge (%) & Mass Needed (ng)
            df["Frac of Cart (%)"] = ((df["Unique Oligos"] * desired_coverage) / cartridge_capacity * 100).round(6)
            df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)

            # Raw volume required (µL)
            df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

            # Per-library dilution (pipette-friendly)
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
            df["Diluted Vol (µL)"] = diluted_vols

            # --- Display Input & Calculations table ---
            display_cols = [
                "Alias", "Library Size", "Unique Oligos", "Qubit Quant (ng/µL)",
                "Frac of Cart (%)", "Mass Needed (ng)", "Volume Needed (µL)",
                "Dilution Factor", "Diluted Vol (µL)"
            ]
            st.subheader("📊 Input and Calculations")
            st.dataframe(df[display_cols].style.format({
                "Qubit Quant (ng/µL)": "{:.2f}",
                "Frac of Cart (%)": "{:.6f}",
                "Mass Needed (ng)": "{:.6f}",
                "Volume Needed (µL)": "{:.6f}",
                "Dilution Factor": "{:.2f}",
                "Diluted Vol (µL)": "{:.2f}"
            }))

            # --- Pool concentration ---
            total_mass_ng = df["Mass Needed (ng)"].sum()
            total_pooled_volume_uL = df["Diluted Vol (µL)"].sum()
            if total_pooled_volume_uL <= 0:
                st.error("Total pooled volume is zero — check your inputs/Qubit concentrations.")
            else:
                pool_conc_ng_uL = total_mass_ng / total_pooled_volume_uL
                pool_conc_nM = pool_conc_ng_uL * 0.8 * 1e6 / (660 * weighted_avg_size)

                st.subheader("📌 Pool Concentration (summary)")
                st.write(f"**Total mass pooled:** {total_mass_ng:.6f} ng")
                st.write(f"**Total pooled volume:** {total_pooled_volume_uL:.2f} µL")
                st.write(f"**Calculated pool concentration (ng/µL):** {pool_conc_ng_uL:.4f}")
                st.write(f"**Calculated pool concentration (nM):** {pool_conc_nM:.3f}")
                st.write(f"**Weighted average library size:** {weighted_avg_size:.1f} bp")

                measured_pool_conc = st.number_input(
                    "Measured pool concentration (ng/µL) — optional override",
                    min_value=0.0,
                    value=float(round(pool_conc_ng_uL, 4)),
                    step=0.01
                )

                # Compute pool + PhiX volumes
                pool_conc_nM_measured = measured_pool_conc * 0.8 * 1e6 / (660 * weighted_avg_size)
                pool_conc_pM_measured = pool_conc_nM_measured * 1000.0
                lib_target_pM = loading_conc_pM * (100 - phiX_pct) / 100.0
                phix_target_pM = loading_conc_pM * phiX_pct / 100.0
                V_pool_uL = lib_target_pM * final_volume_uL / pool_conc_pM_measured
                V_phix_uL = phix_target_pM * final_volume_uL / (1000 / phix_dilution if phix_input_type == "1 nM stock" else 1)
                total_mix_uL = V_pool_uL + V_phix_uL

                st.subheader("🔢 Computed mixing volumes")
                st.write(f"**Volume of pooled library (µL):** {V_pool_uL:.2f}")
                st.write(f"**PhiX volume (µL):** {V_phix_uL:.2f} (for {phiX_pct:.1f}% spike-in)")
                st.write(f"**Total pre-denature mix volume (µL):** {total_mix_uL:.2f}")

                # --- Step-by-step instructions ---
                try:
                    st.subheader("🧪 Step-by-step (high-level / follow your lab SOP)")
                    pool_phix_mix_uL = total_mix_uL
                    naoh_vol_uL = pool_phix_mix_uL
                    neutralize_vol_uL = pool_phix_mix_uL
                    buffer_vol_uL = final_volume_uL - (pool_phix_mix_uL + naoh_vol_uL + neutralize_vol_uL)
                    if buffer_vol_uL < 0:
                        buffer_vol_uL = 0.0
                        st.warning("Computed loading buffer volume < 0 µL. Check PhiX fraction, final volume, or measured pool concentration.")

                    instructions_md = f"""1. **Prepare individual libraries**: pipette each library at the Diluted Vol (µL) listed above.  
- Total pooled volume: **{total_pooled_volume_uL:.2f} µL**  
- Total mass pooled: **{total_mass_ng:.6f} ng**

2. **Combine pooled libraries** into a single tube. Mix gently.

3. **Measure pooled concentration** (Qubit or equivalent) and update if needed.

4. **Mix pool + PhiX**: transfer **{V_pool_uL:.2f} µL** pooled library and **{V_phix_uL:.2f} µL** PhiX.  
- Achieves ~{phiX_pct:.1f}% PhiX.  
- Total pool+PhiX mixture: **{pool_phix_mix_uL:.2f} µL**

5. **Denature with NaOH**:  
- Add **{naoh_vol_uL:.2f} µL** 0.2 N NaOH, mix gently, spin, incubate 5 min.  
- Add **{neutralize_vol_uL:.2f} µL** of Tris-HCl (pH 7.0), mix, spin.

6. **Add loading buffer** to bring total volume to **{final_volume_uL:.0f} µL**:  
- Volume of buffer needed: **{buffer_vol_uL:.2f} µL**

7. **Load** all **{final_volume_uL:.0f} µL** into the cartridge according to SOP.

> Note: If required pool volume exceeds available pooled volume, prepare more pool or adjust plan.
"""
                    st.markdown(instructions_md)
                except Exception as e:
                    st.error(f"Error generating step-by-step instructions: {e}")

    except Exception as e:
        st.error(f"Error parsing TSV: {e}")
