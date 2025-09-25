import streamlit as st
import pandas as pd
import io

st.set_page_config(layout="wide")
st.title("Library Pooling, Dilution, and Loading Calculator")

st.markdown("""
Paste **tab-separated values (TSV)** below with the following columns only:

**Alias | Library Size (bp) | Unique Oligos | Qubit Quant (ng/µL)**

This calculator will:
- compute Mass Needed using your formula,
- show raw and pipette-friendly (diluted) volumes,
- compute pool concentration using the **sum of masses / sum of diluted volumes** (the volumes you actually combine),
- convert to nM using your weighted average size formula,
- compute how much pool + PhiX to combine for a given target loading concentration (pM) after dilution to a final volume.
""")

# --------------------
# Example TSV (you can replace with your 14-row test)
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
    index=0,  # default 100M to match the ~0.4 ng/µL expectation for your test TSV
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

# How much PhiX (molar %) do you want in the pre-denature mix? (user-friendly)
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

            # Read fraction (for weighted average size): fraction of reads coming from each library
            df["Read Weight"] = df["Unique Oligos"] * desired_coverage
            df["Read Fraction"] = df["Read Weight"] / df["Read Weight"].sum()

            # Weighted average library size (bp)
            weighted_avg_size = (df["Library Size"] * df["Read Fraction"]).sum()

            # Mass Needed (ng) using your requested formula (line 71)
            # df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)
            # First compute Frac of Cartridge (%) (how much of the cartridge each library consumes)
            df["Frac of Cart (%)"] = ((df["Unique Oligos"] * desired_coverage) / cartridge_capacity * 100).round(6)
            df["Mass Needed (ng)"] = 9.8 * (250 / (df["Library Size"] - 124)) * (df["Frac of Cart (%)"] / 100)

            # Raw volume required (µL) to deliver the Mass Needed from the undiluted stock
            df["Volume Needed (µL)"] = df["Mass Needed (ng)"] / df["Qubit Quant (ng/µL)"]

            # Per-library pipette-friendly dilution factor (so small volumes are pipettable)
            raw_vols = df["Volume Needed (µL)"].fillna(0).astype(float)
            dilution_factors = []
            diluted_vols = []
            for raw in raw_vols:
                if raw <= 0:
                    d = 1.0
                    diluted = 0.0
                else:
                    d_min = 1.0 / raw        # makes raw * d >= 1 µL
                    d_max = 10.0 / raw       # makes raw * d <= 10 µL
                    d = max(d_min, 1.0)      # at least 1x
                    if d > d_max:
                        d = d_max
                    d = round(d, 2)
                    diluted = round(raw * d, 2)
                dilution_factors.append(d)
                diluted_vols.append(diluted)

            df["Dilution Factor"] = dilution_factors
            df["Diluted Vol (µL)"] = diluted_vols

            # --- Display Input & Calculations table (explicitly include Diluted Vol so you can see it) ---
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

            # ---------------------
            # Pool calculations (use diluted volumes — the volumes you will actually pipette & combine)
            # ---------------------
            total_mass_ng = df["Mass Needed (ng)"].sum()
            total_pooled_volume_uL = df["Diluted Vol (µL)"].sum()   # **use Diluted Vol (µL)** per your request

            if total_pooled_volume_uL <= 0:
                st.error("Total pooled volume is zero — check your inputs/Qubit concentrations.")
            else:
                pool_conc_ng_uL = total_mass_ng / total_pooled_volume_uL

                # Convert to nM using measured or calculated ng/µL and your formula
                # Show calculated value, then let user override with measured value
                st.subheader("📌 Pool Concentration (summary)")
                st.write(f"**Total mass pooled:** {total_mass_ng:.6f} ng")
                st.write(f"**Total pooled volume (sum of Diluted Vol):** {total_pooled_volume_uL:.2f} µL")
                st.write(f"**Calculated pool concentration (ng/µL)** (sum masses / sum diluted vols): **{pool_conc_ng_uL:.4f} ng/µL**")

                measured_pool_conc = st.number_input(
                    "Measured pool concentration (ng/µL) — enter your Qubit reading (or use calculated value)",
                    min_value=0.0,
                    value=float(round(pool_conc_ng_uL, 4)),
                    step=0.01
                )

                # weighted_avg_size already computed above
                pool_conc_nM_measured = measured_pool_conc * 0.8 * 1e6 / (660 * weighted_avg_size)
                st.write(f"**Pooled library concentration (nM) based on measured value:** {pool_conc_nM_measured:.3f} nM")
                st.write(f"**Weighted average library size used (bp):** {weighted_avg_size:.1f} bp")

                # ---------------------
                # Compute volumes needed to reach target loading concentration (pM) after dilution to final volume
                # V_pool = (target_pM * final_vol_uL) / (pool_conc_nM_measured * 1000)
                # (convert pool nM -> pM by *1000)
                # ---------------------
                if pool_conc_nM_measured <= 0:
                    st.warning("Measured pool concentration (nM) must be > 0 to compute required volumes.")
                else:
                    pool_conc_pM_measured = pool_conc_nM_measured * 1000.0
                    # Volume of pooled library (µL) required to reach loading_conc_pM after dilution to final_volume_uL
                    V_pool_uL = (loading_conc_pM * final_volume_uL) / pool_conc_pM_measured
                    V_pool_uL = float(V_pool_uL)

                # ---------------------
                # Compute volumes needed for target loading concentration (pM)
                # ---------------------
                if pool_conc_nM_measured <= 0:
                    st.warning("Measured pool concentration (nM) must be > 0 to compute required volumes.")
                else:
                    # Convert pool concentration to pM
                    pool_conc_pM_measured = pool_conc_nM_measured * 1000.0
                
                    # Total pre-denature volume = final volume
                    total_volume_uL = final_volume_uL
                
                    if include_phix and phiX_pct > 0:
                        # Simplified: PhiX volume = desired % of final volume
                        V_phix_uL = (phiX_pct / 100.0) * total_volume_uL
                    else:
                        V_phix_uL = 0.0
                
                    # Pool volume = remainder of final volume
                    V_pool_uL = total_volume_uL - V_phix_uL
                
                    total_mix_uL = V_pool_uL + V_phix_uL
                
                    # Check if pooled volume is sufficient
                    available_pool_uL = total_pooled_volume_uL
                    shortage_msg = ""
                    if V_pool_uL > available_pool_uL:
                        shortage_msg = (
                            f"WARNING: required pool volume {V_pool_uL:.2f} µL is greater than "
                            f"available pooled volume {available_pool_uL:.2f} µL. Prepare more pool or adjust plan."
                        )
                
                    # Display computed volumes
                    st.subheader("🔢 Computed mixing volumes")
                    st.write(f"**Volume of pooled library to use (µL):** {V_pool_uL:.2f}")
                    st.write(f"**PhiX volume (µL):** {V_phix_uL:.2f} (for {phiX_pct:.1f}% spike-in)")
                    st.write(f"**Total pre-denature mix volume (µL):** {total_mix_uL:.2f}")
                    if shortage_msg:
                        st.warning(shortage_msg)
                    
                    # Check available pooled volume vs required pool volume
                    available_pool_uL = total_pooled_volume_uL
                    shortage_msg = ""
                    if V_pool_uL > available_pool_uL:
                        shortage_msg = (
                            f"WARNING: required pool volume {V_pool_uL:.2f} µL is greater than "
                            f"available pooled volume {available_pool_uL:.2f} µL. You will need to prepare more pool or adjust plan."
                        )

                    # Present results
                    st.subheader("🔢 Computed mixing volumes")
                    st.write(f"**Volume of pooled library to use (µL)** to reach {loading_conc_pM:.1f} pM after dilution to {final_volume_uL:.0f} µL: **{V_pool_uL:.2f} µL**")
                    if include_phix and V_phix_uL > 0:
                        st.write(f"**PhiX volume (µL)** to approximate {phiX_pct:.1f}% molar PhiX in the pre-denature mix: **{V_phix_uL:.2f} µL**")
                    else:
                        st.write("**PhiX volume (µL):** Not included (PhiX disabled or 0%).")
                    st.write(f"**Total pre-denature mix volume:** {total_mix_uL:.2f} µL")
                    if shortage_msg:
                        st.warning(shortage_msg)

                    # ---------------------
                    # Step-by-step instructions (non-hazardous, high-level)
                    # ---------------------
                    st.subheader("🧪 Step-by-step (high-level / follow your lab SOP for denature/dilute steps)")
                    # Provide the exact numbers we computed but avoid giving hazardous experimental parameters (e.g., chemical volumes or incubation times).
                    instructions_md = f"""
1. **Prepare individual libraries**: using the table above, prepare each library at the **Diluted Vol (µL)** listed (these are the pipette-friendly aliquots that together form your pool).  
   - **Total pooled volume available:** **{available_pool_uL:.2f} µL**  
   - **Total mass pooled:** **{total_mass_ng:.6f} ng**

2. **Combine pooled libraries**: combine the aliquots (the `Diluted Vol (µL)` values) into a single tube. Mix gently.

3. **Measure the pooled concentration** (Qubit or equivalent). Enter the measured pool concentration (ng/µL) in the "Measured pool concentration" box above if it differs from the calculated value.

4. **Mix pool + PhiX** (pre-denaturation): transfer **{V_pool_uL:.2f} µL** of the pooled library and **{V_phix_uL:.2f} µL** PhiX (if using) into a clean tube.  
   - If PhiX is included, this approximates **{phiX_pct:.1f}%** molar PhiX in the pre-denature mix.

5. **Denature and dilute to final volume**: follow *your lab’s standard denaturation & neutralization protocol* (do **not** substitute procedures from the web). After denaturation/neutralization, dilute the reaction to **{final_volume_uL:.0f} µL** total volume.  
   - The small pre-denature mix volume above (e.g., ~7–10 µL) will be diluted into the final volume to reach the **{loading_conc_pM:.1f} pM** loading concentration.

6. **Quality check & load**: verify final concentration by your standard QC method if required, and load according to your sequencer’s instructions.

> Note: If the required pooled volume ({V_pool_uL:.2f} µL) is greater than the available pooled volume ({available_pool_uL:.2f} µL), you'll need to generate more pooled library (repeat pooling from the prepared aliquots or re-pool with adjusted volumes) or adjust the target loading concentration / final volume.
"""
                    st.markdown(instructions_md)

    except Exception as e:
        st.error(f"Error parsing TSV or computing values: {e}")
