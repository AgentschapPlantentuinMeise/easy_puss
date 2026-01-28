"""easy_puss: Easy Package for Understanding of Spectral Signals"""

import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

class CompositionAnalysis:
    pass

class ICPOES(CompositionAnalysis):
    def __init__(
            self, file, type='agilent',
            min_corr_coeff: float =0.995,
            max_conc_rsd: float =5,
            max_bias: float =10,
            use_q_flag: bool =True,
            use_rse: bool =True,
            min_std_count: int =10,
            decimals=2, flag=True
    ):
        if type == 'agilent':
            self.data = pd.read_csv(file, skiprows=6, low_memory=False)
            self.process_agilent(
                min_corr_coeff, max_conc_rsd, max_bias, use_q_flag, use_rse,
                min_std_count
            )
        else:
            raise NotImplementedError
        self.masked_results = self.apply_lod_loq_mask(decimals=decimals, flag=flag)
        
    def process_agilent(self, min_corr_coeff, max_conc_rsd, max_bias, use_q_flag, use_rse,min_std_count):
        self.results, self.decisions, self.wl_scores, self.qc_summary = self.select_wavelengths(
            min_corr_coeff=min_corr_coeff, max_conc_rsd=max_conc_rsd,
            max_bias=max_bias, use_q_flag=use_q_flag, use_rse=use_rse
        )
        self.report = self.diagnostic_report()
        self.lodq = self.calibrate_and_compute_lod(min_std_count=min_std_count)
    
    def select_wavelengths(
        self,
        min_corr_coeff: float =0.995,
        max_conc_rsd: float =5,
        max_bias: float =10,
        use_q_flag: bool =True,
        use_rse: bool =True
    ):
        """
        Production‑ready wavelength selection for Agilent 5800 ICP‑OES.
        Uses:
          - Element Symbol
          - Wavelength
          - Radial flag
        """
        df = self.data
        
        # ---------------------------------------------------------
        # 0. SPLIT ELEMENT COLUMN INTO SYMBOL + WAVELENGTH
        # ---------------------------------------------------------
        # Example: "Ag 224.641_R"
        df[["Element Symbol", "Wavelength"]] = (
            df["Element"]
            .astype(str)
            .str.split(n=1, expand=True)
        )
    
        # ---------------------------------------------------------
        # 1. RADIAL FLAG + CLEAN WAVELENGTH
        # ---------------------------------------------------------
        df["Radial"] = df["Element Label"].str.endswith("_R")
    
        # ---------------------------------------------------------
        # 2. ULTRA‑ROBUST NUMERIC CLEANING
        # ---------------------------------------------------------
        def clean_numeric(series):
            cleaned = (
                series.astype(str)
                      .str.replace(",", ".", regex=False)
                      .str.replace(r"[A-Za-z#]", "", regex=True)
                      .str.replace("<|>|—|–|-", "", regex=True)
                      .str.strip()
            )
    
            # Replace empty strings with NaN without warnings
            cleaned = cleaned.mask(cleaned == "", np.nan)
    
            return pd.to_numeric(cleaned, errors="coerce")
    
    
        numeric_cols = [
            "Concentration",
            "Concentration % RSD",
            "Intensity % RSD",
            "Correlation coefficient",
            "%RSE",
            "%RSE limit",
        ]
    
        for col in numeric_cols:
            df[col] = clean_numeric(df[col])
    
        # ---------------------------------------------------------
        # 3. BASIC FLAGS & QC DETECTION
        # ---------------------------------------------------------
        df["_is_crm"] = df["Type"] == "QC"
        df["_has_q_flag"] = df["Flags"].astype(str).str.contains("Q", na=False)
    
        df_crm = df[df["_is_crm"]]
        df_unk = df[~df["_is_crm"]]
    
        # ---------------------------------------------------------
        # 4. CRM BIAS CALCULATION
        # ---------------------------------------------------------
        crm_bias = None
        if not df_crm.empty:
            crm_bias = (
                df_crm.groupby(["Element Symbol", "Wavelength", "Radial"])["Concentration"]
                      .apply(lambda x: x.mean())
                      .rename("crm_mean")
            )
    
        # ---------------------------------------------------------
        # 5. WAVELENGTH SCORING
        # ---------------------------------------------------------
        def score_wavelength(group):
            g_unk = group[~group["_is_crm"]]
            if g_unk.empty:
                g_unk = group
    
            conc_rsd = g_unk["Concentration % RSD"].median()
            intensity_rsd = g_unk["Intensity % RSD"].median()
            corr_coeff = group["Correlation coefficient"].max()
            rse = group["%RSE"].median()
            rse_limit = group["%RSE limit"].median()
            q_flag = group["_has_q_flag"].any()
    
            key = group.name
    
            bias = crm_bias.loc[key] if crm_bias is not None and key in crm_bias.index else np.nan
    
            all_nan_conc = group["Concentration"].isna().all()
    
            reject_reasons = []
            if all_nan_conc:
                reject_reasons.append("no numeric data")
    
            if use_q_flag and q_flag:
                reject_reasons.append("Q flag")
    
            if corr_coeff < min_corr_coeff:
                reject_reasons.append("low corr coeff")
    
            if use_rse and rse > rse_limit:
                reject_reasons.append("RSE > limit")
            
            if conc_rsd > max_conc_rsd:
                reject_reasons.append("high %RSD")
    
            # TODO improve CRM bias calculation logic
            #if use_q_flag and not np.isnan(bias) and abs(bias) > max_bias:
            #    reject_reasons.append("high CRM bias")
    
            rejected = len(reject_reasons) > 0
    
            score = 0
            if not rejected:
                score += max(0, 10 - conc_rsd)
                score += (corr_coeff - 0.995) * 100
                if not np.isnan(bias):
                    score += max(0, 10 - abs(bias))
    
            return pd.Series({
                "wavelength": key[1],
                "radial": key[2],
                "rejected": rejected,
                "reject_reasons": "; ".join(reject_reasons),
                "score": score,
                "conc_rsd": conc_rsd,
                "intensity_rsd": intensity_rsd,
                "corr_coeff": corr_coeff,
                "rse": rse,
                "rse_limit": rse_limit,
                "crm_bias": bias
            })
    
        wl_scores = (
            df.groupby(["Element Symbol", "Wavelength", "Radial"])
              .apply(score_wavelength, include_groups=False)
              .reset_index()
        )
    
        # ---------------------------------------------------------
        # 6. FINAL SELECTION PER SAMPLE
        # ---------------------------------------------------------
        results = []
        group_cols = ["Label", "Rack:Tube", "Element Symbol"]
    
        for key, g in df_unk.groupby(group_cols):
            element = key[2]
    
            g = g.merge(
                wl_scores[wl_scores["Element Symbol"] == element][
                    ["Wavelength", "Radial", "rejected", "score"]
                ],
                on=["Wavelength", "Radial"],
                how="left"
            )
    
            valid = g[~g["rejected"]]
            if valid.empty:
                results.append({
                    "Label": key[0],
                    "Rack:Tube": key[1],
                    "Element Symbol": element,
                    "final_concentration": np.nan,
                    "decision": "all_rejected",
                    "used_wavelengths": ""
                })
                continue
    
            wl_conc = (
                valid.groupby(["Wavelength", "Radial"])["Concentration"]
                     .apply(lambda x: x.mean())
                     .rename("mean_conc")
                     .reset_index()
            )
    
            wl_conc = wl_conc.merge(
                wl_scores[wl_scores["Element Symbol"] == element][
                    ["Wavelength", "Radial", "score"]
                ],
                on=["Wavelength", "Radial"],
                how="left"
            ).sort_values("score", ascending=False)
    
            top = wl_conc.iloc[0]
    
            similar = [top]
            for _, row in wl_conc.iloc[1:].iterrows():
                if np.isnan(top["mean_conc"]) or np.isnan(row["mean_conc"]) or top["mean_conc"] == 0:
                    continue
                diff = abs(row["mean_conc"] - top["mean_conc"]) / top["mean_conc"] * 100
                if diff <= 10:
                    similar.append(row)
    
            if len(similar) > 1:
                final = np.mean([r["mean_conc"] for r in similar])
                used = [f"{r['Wavelength']}{'_R' if r['Radial'] else ''}" for r in similar]
                decision = "average"
            else:
                final = top["mean_conc"]
                used = [f"{top['Wavelength']}{'_R' if top['Radial'] else ''}"]
                decision = "single"
    
            results.append({
                "Label": g["Label"].iloc[0],
                "Rack:Tube": g["Rack:Tube"].iloc[0],
                "Element Symbol": element,
                "final_concentration": final,
                "decision": decision,
                "used_wavelengths": ", ".join(used)
            })
    
        decisions_df = pd.DataFrame(results)
        final_df = (
            decisions_df
                .pivot_table( 
                    index=["Label", "Rack:Tube"],
                    columns="Element Symbol", 
                    values="final_concentration" 
                )
        )
        # ---------------------------------------------------------
        # 7. QC SUMMARY TABLE
        # ---------------------------------------------------------
        qc_summary = wl_scores.copy()
        qc_summary["status"] = np.where(qc_summary["rejected"], "Rejected", "Accepted")
    
        return final_df, decisions_df, wl_scores, qc_summary
    
    # Recalibration and measuring LOD + LOQ
    import numpy as np
    import pandas as pd
    
    def calibrate_and_compute_lod(self, wl_scores=None, min_std_count=10):
        """
        Adds calibration slope, intercept, R², LOD, and LOQ to wl_scores.
        Calibration is linear with intercept: I = m*C + b.
        LOD and LOQ are computed from blank intensity noise.
        """
        df = self.data
        wl_scores = self.wl_scores if wl_scores is None else wl_scores
        
        # ---------------------------------------------------------
        # 1. Identify valid calibration levels (expected concentrations)
        # ---------------------------------------------------------
        std_counts = (
            df[df["Type"] == "STD"]["Concentration"]
            .value_counts()
        )
        known_levels = std_counts[std_counts >= min_std_count].index.to_list()
    
        # ---------------------------------------------------------
        # 2. Helper: build calibration model
        # ---------------------------------------------------------
        def build_calibration(cal):
            conc = cal["Concentration"].astype(float)
            inten = cal["Intensity"].astype(float)
    
            # Fit linear regression with intercept
            m, b = np.polyfit(conc, inten, 1)
    
            # Compute R²
            pred = m * conc + b
            ss_res = np.sum((inten - pred)**2)
            ss_tot = np.sum((inten - np.mean(inten))**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    
            return m, b, r2
    
        # ---------------------------------------------------------
        # 3. Helper: compute LOD and LOQ from blank intensities
        # ---------------------------------------------------------
        def compute_lod_loq(blank_intensities, slope):
            if slope == 0 or np.isnan(slope) or len(blank_intensities) < 2:
                return np.nan, np.nan
    
            sd_blank = np.std(blank_intensities, ddof=1)
            lod = 3 * sd_blank / slope
            loq = 10 * sd_blank / slope
            return lod, loq
    
        # ---------------------------------------------------------
        # 4. Loop through each wavelength and compute calibration + LOD
        # ---------------------------------------------------------
        slopes = []
        intercepts = []
        r2s = []
        lods = []
        loqs = []
    
        for idx, row in wl_scores.iterrows():
            elem = row["Element Symbol"]
            wl = row["Wavelength"]
            rad = row["Radial"]
    
            # Extract calibration standards for this wavelength
            cal = df[
                (df["Element Symbol"] == elem) &
                (df["Wavelength"] == wl) &
                (df["Radial"] == rad) &
                (df["Type"] == "STD")
            ]
    
            # Keep only valid calibration levels
            cal = cal[cal["Concentration"].isin(known_levels)]
    
            # Extract blanks
            blanks = df[
                (df["Element Symbol"] == elem) &
                (df["Wavelength"] == wl) &
                (df["Radial"] == rad) &
                (df["Type"] == "BLK")
            ]
    
            # Build calibration
            if len(cal) >= 2:
                m, b, r2 = build_calibration(cal)
            else:
                m, b, r2 = np.nan, np.nan, np.nan
    
            # Compute LOD/LOQ
            lod, loq = compute_lod_loq(blanks["Intensity"].dropna().astype(float), m)
    
            slopes.append(m)
            intercepts.append(b)
            r2s.append(r2)
            lods.append(lod)
            loqs.append(loq)
    
        # ---------------------------------------------------------
        # 5. Attach results to wl_scores
        # ---------------------------------------------------------
        wl_scores["slope"] = slopes
        wl_scores["intercept"] = intercepts
        wl_scores["r2"] = r2s
        wl_scores["LOD"] = lods
        wl_scores["LOQ"] = loqs
    
        lodq = wl_scores[~wl_scores.rejected].groupby('Element Symbol')[['LOD','LOQ']].mean()
        
        return lodq
    
    def apply_lod_loq_mask(self, results=None, lodq=None, decimals=2, flag=True):
        """
        Apply LOD/LOQ masking rules to a results dataframe.
    
        Parameters
        ----------
        results : pd.DataFrame
            Analytical results with element symbols as columns.
        lodq : pd.DataFrame
            DataFrame indexed by element symbol with columns 'LOD' and 'LOQ'.
    
        Returns
        -------
        pd.DataFrame
            Masked results where:
            - values < LOD → "< LOD"
            - LOD ≤ values < LOQ → "> LOD"
            - values ≥ LOQ → numeric value
        """
        results = self.results if results is None else results
        lodq = self.lodq if lodq is None else lodq
    
        # Align LOD/LOQ to results columns
        LOD = lodq["LOD"].reindex(results.columns)
        LOQ = lodq["LOQ"].reindex(results.columns)
    
        # Broadcast to match results shape
        lod_mat = pd.DataFrame([LOD] * len(results), index=results.index)
        loq_mat = pd.DataFrame([LOQ] * len(results), index=results.index)
    
        # Format LOD values to 2 decimals
        lod_fmt = lod_mat.round(decimals).astype(str)
        
        # Start with numeric values
        out = results.copy().round(decimals).astype(str)
    
        # Apply masks
        if flag:
            # Below LOD → * 
            out = out.mask(results < lod_mat, out + "*") 
            # Between LOD and LOQ → ** 
            out = out.mask((results >= lod_mat) & (results < loq_mat), out + "**")
        else:
            out = out.mask(results < lod_mat, "< " + lod_fmt)
            out = out.mask((results >= lod_mat) & (results < loq_mat),
                       "> " + lod_fmt)
    
        return out
    
    def diagnostic_report(self, wl_scores=None, results_df=None):
        df = self.data
        wl_scores = self.wl_scores if wl_scores is None else wl_scores
        results_df = self.results if results_df is None else results_df
        
        report = {}
    
        # 1. Elements in raw data
        raw_elements = set(df["Element Symbol"].unique())
        report["raw_elements"] = raw_elements
    
        # 2. Elements with at least one wavelength group
        wl_elements = set(wl_scores["Element Symbol"].unique())
        report["wl_elements"] = wl_elements
    
        # 3. Elements fully rejected at wavelength QC stage
        rejected_elements = set(
            wl_scores.groupby("Element Symbol")["rejected"].all()
            .loc[lambda s: s].index
        )
        report["rejected_elements"] = rejected_elements
    
        # 4. Elements with no unknown samples
        unk_elements = set(df.loc[~df["_is_crm"], "Element Symbol"].unique())
        no_unknowns = raw_elements - unk_elements
        report["no_unknowns"] = no_unknowns
    
        # 5. Elements that passed QC but still missing in final results
        final_elements = set(results_df.columns)
        passed_qc_missing = (wl_elements - rejected_elements) - final_elements
        report["passed_qc_missing"] = passed_qc_missing
    
        # 6. Elements in final results
        report["final_elements"] = final_elements
    
        # 7. Detailed per-element explanation
        details = {}
        for elem in sorted(raw_elements):
            details[elem] = {
                "in_raw_data": elem in raw_elements,
                "has_wavelengths": elem in wl_elements,
                "all_wavelengths_rejected": elem in rejected_elements,
                "has_unknowns": elem in unk_elements,
                "appears_in_final": elem in final_elements,
                "missing_reason": None
            }
    
            if elem not in wl_elements:
                details[elem]["missing_reason"] = "No wavelength groups found"
            elif elem in rejected_elements:
                details[elem]["missing_reason"] = "All wavelengths failed QC"
            elif elem not in unk_elements:
                details[elem]["missing_reason"] = "No unknown samples"
            elif elem not in final_elements:
                details[elem]["missing_reason"] = "No valid concentration after QC"
            else:
                details[elem]["missing_reason"] = "OK"
    
        report["details"] = details
    
        return report
    
