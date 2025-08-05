import pandas as pd
import numpy as np
import logging
from scipy.optimize import minimize
from scipy.stats import betabinom, chi2
from statsmodels.stats.multitest import multipletests
from scipy.special import expit
import numdifftools as nd

logger = logging.getLogger(__name__)
pd.set_option("display.precision", 15)
np.set_printoptions(precision=15, suppress=False)



from scipy.optimize import minimize

# Code provided by Eric Talevich and edited by Vladimir from Stack Overflow
def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def safe_minimize(*args, **kwargs):
    result = minimize(*args, **kwargs)
    if not result.success and 'ABNORMAL_TERMINATION' in result.message:
        logger.warning("[WARN] Retry on ABNORMAL_TERMINATION_IN_LNSRCH")
        kwargs["options"] = kwargs.get("options", {})
        kwargs["options"]["maxiter"] = 1000
        return minimize(*args, **kwargs)
    return result


MAX_EXP = 700
from scipy.special import betaln, gammaln, comb

def beta_binomial_ll(params, X, Z, y, n):
    try:
        eps = 1e-12  # Stability constant

        beta = params[:X.shape[1]]
        gamma = params[X.shape[1]:]

        eta = np.clip(X @ beta, -MAX_EXP, MAX_EXP)
        phi_raw = np.clip(Z @ gamma, -MAX_EXP, MAX_EXP)

        p = np.clip(expit(eta), eps, 1 - eps)
        phi = np.clip(expit(phi_raw), 1e-5, 1 - 1e-5)  # Tight bounds help with stable alpha/beta

        alpha = np.clip(p * (1 - phi) / phi, eps, 1e6)
        beta_param = np.clip((1 - p) * (1 - phi) / phi, eps, 1e6)

        binom_mask = phi < 1e-8
        bb_mask = ~binom_mask

        loglik = np.zeros_like(y, dtype=np.float64)

        if np.any(binom_mask):
            p_binom = np.clip(p[binom_mask], eps, 1 - eps)
            loglik[binom_mask] = (
                y[binom_mask] * np.log(p_binom) +
                (n[binom_mask] - y[binom_mask]) * np.log(1 - p_binom) +
                gammaln(n[binom_mask] + 1) -
                gammaln(y[binom_mask] + 1) -
                gammaln(n[binom_mask] - y[binom_mask] + 1)
            )

        if np.any(bb_mask):
            a = alpha[bb_mask]
            b = beta_param[bb_mask]
            yb = y[bb_mask]
            nb = n[bb_mask]

            loglik[bb_mask] = (
                gammaln(nb + 1) -
                gammaln(yb + 1) -
                gammaln(nb - yb + 1) +
                betaln(a + yb, b + nb - yb) -
                betaln(a, b)
            )

        if not np.all(np.isfinite(loglik)):
            logger.warning("[WARN] Non-finite log-likelihoods encountered")

        return -np.sum(loglik)

    except Exception as e:
        logger.error(f"[ERROR] {e}")
        raise



def circ_test(Circ, Linear, group, alpha=0.05, plotsig=True, CircCoordinates=None, circle_description=None):
    logger.info(f"[INFO] circ_test input rows: {Circ.shape[0]}")  
    logger.info(f"[INFO] Circ shape: {Circ.shape}")
    logger.info(f"[INFO] Linear shape: {Linear.shape}")
    logger.info(f"[INFO] Group length: {len(group)}")

    if Circ.shape != Linear.shape:
        raise ValueError("Circ and Linear must have the same shape.")

    if circle_description is None:
        raise ValueError("circle_description must be provided (list of column indices or names).")

    # Convert index-based to names, if needed
    if all(isinstance(i, int) for i in circle_description):
        circle_description = list(Circ.columns[circle_description])

    if len(group) != Circ.shape[1] - len(circle_description):
        raise ValueError("Length of 'group' must match number of sample columns.")

    group = pd.Series(group)
    logger.debug("Group counts:")
    logger.debug(group.value_counts()) 
    unique_groups = sorted(np.unique(group))
    group_lookup = {g: i for i, g in enumerate(unique_groups)}
    group_ids = group.map(group_lookup).values

    X = pd.get_dummies(group, drop_first=False).astype(float)
    Z = np.ones((len(group), 1))  # phi ~ 1

    X_null = np.ones((len(group), 1))  # Intercept-only for Null

    logger.info(f"[INFO] Detected groups: {list(unique_groups)}")

    p_vals = []
    ratio_means_df = pd.DataFrame(index=Circ.index)
    for grp in unique_groups:
        ratio_means_df[f"group_{grp}_ratio_mean"] = np.nan
        
    options = {
        'maxiter': 2000,
        'xatol': 1e-8,     # Parameter convergence tolerance
        'fatol': 1e-8,     # Function convergence tolerance
        'disp': False      # Turn off verbose output
    }


    for row_num, (idx, row) in enumerate(Circ.iterrows(), start=1):
        gene_name = CircCoordinates.loc[idx, "Gene"] if CircCoordinates is not None and idx in CircCoordinates.index else "Unknown"
        logger.info(f"[INFO] Processing row {row_num}: index = {idx}, gene = {gene_name}")

        try:
            circ_vals = pd.to_numeric(row.drop(labels=circle_description), errors='coerce').astype(float).values
            linear_vals = pd.to_numeric(Linear.loc[idx].drop(labels=circle_description), errors='coerce').astype(float).values

            # Clean up NaNs to ensure sum works
            circ_vals = np.nan_to_num(circ_vals, nan=0.0)
            linear_vals = np.nan_to_num(linear_vals, nan=0.0)

            # Now safely sum
            total = circ_vals + linear_vals

            # Fix 0s like R script
            total[total == 0] = 1

            # Compute ratios
            ratios = circ_vals / total

            



            # Debugging output
            logger.debug(f"\n[DEBUG] Row: {idx}")
            logger.debug(f"[DEBUG] circ_vals: {circ_vals}")
            logger.debug(f"[DEBUG] linear_vals: {linear_vals}")
            logger.debug(f"[DEBUG] total: {total}")
            logger.debug(f"[DEBUG] ratios: {ratios}")
            logger.debug(f"[DEBUG] group: {group.tolist()}")

            df = pd.DataFrame({'circ': circ_vals, 'total': total, 'group': group, 'ratio': ratios})

            for grp in unique_groups:
                group_mask = df['group'] == grp
                mean_ratio = df.loc[group_mask, 'ratio'].mean(skipna=True)

                # fallback if the mean becomes NaN (i.e., no valid ratios)
                if pd.isna(mean_ratio):
                    mean_ratio = 0.0
                logger.debug(f"[DEBUG] Group: {grp}, N={group_mask.sum()}, Mean Ratio={mean_ratio}")
                ratio_means_df.loc[idx, f"group_{grp}_ratio_mean"] = mean_ratio


            y = df['circ'].astype(int).values
            n = df['total'].astype(int).values

            init_alt = np.concatenate([np.zeros(X.shape[1]), np.full(Z.shape[1], 0.1)])
            init_null = np.concatenate([np.zeros(X_null.shape[1]), np.full(Z.shape[1], 0.1)])


            
            bounds_alt = [(-10, 10)] * (X.shape[1] + Z.shape[1])
            bounds_null = [(-10, 10)] * (X_null.shape[1] + Z.shape[1])


            # Try optimization with retries
            alt_attempts = 0
            res_alt = None
            while alt_attempts < 10:
                if alt_attempts == 0:
                    init_try = init_alt
                elif alt_attempts == 1:
                    init_try = np.random.normal(0, 0.5, X.shape[1] + Z.shape[1])
                else:
                    init_try = np.random.uniform(-1, 1, X.shape[1] + Z.shape[1])

                res_alt = safe_minimize(
                    beta_binomial_ll,
                    init_try,
                    args=(X.values, Z, y, n),
                    method='Nelder-Mead',
                    options={
                        'maxiter': 2000,
                        'xatol': 1e-8,
                        'fatol': 1e-8,
                        'disp': False
                    }
                )



                if res_alt.success:
                    logger.info(f"[INFO] res alt successful on attempt {alt_attempts + 1}")
                    try:
                        hess_alt = nd.Hessian(lambda p: beta_binomial_ll(p, X.values, Z, y, n))(res_alt.x)
                        logger.debug(f"[DEBUG] Hessian (alt):\n{hess_alt}")
                    except Exception as e:
                        logger.warning(f"[WARN] Failed to compute Hessian for alt model: {e}")
                        hess_alt = None
                    break
                else:
                    logger.warning(f"[WARN] res alt failed on attempt {alt_attempts + 1} with message: {res_alt.message}")
                    alt_attempts += 1

            # Try optimization for null with retries
            null_attempts = 0
            res_null = None
            while null_attempts < 10:
                if null_attempts == 0:
                    init_try_null = init_null
                elif null_attempts == 1:
                    init_try_null = np.random.normal(0, 0.5, X_null.shape[1] + Z.shape[1])
                else:
                    init_try_null = np.random.uniform(-1, 1, X_null.shape[1] + Z.shape[1])

                res_null = safe_minimize(
                    beta_binomial_ll,
                    init_try_null,
                    args=(X_null, Z, y, n),
                    method='Nelder-Mead',
                    options={
                        'maxiter': 2000,
                        'xatol': 1e-8,
                        'fatol': 1e-8,
                        'disp': False
                    }
                )


                if res_null.success:
                    logger.info(f"[INFO] res null successful on attempt {null_attempts + 1}")
                    try:
                        hess_null = nd.Hessian(lambda p: beta_binomial_ll(p, X_null, Z, y, n))(res_null.x)
                        logger.debug(f"[DEBUG] Hessian (null):\n{hess_null}")
                    except Exception as e:
                        logger.warn(f"[WARN] Failed to compute Hessian for null model: {e}")
                        hess_null = None
                    break
                else:
                    logger.warning(f"[WARN] res null failed on attempt {null_attempts + 1} with message: {res_null.message}")
                    null_attempts += 1

            logger.debug(f"res null {res_null}")
            
            try:
                se_alt = np.sqrt(np.diag(np.linalg.inv(hess_alt)))
                logger.info(f"[INFO] Standard errors (alt): {se_alt}")
            except:
                se_alt = None
                logger.warning("[WARN] Failed to compute standard errors for alt model")

            try:
                se_null = np.sqrt(np.diag(np.linalg.inv(hess_null)))
                logger.info(f"[INFO] Standard errors (null): {se_null}")
            except:
                se_null = None
                logger.warning("[WARN] Failed to compute standard errors for null model")


            


            if not res_alt.success or not res_null.success:
                raise RuntimeError("Optimization failed.")

            ll_alt = -res_alt.fun
            logger.debug(f"ll alt{ll_alt}")
            ll_null = -res_null.fun
            logger.debug(f"ll null{ll_null}")
            
            neg_ll_alt = -res_alt.fun
            neg_ll_null = -res_null.fun
            logger.debug(f"ll alt neg {neg_ll_alt}")
            logger.debug(f"ll null neg{neg_ll_null}")
            llr_stat = 2 * (ll_alt - ll_null)

            logger.debug(f"ll stat{llr_stat}")
           
            df_diff = X.shape[1] - X_null.shape[1]
            logger.debug(f"df diff{df_diff}")
            p_val = chi2.sf(llr_stat, df_diff)
            logger.debug(f"pval{p_val}")

        except Exception as e:
            logger.error(f"Error at row {idx}: {e}", exc_info=True)
            p_val = np.nan
            logger.info(f"p-value: {p_val}")
        

        p_vals.append(p_val)

    # Ensure p_vals is a NumPy array
    p_vals = np.array(p_vals, dtype=np.float64)
    logger.debug(f"numpy array for p vals: {p_vals}")

    # Create a mask for valid (non-NaN) p-values
    valid_mask = ~np.isnan(p_vals)

    # Adjust using custom Benjamini-Hochberg function (exact match to R)
    p_adj = p_adjust_bh(p_vals)

    logger.info(f"Adjusted p-values (BH): {p_adj}")

    # Now filter on adjusted values
    sig_mask = (p_adj <= alpha) & (~np.isnan(p_adj))
    logger.info(f"sig_mask: {sig_mask}")

    sig_dat = Circ.loc[sig_mask].copy()
    sig_dat["raw_p"] = p_vals[sig_mask]  # ← add raw p-values directly
    sig_ratios = ratio_means_df.loc[sig_mask]
    sig_p = p_adj[sig_mask]

    logger.info(f"sig p_values: {sig_p}")



    logger.info(f"[INFO] Total candidates tested: {len(p_vals)}")
    logger.info(f"[INFO] Significant candidates found: {sig_dat.shape[0]}")

    if CircCoordinates is not None:
        coords = CircCoordinates.loc[sig_dat.index]
        summary_table = pd.concat([
            coords,
            pd.Series(sig_p, index=sig_dat.index, name="sig_p"),
            sig_dat["raw_p"],
            sig_ratios
        ], axis=1)
    else:
        summary_table = pd.concat([
            Circ.loc[sig_mask, circle_description],
            pd.Series(sig_p, index=sig_dat.index, name="sig_p"),
            sig_dat["raw_p"],
            sig_ratios
        ], axis=1)



    return {
        "summary_table": summary_table,
        "sig.dat": sig_dat,
        "p.val": p_vals.tolist(),
        "p.adj": p_adj.tolist(),
        "sig_p": sig_p.tolist(),
        "ratios": sig_ratios,
        "hessian_alt": hess_alt,
        "hessian_null": hess_null,
        "standard_errors_alt": se_alt if 'se_alt' in locals() else None,
        "standard_errors_null": se_null if 'se_null' in locals() else None
    }




def circ_max_perc(circ_vals, linear_vals, Nreplicates):
    """
    Calculate the maximum proportion of circRNA reads over total (circ + linear) across all groups.
    """
    circ = np.array(circ_vals, dtype=float)
    linear = np.array(linear_vals, dtype=float)

    if len(circ) != len(linear):
        raise ValueError("Number of circRNA and linear samples do not match.")

    Ngroups = len(circ) // Nreplicates
    circ_sums = [sum(circ[i*Nreplicates:(i+1)*Nreplicates]) for i in range(Ngroups)]
    linear_sums = [sum(linear[i*Nreplicates:(i+1)*Nreplicates]) for i in range(Ngroups)]
    total = np.array(circ_sums) + np.array(linear_sums)

    with np.errstate(divide='ignore', invalid='ignore'):
        perc_vector = np.true_divide(circ_sums, total)
        perc = np.nanmax(perc_vector)

    logger.debug(f"[DEBUG] circ_sums: {circ_sums}")
    logger.debug(f"[DEBUG] linear_sums: {linear_sums}")
    logger.debug(f"[DEBUG] percentage_vector: {np.round(perc_vector * 100, 2)}")
    logger.debug(f"[DEBUG] max_percentage: {round(perc * 100, 2) if not np.isnan(perc) else 0}%")

    return perc if not np.isnan(perc) else 0


def circ_filter(circ, linear, Nreplicates=3, filter_sample=4, filter_count=5, percentage=1, circle_description=[0,1,2]):
    """
    Filter circRNA based on count and percentage thresholds.
    """
    if circ.shape != linear.shape:
        raise ValueError("Circ and linear data must have the same dimensions.")

    sample_cols = [col for i, col in enumerate(circ.columns) if i not in circle_description]
    keep_rows = []

    logger.debug("[DEBUG] Starting circ_filter...")
    for idx in circ.index:
        circ_counts = circ.loc[idx, sample_cols].astype(float).values
        linear_counts = linear.loc[idx, sample_cols].astype(float).values

        n_above_threshold = np.sum(circ_counts >= filter_count)
        max_perc = circ_max_perc(circ_counts, linear_counts, Nreplicates)

        logger.debug(f"\n[DEBUG] Row {idx}")
        logger.debug(f"    circ_counts >= {filter_count}: {n_above_threshold} samples")
        logger.debug(f"    max_percentage: {round(max_perc * 100, 2)}%")

        if n_above_threshold < filter_sample or max_perc < (percentage / 100):
            logger.debug(f"    -> Row {idx} removed")
        else:
            logger.debug(f"    -> Row {idx} kept")
            keep_rows.append(idx)


    logger.debug(f"\n[DEBUG] Final: Kept {len(keep_rows)} / {circ.shape[0]} rows")
    

    return circ.loc[keep_rows] if keep_rows else circ.iloc[[], :]


    

import pandas as pd
import pysam
from multiprocessing import Pool
import os
from functools import partial

def read_one_region(bamfile, chrom, start, end, strand, spliced_only=True, lib_stranded="unstranded"):
    """
    Count spliced reads overlapping a region, considering strand if specified.
    """
    samfile = pysam.AlignmentFile(bamfile, "rb")

    count = 0
    for read in samfile.fetch(chrom, start, end):
        # Quality filters
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail:
            continue

        # Spliced read: 'N' in CIGAR
        if spliced_only and not any(c[0] == 3 for c in read.cigartuples):  # 3 = 'N' in CIGAR
            continue

        # Strand filtering
        if lib_stranded == "fs":
            read_strand = '-' if read.is_read1 else '+'
        elif lib_stranded == "ss":
            read_strand = '+' if read.is_read1 else '-'
        else:
            read_strand = '*'

        if strand != '*' and read_strand != strand:
            continue

        count += 1

    samfile.close()
    return count

def host_gene_count(CircCoordinates, bamfiles, ncore=2, stranded='unstranded', spliced_only=True):
    """
    Count host gene reads flanking circRNA junctions using BAM files.
    """
    if isinstance(CircCoordinates, str):
        CircCoordinates = pd.read_csv(CircCoordinates, sep='\t', header=None)

    # Assign strand column
    if stranded == 'fs':
        strand_info = CircCoordinates.iloc[:, 5]
    elif stranded == 'ss':
        strand_info = CircCoordinates.iloc[:, 5].apply(lambda x: '-' if x == '+' else '+')
    else:
        strand_info = ['*'] * len(CircCoordinates)

    results = []

    for bamfile in bamfiles:
        if not os.path.exists(bamfile) or not os.path.exists(bamfile + '.bai'):
            raise FileNotFoundError(f"BAM or BAI file not found: {bamfile}")

        left_counts = []
        right_counts = []

        args = zip(
            CircCoordinates.iloc[:, 0],                      # chrom
            CircCoordinates.iloc[:, 1],                      # start
            CircCoordinates.iloc[:, 2],                      # end
            strand_info
        )

        with Pool(ncore) as p:
            left_args = [(bamfile, chrom, start, start + 1, strand, spliced_only, stranded) for chrom, start, _, strand in args]
            right_args = [(bamfile, chrom, end - 1, end, strand, spliced_only, stranded) for chrom, _, end, strand in args]

            left_counts = p.starmap(read_one_region, left_args)
            right_counts = p.starmap(read_one_region, right_args)

        avg_counts = [int((l + r) / 2) for l, r in zip(left_counts, right_counts)]
        results.append(avg_counts)

    count_df = pd.DataFrame(results).T
    count_df.columns = [os.path.basename(b).split('.')[0] for b in bamfiles]
    final_df = pd.concat([CircCoordinates.reset_index(drop=True), count_df], axis=1)

    return final_df


import pandas as pd
import numpy as np
from scipy import stats

def summary_se(data, measurevar, groupvars=None, na_rm=False, conf_interval=0.95):
    """
    Summary statistics: count, mean, standard deviation, standard error, and confidence interval.

    Parameters:
        data (pd.DataFrame): The data frame.
        measurevar (str): The name of the column to summarize.
        groupvars (list of str): Column names to group by.
        na_rm (bool): Whether to remove NA/null values.
        conf_interval (float): Confidence interval (e.g. 0.95).

    Returns:
        pd.DataFrame: Summary table with N, mean, sd, se, and ci.
    """
    if groupvars is None:
        groupvars = []

    def summarize(group):
        vals = group[measurevar]
        if na_rm:
            vals = vals.dropna()

        N = len(vals)
        mean = np.mean(vals)
        sd = np.std(vals, ddof=1)
        se = sd / np.sqrt(N) if N > 0 else np.nan
        ci_mult = stats.t.ppf(conf_interval/2 + 0.5, df=N-1) if N > 1 else np.nan
        ci = se * ci_mult if N > 1 else np.nan

        return pd.Series({'N': N, measurevar: mean, 'sd': sd, 'se': se, 'ci': ci})

    summary = data.groupby(groupvars, dropna=True).apply(summarize).reset_index()

    return summary
