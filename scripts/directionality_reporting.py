"""Descriptive donor support, independent of effect size or significance tests."""

import numpy as np


REPORTING_VERSION = "donor-support-v1"


def summarize_donors(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("Donor contrasts must be a finite one-dimensional vector")
    n = len(values)
    positive, negative, zero = (int(np.sum(values > 0)), int(np.sum(values < 0)), int(np.sum(values == 0)))
    mean = float(np.mean(values)) if n else np.nan
    median = float(np.median(values)) if n else np.nan
    if not n:
        support = "not estimable"
    elif n == 1:
        support = "single donor pair; no across-donor replication"
    elif positive and negative:
        support = "mixed donor directions"
    elif zero == n:
        support = "no observed contrast at this scale"
    else:
        support = f"{positive or negative} of {n} {'positive' if positive else 'negative'}"
    agree = positive if mean > 0 else negative if mean < 0 else np.nan
    return {
        "mean_log2_contrast": mean, "median_log2_contrast": median,
        "n_pairs": n, "n_positive": positive, "n_negative": negative, "n_zero": zero,
        "n_agree_mean": agree, "agreement_with_mean": agree / n if n else np.nan,
        "mean_direction": "not estimable" if not n else "positive" if mean > 0 else "negative" if mean < 0 else "zero",
        "donor_support": support,
    }


def stimulation_summary(values):
    result = summarize_donors(values)
    result["mean_log2fc_stim_vs_ctrl"] = result.pop("mean_log2_contrast")
    result["median_log2fc_stim_vs_ctrl"] = result.pop("median_log2_contrast")
    result["direction"] = {"positive": "up_in_stim", "negative": "down_in_stim",
                           "zero": "zero_mean", "not estimable": "not_tested"}[result["mean_direction"]]
    result["normalization_status"] = "not assessed; run normalization sensitivity"
    return result
