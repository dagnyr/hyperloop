import pandas as pd

def parse_pos_list(s):
    # expects "91;92;93;..."
    return [int(x) for x in str(s).split(";") if x != ""]

def align_by_uniprot_positions(df, pos_col="cdr3_uniprot_positions",
                              aa_col="cdr3_aa", mask_col="cdr3_saltbridge_mask",
                              window=None):
    # Parse positions
    df = df.copy()
    df["_pos"] = df[pos_col].apply(parse_pos_list)

    # Sanity checks
    bad = df[(df["_pos"].str.len() != df[aa_col].str.len()) | (df["_pos"].str.len() != df[mask_col].str.len())]
    if len(bad) > 0:
        raise ValueError("Some rows have mismatched lengths between positions / aa / mask:\n"
                         + bad[["pdb_id","chain",pos_col,aa_col,mask_col]].to_string(index=False))

    # Pick window
    all_pos = [p for row in df["_pos"] for p in row]
    if window is None:
        lo, hi = min(all_pos), max(all_pos)
    else:
        lo, hi = window  # inclusive

    axis = list(range(lo, hi + 1))

    def build_aligned(row):
        pos = row["_pos"]
        aa = row[aa_col]
        mask = row[mask_col]
        aa_by_pos = {pos[i]: aa[i] for i in range(len(pos))}
        m_by_pos  = {pos[i]: mask[i] for i in range(len(pos))}
        aligned_aa = "".join(aa_by_pos.get(p, "-") for p in axis)
        aligned_m  = "".join(m_by_pos.get(p, " ") for p in axis)
        return pd.Series({"uniprot_lo": lo, "uniprot_hi": hi,
                          "aligned_cdr3_aa": aligned_aa,
                          "aligned_saltbridge_mask": aligned_m})

    out = df.join(df.apply(build_aligned, axis=1))
    return out, axis

# Example usage:
# df = pd.read_csv("cdr3_saltbridges.csv", sep="\t")  # or comma
# aligned, axis = align_by_uniprot_positions(df, window=(80, 110))
# aligned.to_csv("cdr3_aligned_by_uniprot.tsv", sep="\t", index=False)

import numpy as np

def saltbridge_frequency(aligned_df, axis, mask_col="aligned_saltbridge_mask"):
    # mask is same length as axis; '*' = salt bridge
    masks = aligned_df[mask_col].astype(str).to_list()
    M = np.array([list(m) for m in masks])  # shape: (n_rows, n_positions)
    freq = (M == "*").mean(axis=0)
    return pd.DataFrame({"uniprot_pos": axis, "saltbridge_freq": freq})

# freq_df = saltbridge_frequency(aligned, axis)
# freq_df.sort_values("saltbridge_freq", ascending=False).head(10)

pip install logomaker matplotlib pandas numpy

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

def logo_from_aligned(aligned_seqs, axis_positions, title="CDR3 logo (UniProt-aligned)",
                      drop_gap_cols=False):
    # aligned_seqs: list of same-length strings, using '-' for gaps
    L = len(axis_positions)
    if any(len(s) != L for s in aligned_seqs):
        raise ValueError("All aligned sequences must have same length as axis_positions.")

    # Optionally drop columns that are mostly gaps
    cols = np.arange(L)
    if drop_gap_cols:
        gap_frac = np.mean(np.array([list(s) for s in aligned_seqs]) == "-", axis=0)
        cols = cols[gap_frac < 0.9]
        axis_positions = [axis_positions[i] for i in cols]
        aligned_seqs = ["".join(s[i] for i in cols) for s in aligned_seqs]

    # Count matrix: rows = positions, cols = amino acids
    counts = pd.DataFrame(0, index=axis_positions, columns=AA20, dtype=float)

    for s in aligned_seqs:
        for pos, aa in zip(axis_positions, s):
            if aa in counts.columns:
                counts.loc[pos, aa] += 1
            # ignore gaps '-' and any weird letters

    # Convert counts to probabilities (frequency logo)
    probs = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    # Plot
    plt.figure(figsize=(max(8, len(axis_positions) * 0.25), 3))
    logo = logomaker.Logo(probs, shade_below=0.5, fade_below=0.5)
    plt.title(title)
    plt.xlabel("UniProt position")
    plt.ylabel("Frequency")
    plt.xticks(
        ticks=range(len(axis_positions)),
        labels=axis_positions,
        rotation=90
    )
    plt.tight_layout()
    return logo

# Example usage:
# aligned_seqs = aligned_df["aligned_cdr3_aa"].tolist()
# axis_positions = axis  # list of UniProt positions used in alignment
# logo_from_aligned(aligned_seqs, axis_positions)
# plt.show()

# ------------------------------ logo w/ heatmap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

def make_prob_matrix(aligned_seqs, axis_positions):
    counts = pd.DataFrame(0, index=axis_positions, columns=AA20, dtype=float)
    for s in aligned_seqs:
        for pos, aa in zip(axis_positions, s):
            if aa in counts.columns:
                counts.loc[pos, aa] += 1
    probs = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return probs

def saltbridge_freq(aligned_masks, axis_positions):
    # aligned_masks: list of strings same length as axis_positions; '*' indicates bridge
    M = np.array([list(m) for m in aligned_masks])
    freq = (M == "*").mean(axis=0)  # per column
    return pd.Series(freq, index=axis_positions)

def plot_logo_with_bridge_track(aligned_seqs, aligned_masks, axis_positions,
                                title="CDR3 logo + salt-bridge frequency"):
    probs = make_prob_matrix(aligned_seqs, axis_positions)
    sb = saltbridge_freq(aligned_masks, axis_positions).values

    L = len(axis_positions)

    fig = plt.figure(figsize=(max(8, L * 0.25), 4))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4, 1], hspace=0.05)

    # --- Logo ---
    ax1 = fig.add_subplot(gs[0, 0])
    logomaker.Logo(probs, ax=ax1)
    ax1.set_title(title)
    ax1.set_ylabel("Frequency")
    ax1.set_xticks([])  # x labels on bottom track

    # --- Bridge frequency track (bar) ---
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    x = np.arange(L)
    ax2.bar(x, sb)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("SB")
    ax2.set_xlabel("UniProt position")
    ax2.set_xticks(x)
    ax2.set_xticklabels(axis_positions, rotation=90)

    plt.tight_layout()
    return fig

# Usage:
# aligned_seqs  = aligned_df["aligned_cdr3_aa"].tolist()
# aligned_masks = aligned_df["aligned_saltbridge_mask"].tolist()
# axis_positions = axis  # list of UniProt positions
# plot_logo_with_bridge_track(aligned_seqs, aligned_masks, axis_positions)
# plt.show()
