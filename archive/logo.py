#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")

def parse_pos_list(s):
    if pd.isna(s):
        return []
    s = str(s).strip()
    if not s:
        return []
    return [int(x) for x in s.split(";") if x != ""]

def sanity_check(df, pos_col, aa_col, mask_col):
    missing = [c for c in [pos_col, aa_col, mask_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}\nFound: {list(df.columns)}")

    df = df.copy()
    df["_pos"] = df[pos_col].apply(parse_pos_list)
    df[aa_col] = df[aa_col].astype(str)
    df[mask_col] = df[mask_col].astype(str)

    bad = df[
        (df["_pos"].str.len() != df[aa_col].str.len()) |
        (df["_pos"].str.len() != df[mask_col].str.len())
    ]
    if len(bad) > 0:
        cols = [c for c in ["pdb_id", "chain", pos_col, aa_col, mask_col] if c in bad.columns]
        raise ValueError("Length mismatch rows:\n" + bad[cols].to_string(index=False))

    return df

def align_rows(df, pos_col, aa_col, mask_col, window=None):
    all_pos = [p for row in df["_pos"] for p in row]
    if not all_pos:
        raise ValueError("No UniProt positions parsed.")

    if window is None:
        lo, hi = min(all_pos), max(all_pos)
    else:
        lo, hi = window
    axis_positions = list(range(lo, hi + 1))

    def build_aligned(row):
        pos = row["_pos"]
        aa  = row[aa_col]
        ms  = row[mask_col]
        aa_by_pos = {pos[i]: aa[i] for i in range(len(pos))}
        ms_by_pos = {pos[i]: ms[i] for i in range(len(pos))}
        aligned_aa = "".join(aa_by_pos.get(p, "-") for p in axis_positions)
        aligned_ms = "".join(ms_by_pos.get(p, " ") for p in axis_positions)
        return pd.Series({"aligned_cdr3_aa": aligned_aa, "aligned_saltbridge_mask": aligned_ms})

    out = df.join(df.apply(build_aligned, axis=1))
    return out, axis_positions

def drop_gap_columns(aligned_seqs, aligned_masks, axis_positions, max_gap_frac=0.90):
    # Build char matrix for gap fraction
    S = np.array([list(s) for s in aligned_seqs])  # (n, L)
    gap_frac = (S == "-").mean(axis=0)
    keep = gap_frac <= max_gap_frac
    keep_idx = np.where(keep)[0].tolist()

    axis2 = [axis_positions[i] for i in keep_idx]
    seqs2 = ["".join(s[i] for i in keep_idx) for s in aligned_seqs]
    masks2 = ["".join(m[i] for i in keep_idx) for m in aligned_masks]
    return seqs2, masks2, axis2

def make_prob_matrix(aligned_seqs):
    # index will be 0..L-1 for logomaker compatibility with our heat strip
    L = len(aligned_seqs[0])
    counts = pd.DataFrame(0.0, index=list(range(L)), columns=AA20)
    for s in aligned_seqs:
        for i, aa in enumerate(s):
            if aa in counts.columns:
                counts.loc[i, aa] += 1.0
    denom = counts.sum(axis=1).replace(0, np.nan)
    probs = counts.div(denom, axis=0).fillna(0.0)
    return probs

def saltbridge_frequency(aligned_masks):
    M = np.array([list(m) for m in aligned_masks])  # (n, L)
    return (M == "*").mean(axis=0)                  # (L,)

def plot_logo_with_heatstrip(aligned_seqs, aligned_masks, axis_positions,
                             out_png, out_pdf=None,
                             title=None, tick_step=2):

    probs = make_prob_matrix(aligned_seqs)
    sb = saltbridge_frequency(aligned_masks)

    L = len(axis_positions)
    x = np.arange(L)

    L = len(axis_positions)

    assert all(len(s) == L for s in aligned_seqs), "AA alignment length mismatch"
    assert all(len(m) == L for m in aligned_masks), "Salt-bridge mask length mismatch"

    if title is None:
        title = f"CDR3 logo (IGMT-aligned)   n={len(aligned_seqs)}"

    fig = plt.figure(figsize=(max(8, L * 0.28), 4.8))

    # 👇 explicit vertical spacing
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[4.2, 1.0],  # give heatmap breathing room
        hspace=0.35                # THIS is the key
    )

    # ── Top: logo ─────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    logomaker.Logo(probs, ax=ax1)
    ax1.set_title(title, pad=12)
    ax1.set_ylabel("AA frequency")

    ticks = np.arange(0, L, tick_step)
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(
        [axis_positions[i] for i in ticks],
        rotation=90
    )

    # add padding so labels don’t touch the heatmap
    ax1.tick_params(axis="x", pad=6)

    # ── Bottom: heat strip ─────────────────
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    strip = sb.reshape(1, -1)
    im = ax2.imshow(strip, aspect="auto")

    ax2.set_yticks([0])
    ax2.set_yticklabels(["SB freq"])

    ax2.set_xticks(ticks)
    ax2.set_xticklabels(
        [axis_positions[i] for i in ticks],
        rotation=90
    )

    ax2.set_xlabel("UniProt position", labelpad=10)
    ax2.tick_params(axis="x", pad=4)

    # colorbar (anchored to bottom panel only)
    cbar = fig.colorbar(
        im,
        ax=ax2,
        orientation="vertical",
        fraction=0.04,
        pad=0.02
    )
    cbar.set_label("P(salt bridge)")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    if out_pdf:
        fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--pos_col", default="cdr3_uniprot_positions")
    ap.add_argument("--aa_col", default="cdr3_aa")
    ap.add_argument("--mask_col", default="cdr3_saltbridge_mask")

    ap.add_argument("--window_lo", type=int, default=None)
    ap.add_argument("--window_hi", type=int, default=None)

    ap.add_argument("--filter_expr", default=None,
                    help='Optional pandas query, e.g. \'chain == "D"\'')

    ap.add_argument("--max_gap_frac", type=float, default=0.90,
                    help="Drop columns where gap fraction > this (default 0.90)")
    ap.add_argument("--tick_step", type=int, default=2,
                    help="Label every Nth UniProt position (default 2)")

    ap.add_argument("--out_prefix", default="cdr3_logo_fixed")
    args = ap.parse_args()

    df = read_table(args.infile)
    if args.filter_expr:
        df = df.query(args.filter_expr)

    df = sanity_check(df, args.pos_col, args.aa_col, args.mask_col)

    window = None
    if args.window_lo is not None and args.window_hi is not None:
        window = (args.window_lo, args.window_hi)

    aligned_df, axis_positions = align_rows(df, args.pos_col, args.aa_col, args.mask_col, window=window)
    aligned_seqs = aligned_df["aligned_cdr3_aa"].tolist()
    aligned_masks = aligned_df["aligned_saltbridge_mask"].tolist()

    # Drop mostly-gap columns so the logo isn’t mostly empty space
    aligned_seqs, aligned_masks, axis_positions = drop_gap_columns(
        aligned_seqs, aligned_masks, axis_positions, max_gap_frac=args.max_gap_frac
    )

    out_png = f"{args.out_prefix}.png"
    out_pdf = f"{args.out_prefix}.pdf"
    plot_logo_with_heatstrip(
        aligned_seqs, aligned_masks, axis_positions,
        out_png=out_png, out_pdf=out_pdf,
        tick_step=args.tick_step
    )

    # write aligned table for debugging
    aligned_df.to_csv(f"{args.out_prefix}.aligned.tsv", sep="\t", index=False)

    print("Wrote:", out_png)
    print("Wrote:", out_pdf)
    print("Wrote:", f"{args.out_prefix}.aligned.tsv")

if __name__ == "__main__":
    main()





#python logo_heatstrip_fixed.py your_table.csv \
#  --window_lo 80 --window_hi 110 \
#  --tick_step 5 \
#  --out_prefix cdr3_80_110_fixed_ticks5

#python logo_heatstrip_fixed.py your_table.csv \
#  --window_lo 80 --window_hi 110 \
#  --max_gap_frac 0.70 \
#  --out_prefix cdr3_80_110_fixed_gap70

#python logo.py out/cdr3_extracted_tpmhc.csv \
#  --window_lo 80 --window_hi 110 \
#  --out_prefix cdr3_80_110_fixed_pmhc
