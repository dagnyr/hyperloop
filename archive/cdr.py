#!/usr/bin/env python3
"""
cdr3_from_uniprot_fasta.py

Reads a FASTA containing TCR chain sequences (e.g. your out/tcr_sequences.fa that
includes UniProt sequences) and extracts CDR3 using ANARCI with a window-rescue
fallback. If a matching mask record follows the sequence record, it also slices
the mask down to CDR3.

Input FASTA format expected (pairs, mask optional):
  >PDBID|TCR_chain=X|uniprot=ACC
  SEQUENCE...
  >PDBID|TCR_chain=X|uniprot=ACC|saltbridge_mask_uniprot
  ....*..*....

Usage:
  python cdr3_from_uniprot_fasta.py \
      --fasta out/tcr_sequences.fa \
      --out out/cdr3_extracted.csv
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Optional, Tuple, List

from Bio import SeqIO
from anarci import anarci


CDR3_START = 105
CDR3_END = 117


# ---------- ANARCI helpers ----------

def _extract_cdr3_from_numbering(numbering, offset0: int):
    """
    numbering: [((imgt_pos, ins), aa), ...] for the *window* sequence
    offset0: 0-index start of window in the full sequence
    Returns (cdr3_seq, cdr3_full_positions_1idx)
    """
    cdr3_chars = []
    cdr3_pos_1idx = []
    seq_i = -1  # 0-index within the window (increments on non-gap)

    for (imgt_pos, ins), aa in numbering:
        if aa != "-":
            seq_i += 1
        if aa == "-":
            continue

        if CDR3_START <= imgt_pos <= CDR3_END:
            cdr3_chars.append(aa)
            cdr3_pos_1idx.append(offset0 + seq_i + 1)  # full-seq position (1-indexed)

    cdr3_seq = "".join(cdr3_chars) if cdr3_chars else None
    return cdr3_seq, cdr3_pos_1idx


def _get_numbering_from_anarci(seq_fragment: str):
    """
    Returns numbering list or None.
    """
    numbered, _, _ = anarci([("q", seq_fragment)], scheme="imgt")
    if not numbered or numbered[0] is None or len(numbered[0]) == 0:
        return None

    domain = numbered[0][0]

    numbering = None
    # numbering is the list shaped like [((pos,ins), aa), ...]
    for x in reversed(domain):
        if isinstance(x, list) and x and isinstance(x[0], tuple) and len(x[0]) == 2:
            numbering = x
            break

    return numbering


def anarci_cdr3_positions_with_window_rescue(full_seq: str):
    """
    Try ANARCI on the full sequence. If it fails (or yields no CDR3),
    scan windows and return the best plausible hit.

    Returns:
      cdr3_seq (str|None)
      cdr3_full_positions_1idx (list[int])  # positions in full_seq (1-indexed)
      method (str)  # "full" or "window:start:end" or "fail"
    """
    full_seq = full_seq.strip().replace(" ", "").replace("\n", "").upper()

    # 1) Full sequence attempt
    numbering = _get_numbering_from_anarci(full_seq)
    if numbering is not None:
        cdr3_seq, cdr3_pos = _extract_cdr3_from_numbering(numbering, offset0=0)
        if cdr3_seq and 6 <= len(cdr3_seq) <= 30:
            return cdr3_seq, cdr3_pos, "full"

    # 2) Window rescue
    window_sizes = [90, 100, 110, 120, 130, 140, 150]
    step = 10

    best = None  # (len_cdr3, cdr3_seq, cdr3_pos, method)
    for w in window_sizes:
        if w > len(full_seq):
            continue
        for start in range(0, len(full_seq) - w + 1, step):
            frag = full_seq[start:start + w]
            numbering = _get_numbering_from_anarci(frag)
            if numbering is None:
                continue

            cdr3_seq, cdr3_pos = _extract_cdr3_from_numbering(numbering, offset0=start)
            if not cdr3_seq:
                continue
            if not (6 <= len(cdr3_seq) <= 30):
                continue

            method = f"window:{start+1}:{start+w}"  # 1-indexed display
            cand = (len(cdr3_seq), cdr3_seq, cdr3_pos, method)

            # prefer longer CDR3 if multiple hits appear
            if best is None or cand[0] > best[0]:
                best = cand

    if best is None:
        return None, [], "fail"

    return best[1], best[2], best[3]


# ---------- FASTA parsing (sequence + optional following mask) ----------

def iterate_fasta_pairs(fasta_path: str):
    """
    Yields:
      header, sequence, mask

    Assumes mask immediately follows its sequence record.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    i = 0
    while i < len(records):
        header = records[i].id
        seq = str(records[i].seq)

        mask = None
        if i + 1 < len(records):
            nxt = records[i + 1].id
            if ("saltbridge_mask" in nxt) or (nxt.endswith("saltbridge_mask_uniprot")):
                mask = str(records[i + 1].seq)
                i += 2
            else:
                i += 1
        else:
            i += 1

        yield header, seq, mask


def parse_header_fields(header: str):
    """
    Extract pdb_id, chain_id, uniprot_acc from headers like:
      1ABC|TCR_chain=A|uniprot=P12345
    """
    parts = header.split("|")
    pdb_id = parts[0] if parts else header
    chain_id = None
    uniprot_acc = None

    for p in parts[1:]:
        if p.startswith("TCR_chain="):
            chain_id = p.split("=", 1)[1]
        elif p.startswith("uniprot="):
            uniprot_acc = p.split("=", 1)[1]

    return pdb_id, chain_id, uniprot_acc


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Input FASTA (e.g. out/tcr_sequences.fa)")
    ap.add_argument("--out", required=True, help="Output CSV (e.g. out/cdr3_extracted.csv)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "pdb_id",
            "chain",
            "uniprot_acc",
            "seq_len",
            "cdr3_aa",
            "cdr3_len",
            "cdr3_uniprot_positions",
            "cdr3_saltbridge_mask",
            "method",
            "status",
        ])

        for header, seq, mask in iterate_fasta_pairs(args.fasta):
            pdb_id, chain_id, unp = parse_header_fields(header)

            seq_clean = seq.strip().replace(" ", "").replace("\n", "").upper()
            mask_clean = None
            if mask is not None:
                mask_clean = mask.strip().replace(" ", "").replace("\n", "")
                if len(mask_clean) != len(seq_clean):
                    # still proceed without mask; record mismatch
                    mask_clean = None

            try:
                cdr3_seq, cdr3_pos, method = anarci_cdr3_positions_with_window_rescue(seq_clean)
                if not cdr3_seq:
                    w.writerow([pdb_id, chain_id, unp, len(seq_clean), None, 0, "", None, method, "no_cdr3"])
                    continue

                # slice mask if available
                cdr3_mask = None
                if mask_clean is not None and cdr3_pos:
                    # cdr3_pos are 1-indexed
                    cdr3_mask = "".join(mask_clean[p - 1] for p in cdr3_pos if 1 <= p <= len(mask_clean))

                w.writerow([
                    pdb_id,
                    chain_id,
                    unp,
                    len(seq_clean),
                    cdr3_seq,
                    len(cdr3_seq),
                    ";".join(map(str, cdr3_pos)),
                    cdr3_mask,
                    method,
                    "ok",
                ])
            except Exception as e:
                w.writerow([pdb_id, chain_id, unp, len(seq_clean), None, 0, "", None, "error", f"error:{type(e).__name__}"])

    print("✓ Wrote:", args.out)


if __name__ == "__main__":
    main()
