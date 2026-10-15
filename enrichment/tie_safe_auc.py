#!/usr/bin/env python3
"""Shared tie-safe AUC/logAUC implementation, factored out of
R1/benchmark/scripts/final_analysis.py so the same, verified algorithm can be
reused for Boltz-2 (and anything else added later) instead of going back
through the manuscript's raw enrich.sh -> DOCKBASE enrich.py path.

Why this exists: enrich.sh ranks compounds with `sort -k22,22n`, which sorts
only on the score column. When two compounds land on the exact same score,
GNU sort falls back to comparing the whole line as text, which happens to put
active-compound IDs before decoy IDs -- so every tied pair gets silently
credited as "correctly separated" even though the method produced identical
scores for both. This is invisible for continuous scores (Vina, Glide, AF3
ipTM/mPAE) where ties are rare, but inflated DiffDock-L's D4/S2 AUC by >13
points, because its confidence is only reported to 2 decimal places.

auc_logauc() below fixes this by grouping compounds by their exact score
value BEFORE ranking, so a tied group of N actives and M decoys moves the
ROC curve diagonally (both true-positive and false-positive rate advance
together) instead of being arbitrarily ordered by ID text.

AUC and logAUC are otherwise the standard DOCKBASE definitions: logAUC is the
Mysinger-Shoichet adjusted logAUC (log-uniform x-axis between false-positive
rate 0.001 and 1.0, minus the random-classifier baseline so random = 0).
"""
import csv
import pandas as pd
import math
from collections import defaultdict


def auc_logauc(scores, lig: list, dec:list):
    """scores: {compound_id: value}, SMALLER = BETTER. Missing compounds are
    ranked worst. Returns (auc, logauc) as percentages, or (None, None) if
    either class is empty."""
    if not lig or not dec:
        return None, None
    WORST = float('inf')
    groups = defaultdict(lambda: [0, 0])
    ligset = set(lig)
    for c in list(lig) + list(dec):
        g = groups[scores.get(c, WORST)]
        g[0 if c in ligset else 1] += 1
    nL, nD = len(lig), len(dec)
    tp = fp = 0
    pts = [(0.0, 0.0)]
    for s in sorted(groups):          # sorted on the score VALUE, not on text
        a, d = groups[s]              # everyone tied at this score, together
        tp += a
        fp += d
        pts.append((100.0 * fp / nD, 100.0 * tp / nL))
    auc = sum((b[0] - a[0]) / 100 * (a[1] + (b[1] - a[1]) / 2) / 100
              for a, b in zip(pts[:-1], pts[1:]))
    LO, HI = 0.001, 1.0
    RAND = (HI - LO) / math.log(10) / math.log10(HI / LO)
    npts, prev = [], None
    for x, y in pts:
        xf, yf = x / 100, y / 100
        if prev and prev[0] < LO <= xf:
            x0, y0 = prev
            npts.append((LO, y0 + (yf - y0) * (LO - x0) / (xf - x0) if xf > x0 else yf))
        if LO <= xf <= HI:
            npts.append((xf, yf))
        prev = (xf, yf)
    area = 0.0
    for a, b in zip(npts[:-1], npts[1:]):
        if b[0] - a[0] < 1e-9:
            continue
        m = (b[1] - a[1]) / (b[0] - a[0])
        area += (b[1] - a[1]) / math.log(10) + (b[1] - m * b[0]) * (
            math.log10(b[0]) - math.log10(a[0]))
    return 100 * auc, 100 * (area / math.log10(HI / LO) - RAND)


def read_csv_scores(path, higher_better=False, id_col=None, score_col=None):
    """Read {compound_id,score} from a CSV (docking scores.csv, Boltz-2
    scores.csv, or a wider split.csv), sign-flipped if higher_better so
    'smaller = better' holds for auc_logauc(). id_col/score_col are column
    names; they default to the 1st and 2nd column."""
    out = {}
    if not path or not __import__('os').path.exists(path):
        return out
    with open(path) as f:
        rd = csv.DictReader(f)
        id_col = id_col or rd.fieldnames[0]
        score_col = score_col or rd.fieldnames[1]
        for r in rd:
            v = r.get(score_col, '')
            if v in ('', 'NA'):
                continue
            try:
                v = float(v)
            except ValueError:
                continue
            out[r[id_col]] = -v if higher_better else v
    return out


def read_paper_extract(path, higher_better, sentinel):
    """Read the original paper's extract_all.sort.uniq.txt (DOCKBASE's
    22-column format; id in column 3, value in column 22). Rows equal to
    `sentinel` mean 'not scored' and are dropped rather than kept as a real
    value. Matches R1/benchmark/scripts/final_analysis.py:read_extract."""
    out = {}
    if not path or not __import__('os').path.exists(path):
        return out
    for line in open(path):
        p = line.split()
        if len(p) < 22:
            continue
        try:
            v = float(p[21])
        except ValueError:
            continue
        if sentinel is not None and abs(v - sentinel) < 1e-9:
            continue
        out[p[2]] = -v if higher_better else v
    return out


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Compute tie-safe AUC/logAUC from a CSV or DOCKBASE extract')
    parser.add_argument('scores', help='CSV or extract_all.sort.uniq.txt')
    parser.add_argument('--scores-col-name', required=True, help='Column name in CSV to use as score')
    parser.add_argument('--id-col-name', default='compound_id',
                        help='Column name in CSV holding the compound ID')
    parser.add_argument('--active-col-name', default='is_active',
                        help='Column name in CSV flagging actives (1) vs decoys (0)')
    parser.add_argument('--lig', help='Actives list (no header, 1 ID per line); '
                        'only needed for non-CSV (DOCKBASE extract) input')
    parser.add_argument('--dec', help='Decoys list; only needed for non-CSV input')
    parser.add_argument('--higher-better', action='store_true',
                        help='Sign-flip scores so larger = better')
    parser.add_argument('--sentinel', type=float, default=None,
                        help='DOCKBASE sentinel value for "not scored"')
    args = parser.parse_args()

    if args.scores.endswith('.csv'):
        df = pd.read_csv(args.scores)
        scores = read_csv_scores(args.scores, args.higher_better,
                                 args.id_col_name, args.scores_col_name)
        lig_list = df.loc[df[args.active_col_name] == 1, args.id_col_name].tolist()
        dec_list = df.loc[df[args.active_col_name] == 0, args.id_col_name].tolist()
    else:
        if not (args.lig and args.dec):
            parser.error('--lig and --dec are required for non-CSV input')
        scores = read_paper_extract(args.scores, args.higher_better, args.sentinel)
        lig_list = pd.read_csv(args.lig, header=None).iloc[:, 0].tolist()
        dec_list = pd.read_csv(args.dec, header=None).iloc[:, 0].tolist()

    n_hit = sum(c in scores for c in lig_list + dec_list)
    if n_hit == 0:
        raise SystemExit('No ligand/decoy IDs matched the score file; '
                         'check --id-col-name / --scores-col-name')
    print(f'scored {n_hit}/{len(lig_list) + len(dec_list)} compounds '
          f'({len(lig_list)} actives, {len(dec_list)} decoys)')

    auc, logauc = auc_logauc(scores, lig_list, dec_list)
    print(f'AUC={auc:.2f} logAUC={logauc:.2f}')
