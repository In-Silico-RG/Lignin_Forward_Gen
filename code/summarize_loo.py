#!/usr/bin/env python
"""
summarize_loo.py

Collect leave-one-sample-out fold logs written by calibrate_ref.py / calibrate_kmc_ref.py
(--hold-out) and print, per held-out sample, the held-out chi2 beside the in-sample chi2
of the corresponding full fit, plus the sums and their ratio.

    python code/summarize_loo.py 'results/loo16/calibration_ref_loo-*.log' results/calibration_ref.log
    python code/summarize_loo.py 'results/loo60/calibration_ref_loo_sw-*.log' results/calibration_ref_v60_sw.log
"""
import glob
import re
import sys


def insample_chi2(full_log):
    out = {}
    for l in open(full_log):
        m = re.match(r"(\S+)\s+([\d.]+)\s+([\d.]+)\s+(b5:|bb:|b1:|c55:|c4o5:|dc:)", l)
        if m:
            out[m.group(1)] = float(m.group(3))
    return out


def main(pattern, full_log):
    ins = insample_chi2(full_log)
    files = sorted(glob.glob(pattern))
    tot_ho = tot_in = 0.0
    print("%-36s %5s %9s %9s  %s" % ("held-out sample", "S/G", "held-out", "in-sample", "z per ratio"))
    for f in files:
        lines = open(f).read().split("\n")
        i = [k for k, l in enumerate(lines) if l.startswith("held-out samples")]
        if not i:
            continue
        parts = lines[i[0] + 1].split()
        sid, sg, c = parts[0], parts[1], float(parts[2])
        tot_ho += c
        tot_in += ins.get(sid, 0.0)
        print("%-36s %5s %9.2f %9.2f  %s" % (sid, sg, c, ins.get(sid, float("nan")), " ".join(parts[3:])))
    print("%d folds; SUM held-out %.1f vs in-sample %.1f (ratio %.2f)"
          % (len(files), tot_ho, tot_in, tot_ho / tot_in if tot_in else float("nan")))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
