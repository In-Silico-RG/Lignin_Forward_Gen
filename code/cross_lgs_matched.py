#!/usr/bin/env python
"""
cross_lgs_matched.py

Cross-algorithm validation, the LGS arm: run the generator that made the original
claim, on the matched SW/HW/GR parameterization.

WHAT LGS ACTUALLY DOES (read from the source, not from our notes)
-----------------------------------------------------------------
docs/03 recorded LGS as "whole-structure rejection sampling: build a complete
candidate, reject the whole thing if its bond profile misses target, redraw from
scratch". **That is wrong**, and the error mattered enough to send an earlier
reconstruction down a blind alley. From
`Sourcecode/java/src/main/java/com/mol/construct/GeneratePolymer.java`
(`getPossiblePolymer` -> `getList` -> `permute`), LGS:

  1. `getList()` builds a FIXED MULTISET of bonds: each bond type appears
     `calPer(pct, dp) = round(dp * pct/100)` times. The target bond composition is
     satisfied BY CONSTRUCTION -- it is never sampled and never rejected.
  2. the same is done for the monomer multiset (the S/G composition), and
  3. `permute()` enumerates ORDERINGS of those multisets (capped by noOfPer);
     the cross-product of monomer x bond orderings gives candidate sequences,
     which are then built and kept only if chemically constructible.

So LGS is a **sequence-controlled enumerator** (its own README's words), not a
rejection sampler. It does not search until it meets the target; it starts at the
target and enumerates the valid arrangements of it. What it discards are
chemically-invalid ARRANGEMENTS (e.g. a 5-5 on a syringyl unit, which has no free
C5) -- never bond profiles.

The consequence for the dead manuscript's claim is direct: LGS's cost measures
the size of the valid-arrangement space of *the target you typed in*, and the
count of structures it emits is a combinatorial property of that target crossed
with the S-blocks-C5 rule. It is not a difficulty of "hitting" anything.

This also means LGS is even more strongly BACKWARD than we described: its target
is imposed deterministically, not merely searched toward. That strengthens the
paper's taxonomy while correcting its mechanism.

METRICS
-------
  structures  -- valid structures emitted (LGS's own count; the enumeration size)
  lgs_time_s  -- the runtime LGS reports itself
  wall_s      -- our wall-clock, for comparison with cross_wallclock_stability

Determinism: LGS enumerates, so `structures` should be identical across repeats
of one config. We repeat anyway to check -- CombinationAlg.getRandomNumberInRange
exists in the source and we want to know whether it reaches the counts.

Requires: java (tested on OpenJDK 21) and the LGS jar. Path is auto-detected;
override with --jar.

Usage:
    python code/cross_lgs_matched.py            # matched SW/HW/GR
    python code/cross_lgs_matched.py --repeats 3
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_common import CONFIGS, DP, results_dir, write_csv, describe

DEFAULT_JAR = os.path.expanduser(
    "~/SOFT/lignin-structure-generator-main/Executable_Jar_and_Config_0225/"
    "lgs-generator-1.0.jar")
LOG4J = os.path.expanduser(
    "~/SOFT/lignin-structure-generator-main/Executable_Jar_and_Config/resources/log4j.xml")

# LGS bond vocabulary: BO4, BB, B5, _4O5, _55, DBDO -- note there is NO beta-1.
# We take each config's beta-O-4 target from the old manuscript and distribute the
# remainder using cross_common's chemistry-appropriate residual for that S/G, with
# the beta-1 share dropped (LGS cannot express it) and the rest renormalized, so
# beta-O-4 still lands exactly on the manuscript target. DBDO is left at 0: it is
# absent from every other generator here, so giving it weight would make LGS run
# different chemistry. Integers, summing to 100, as LGS expects.
LGS_BONDCONFIG = {
    #        BO4  BB  B5  _4O5  _55  DBDO
    "SW": dict(BO4=54, BB=4, B5=17, _4O5=6, _55=19, DBDO=0),
    "HW": dict(BO4=69, BB=10, B5=12, _4O5=3, _55=6, DBDO=0),
    "GR": dict(BO4=77, BB=5, B5=9, _4O5=3, _55=6, DBDO=0),
}

CONFIG_TEMPLATE = """---
dp: {dp}
min_dp: 3
max_dp: 5
bondconfig:
  BO4: {BO4}
  BB: {BB}
  B5: {B5}
  _4O5: {_4O5}
  _55: {_55}
  DBDO: {DBDO}
s_g_ratio: {sg}
monoconfig:
  G: 0
  S: 0
  H: 0
png: false
matrices: false
sdf: false
json: false
"""

RE_STRUCTS = re.compile(
    r"Number of polymers generated=(\d+),\s*No of linear structures\s*:(\d+),"
    r"\s*No\. of branched structures:(\d+)")
RE_TIME = re.compile(r"Time on execution=(\d+):(\d+):(\d+)")


def run_lgs(jar, label, sg_ratio, workdir):
    """One LGS run. Returns (structures, linear, branched, lgs_time_s, wall_s)."""
    res_dir = os.path.join(workdir, "resources")
    if not os.path.isdir(res_dir):
        os.makedirs(res_dir)
    if os.path.isfile(LOG4J):
        shutil.copy(LOG4J, res_dir)

    cfg = dict(LGS_BONDCONFIG[label])
    cfg["dp"] = DP
    cfg["sg"] = sg_ratio
    with open(os.path.join(res_dir, "project_config.yaml"), "w") as fh:
        fh.write(CONFIG_TEMPLATE.format(**cfg))

    # LGS resolves resources/ and output/ relative to the working directory
    t0 = time.time()
    proc = subprocess.Popen(["java", "-jar", jar], cwd=workdir,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _ = proc.communicate()
    wall = time.time() - t0
    out = out.decode("utf-8", "replace")

    m = RE_STRUCTS.search(out)
    if not m:
        sys.stderr.write(out[-2000:])
        raise RuntimeError("could not parse LGS structure count for %s" % label)
    structs, linear, branched = int(m.group(1)), int(m.group(2)), int(m.group(3))

    t = RE_TIME.search(out)
    lgs_t = (int(t.group(1)) * 3600 + int(t.group(2)) * 60 + int(t.group(3))) if t else float("nan")
    return structs, linear, branched, lgs_t, wall


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jar", default=DEFAULT_JAR)
    ap.add_argument("--repeats", type=int, default=3,
                    help="repeats per config (LGS enumerates, so the structure "
                         "count is expected to be identical across repeats)")
    args = ap.parse_args()

    if not os.path.isfile(args.jar):
        sys.exit("LGS jar not found: %s (pass --jar)" % args.jar)

    out_dir = results_dir()
    log_path = os.path.join(out_dir, "cross_lgs_matched.log")
    lines = []

    def log(s):
        lines.append(s)
        print(s)
        sys.stdout.flush()

    log("Cross-algorithm validation -- LGS, matched params")
    log("jar: %s" % args.jar)
    log("LGS is a sequence-controlled ENUMERATOR, not a rejection sampler")
    log("(see docstring; read from GeneratePolymer.java). docs/03's description")
    log("of it as whole-structure rejection sampling is wrong.")
    log("params: DP=%d  repeats=%d" % (DP, args.repeats))
    log("")
    log("%-3s %6s %8s | %11s %9s %10s | %10s %10s"
        % ("cfg", "S/G", "BO4 tgt", "structures", "linear", "branched",
           "lgs_t(s)", "wall_s"))

    rows = []
    workdir = tempfile.mkdtemp(prefix="lgs_matched_")
    try:
        for label, sg, target in CONFIGS:
            for rep in range(args.repeats):
                wd = os.path.join(workdir, "%s_%d" % (label, rep))
                os.makedirs(wd)
                structs, lin, br, lgs_t, wall = run_lgs(args.jar, label, sg, wd)
                rows.append({
                    "config": label, "sg_ratio": sg, "target_bo4": target,
                    "rep": rep, "structures": structs, "linear": lin,
                    "branched": br, "lgs_time_s": lgs_t, "wall_s": wall,
                })
                log("%-3s %6.2f %8d | %11d %9d %10d | %10.1f %10.2f"
                    % (label, sg, LGS_BONDCONFIG[label]["BO4"], structs, lin, br,
                       lgs_t, wall))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    write_csv(os.path.join(out_dir, "cross_lgs_matched.csv"),
              ["config", "sg_ratio", "target_bo4", "rep", "structures",
               "linear", "branched", "lgs_time_s", "wall_s"], rows)

    log("")
    log("=" * 72)
    log("SUMMARY")
    for label, sg, target in CONFIGS:
        s = [r["structures"] for r in rows if r["config"] == label]
        w = [r["wall_s"] for r in rows if r["config"] == label]
        det = "DETERMINISTIC" if len(set(s)) == 1 else "VARIES: %s" % sorted(set(s))
        log("  %-3s structures=%s (%s)   wall_s mean=%.2f sd=%.2f"
            % (label, s[0] if len(set(s)) == 1 else s, det, np.mean(w), np.std(w)))

    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote %s" % os.path.join(out_dir, "cross_lgs_matched.csv"))
    print("Wrote %s" % log_path)


if __name__ == "__main__":
    main()
