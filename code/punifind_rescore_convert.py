# -*- coding: utf-8 -*-
# punifind_rescore_convert.py
# 把 pUniFind rescore 的打分 pkl（pscore_res_rankloss_*.pkl 的 d[2] 谱级结果）转换为
# pFind 格式的 fdr0.01 .spectra（File_Name / Sequence / Modification），供 isoPSV
# 流程作为策略 4 输入。口径：pscore_is_target 且 pscore_fdr <= 0.01 的谱级最佳候选
# 各一行；修饰用 pFind 格式（modification_pfind，如 "2,Oxidation[M];"），并归一
# pUniFind 多修饰位点前的多余 "_"（如 "0,Acetyl[ProteinN-term];_3,Oxidation[M];"
# -> "0,Acetyl[ProteinN-term];3,Oxidation[M];"）。
# Convert the pUniFind rescore pkl (spectrum-level results in d[2] of
# pscore_res_rankloss_*.pkl) into a pFind-style fdr0.01 .spectra
# (File_Name / Sequence / Modification) to feed the isoPSV pipeline as strategy 4.
# Rule: one row per spectrum, the best candidate with pscore_is_target and
# pscore_fdr <= 0.01; the modification uses the pFind format
# (modification_pfind, e.g. "2,Oxidation[M];"), normalizing the stray "_"
# pUniFind puts before later modification sites of multi-modified peptides
# (e.g. "0,Acetyl[ProteinN-term];_3,Oxidation[M];"
#   -> "0,Acetyl[ProteinN-term];3,Oxidation[M];").
import argparse
import pickle
import re

def norm_mod(mp):
    # pUniFind 对多修饰肽段的 modification_pfind 会在后续修饰位点前多带 "_"
    # （如 "0,Acetyl[ProteinN-term];_3,Oxidation[M];"），pFind/pDeep2 格式应为
    # "3,Oxidation[M];"——把 ";_" 归一为 ";"、去掉开头 "_"（2026-09-06）
    # for multi-modified peptides pUniFind's modification_pfind prefixes the
    # later sites with "_" (e.g. "0,Acetyl[ProteinN-term];_3,Oxidation[M];");
    # the pFind/pDeep2 format is "3,Oxidation[M];" - normalize ";_" to ";"
    # and drop any leading "_" (2026-09-06)
    mp = mp.replace(";_", ";")
    mp = re.sub(r"^_(?=\d)", "", mp)
    return mp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pkl", help="pscore_res_rankloss_*.pkl path")
    ap.add_argument("out", help="output .spectra path")
    ap.add_argument("--fdr", type=float, default=0.01)
    args = ap.parse_args()

    with open(args.pkl, "rb") as f:
        d = pickle.load(f)
    # d = (raw_pep_stats, {}, raw_spectra) —— 谱级结果在 d[2]
    rows = []
    for raw, lst in d[2].items():
        for e in lst:
            fdr = e.get("pscore_fdr")
            is_target = e.get("pscore_is_target")
            if not is_target:
                continue
            if not isinstance(fdr, (int, float)):
                continue
            if fdr > args.fdr:
                continue
            rows.append((e["spectrum_name"],
                         e.get("peptide_pfind", ""),
                         norm_mod(e.get("modification_pfind", ""))))
    with open(args.out, "w", encoding="utf-8") as f:
        f.write("File_Name\tSequence\tModification\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print("wrote %d rows -> %s" % (len(rows), args.out))

if __name__ == "__main__":
    main()
