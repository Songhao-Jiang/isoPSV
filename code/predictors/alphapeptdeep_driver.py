# -*- coding: utf-8 -*-
"""alphapeptdeep_driver.py — isoPSV 的 AlphaPeptDeep(peptdeep) 调用桥
。

用法：
  python alphapeptdeep_driver.py --in <spec_id/peptide/modinfo/charge TSV> \
      --out <spec_id/b1/b2/y1/y2 TSV> [--instrument Lumos] [--ce 30]

输入列：spec_id, peptide, modinfo（pFind 族 "9,Carbamidomethyl[C]" 格式，可为
空），charge。输出列：spec_id, b1, b2, y1, y2——强度逗号列表（0-1 归一）：
  b1/b2 ＝ b1..b_{n-1} 正序（peptdeep b_z1/b_z2 行本就正序）；
  y1/y2 ＝ y1..y_{n-1} 正序（peptdeep y_z1/y_z2 行为位点序 y_{n-1}..y1，
          此处翻转，R 侧组块时再翻回位点序——与 pDeep3 桥一致）。

修饰映射：pFind 名→alphabase 名（Name@AA；AnyN-term→Protein_N-term@site 0）；
未知名跳过并 stderr 警告（与 ms2pip 适配器同纪律）。
"""
import argparse
import csv
import sys

MOD_MAP = {
    "Oxidation": ("Oxidation", "M"),
    "Carbamidomethyl": ("Carbamidomethyl", "C"),
    "Acetyl": ("Acetyl", "Protein_N-term"),
    "Phospho": ("Phospho", "STY"),
    "Deamidated": ("Deamidated", "NQ"),
    "Methyl": ("Methyl", "KR"),
    "Dimethyl": ("Dimethyl", "KR"),
    "Trimethyl": ("Trimethyl", "KR"),
    "Gln->pyro-Glu": ("Gln->pyro-Glu", "Q"),
    "Glu->pyro-Glu": ("Glu->pyro-Glu", "E"),
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--instrument", default="Lumos")
    ap.add_argument("--ce", default="30")
    ap.add_argument("--batch", default="1024")
    return ap.parse_args()


def pfind_mods_to_peptdeep(modinfo):
    """'9,Carbamidomethyl[C];0,Acetyl[AnyN-term]' -> (mods, mod_sites)"""
    mods, sites = "", ""
    if not modinfo:
        return mods, sites
    out_m, out_s = [], []
    for tok in modinfo.strip(";").split(";"):
        if not tok:
            continue
        pos, name_site = tok.split(",", 1)
        name = name_site[: name_site.find("[")] if "[" in name_site else name_site
        site = name_site[name_site.find("[") + 1: name_site.rfind("]")] if "[" in name_site else ""
        if name not in MOD_MAP:
            sys.stderr.write("[alphapeptdeep_driver] 未知修饰（跳过）: %s\n" % tok)
            continue
        ab_name, ab_site = MOD_MAP[name]
        if site in ("AnyN-term", "ProteinN-term"):
            out_m.append("%s@Protein_N-term" % ab_name)
            out_s.append("0")
        elif site == "AnyC-term":
            out_m.append("%s@Protein_C-term" % ab_name)
            out_s.append("-1")
        else:
            out_m.append("%s@%s" % (ab_name, site or ab_site))
            out_s.append(pos)
    return ";".join(out_m), ";".join(out_s)


def main():
    args = parse_args()
    import pandas as pd
    from peptdeep.pretrained_models import ModelManager

    rows = []
    with open(args.inp, "r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    if not rows:
        sys.stderr.write("[alphapeptdeep_driver] 输入为空：%s\n" % args.inp)
        return 3

    prec = {
        "spec_id": [], "sequence": [], "mods": [], "mod_sites": [],
        "charge": [], "nce": [], "instrument": [],
    }
    for r in rows:
        m, s = pfind_mods_to_peptdeep(r["modinfo"] or "")
        prec["spec_id"].append(r["spec_id"])
        prec["sequence"].append(r["peptide"])
        prec["mods"].append(m)
        prec["mod_sites"].append(s)
        prec["charge"].append(int(r["charge"]))
        prec["nce"].append(float(args.ce))
        prec["instrument"].append(args.instrument)
    precursor_df = pd.DataFrame(prec)

    mgr = ModelManager(mask_modloss=False, device="cpu")
    mgr.load_installed_models()
    frag = mgr.predict_ms2(precursor_df, batch_size=int(args.batch))

    # predict_ms2 会按肽长/序列重排 precursor_df——一律以 spec_id 列对齐，
    # 绝不按行序对齐
    pep_len = {r["spec_id"]: len(r["peptide"]) for r in rows}
    out_rows, skip = [], 0
    for _, prow in precursor_df.iterrows():
        sid = prow["spec_id"]
        n_frag = pep_len[sid] - 1
        if n_frag < 1:
            skip += 1
            continue
        start = int(prow["frag_start_idx"])
        stop = int(prow["frag_stop_idx"])
        blk = frag.iloc[start:stop]
        if len(blk) < n_frag:
            skip += 1
            continue
        b1 = blk["b_z1"].to_list()[:n_frag]
        b2 = blk["b_z2"].to_list()[:n_frag]
        y1 = blk["y_z1"].to_list()[:n_frag][::-1]   # 位点序 -> 正序
        y2 = blk["y_z2"].to_list()[:n_frag][::-1]
        fmt = lambda v: ",".join("%.5f" % max(0.0, float(x)) for x in v)
        out_rows.append({
            "spec_id": sid,
            "b1": fmt(b1), "b2": fmt(b2),
            "y1": fmt(y1), "y2": fmt(y2),
        })

    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["spec_id", "b1", "b2", "y1", "y2"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    sys.stderr.write("[alphapeptdeep_driver] 预测 %d/%d 条（跳过 %d；instrument=%s ce=%s）\n"
                     % (len(out_rows), len(rows), skip, args.instrument, args.ce))
    return 0


if __name__ == "__main__":
    sys.exit(main())
