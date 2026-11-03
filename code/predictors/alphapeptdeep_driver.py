# -*- coding: utf-8 -*-
"""alphapeptdeep_driver.py — isoPSV 的 AlphaPeptDeep(peptdeep) 调用桥
（工作站自写；契约同 pdeep3_driver.py； 修饰名映射层）。

用法：
  python alphapeptdeep_driver.py --in <spec_id/peptide/modinfo/charge TSV> \
      --out <spec_id/b1/b2/y1/y2 TSV> [--instrument Lumos] [--ce 30]

输入列：spec_id, peptide, modinfo（pFind 族 "9,Carbamidomethyl[C]" 格式，可为
空），charge。输出列：spec_id, b1, b2, y1, y2——强度逗号列表（0-1 归一）：
  b1/b2 ＝ b1..b_{n-1} 正序（peptdeep b_z1/b_z2 行本就正序）；
  y1/y2 ＝ y1..y_{n-1} 正序（peptdeep y_z1/y_z2 行为位点序 y_{n-1}..y1，
          此处翻转，R 侧组块时再翻回位点序——与 pDeep3 桥一致）。

修饰名映射层（外部机器 KeyError 修复；合并步骤输出的统一 pFind
族格式在此翻译，与 ms2pip 适配器同"未知名跳过并计数"纪律）：
  1) 位置令牌正规化——pFind 终端位点常带残基后缀（AnyN-termQ / AnyN-termK
     / ProteinN-term 等）。早期实现只认裸 "AnyN-term"，带后缀令牌漏进残基分支，
     拼出 "Gln->pyro-Glu@AnyN-termQ" 之类名字交给 peptdeep，featurize 直接
     KeyError 整批崩。现按前缀识别终端令牌（残基后缀剥离并利用）。
  2) 名字三层解析——别名表 → 运行时 alphabase 词表直查（MOD_INFO_DICT，
     "Name@Site" 复合键，约 3000 条）→ 均失败按"未知名"跳过该肽并计数；
     决不把词表外的名字交给 peptdeep。
  3) 终端位点优先 "残基^Any_N-term" 形（如 Gln->pyro-Glu@Q^Any_N-term，
     mod_sites＝残基 1 基位置），无残基锚时退 Any_N-term/Any_C-term
     （mod_sites 0/-1），再退 Protein_N-term/Protein_C-term。
  4) [site] 缺失时用 pos 处残基或默认位点表补；predict_ms2 兜底
     try/except——异常时回贴已解析修饰键清单后 exit 4（R 侧按预测器
     失败降级该列，不阻断整轮）。
"""
import argparse
import csv
import re
import sys

# 别名表：pFind/引擎槽名 -> alphabase(Unimod) 名。同名条目不必列。
#  首批覆盖现场点名的名字与常见变体；运行时词表直查承担其余。
MOD_ALIAS = {
    "Carboxymethylation": "Carbamidomethyl",
    "Carboxyamidomethylation": "Carbamidomethyl",
    "Oxidation(M)": "Oxidation",
    "Phosphorylation": "Phospho",
    "Deamidation": "Deamidated",
    "Acetylation": "Acetyl",
}

# [site] 缺失且 pos 处残基解析失败时的默认位点（承接早期实现 MOD_MAP 的默认位）
MOD_DEFAULT_SITE = {
    "Oxidation": "M",
    "Carbamidomethyl": "C",
    "Phospho": "S",
    "Deamidated": "N",
    "Methyl": "K",
    "Dimethyl": "K",
    "Trimethyl": "K",
    "Gln->pyro-Glu": "Q",
    "Glu->pyro-Glu": "E",
    "Carbamyl": "K",
}

# pFind 终端位点令牌：AnyN-term / AnyN-termQ / ProteinN-termK / Any_C-term …
TERMINAL_RE = re.compile(r"^(any|protein)[-_]?(n|c)[-_]?term(.*)$", re.I)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--instrument", default="Lumos")
    ap.add_argument("--ce", default="30")
    ap.add_argument("--batch", default="1024")
    return ap.parse_args()


def load_mod_vocab():
    """运行时取 alphabase 词表键集（"Name@Site" 复合键）；取不到返回 None
    （此时退化为别名表＋默认位点表的保守模式）。"""
    import importlib
    for mod_path in ("alphabase.constants.modification", "alphabase.modification"):
        try:
            mod = importlib.import_module(mod_path)
            d = getattr(mod, "MOD_INFO_DICT", None)
            if d:
                return set(d.keys())
        except Exception:
            pass
    return None


def pfind_mods_to_peptdeep(modinfo, peptide, vocab, stats):
    """pFind 族 '0,Gln->pyro-Glu[AnyN-termQ];9,Carbamidomethyl[C]'
    -> (mods, mod_sites)；任一令牌解析失败返回 None（整肽跳过并计数）。
    stats: dict，统计 direct/alias/unknown 名字与未知名清单。"""
    if not modinfo:
        return "", ""
    out_m, out_s = [], []
    for tok in modinfo.strip(";").split(";"):
        if not tok:
            continue
        pos, name_site = tok.split(",", 1)
        name = name_site[: name_site.find("[")] if "[" in name_site else name_site
        site = (name_site[name_site.find("[") + 1: name_site.rfind("]")]
                if "[" in name_site else "")
        base = MOD_ALIAS.get(name, name)
        if base != name:
            stats["alias_names"] += 1
        # ---- 位点分类：终端（可带残基后缀）or 残基 ----
        m = TERMINAL_RE.match(site.strip()) if site else None
        aa_at_pos = ""
        try:
            aa_at_pos = peptide[int(pos) - 1] if pos.isdigit() and 0 < int(pos) <= len(peptide) else ""
        except Exception:
            aa_at_pos = ""
        resolved = None
        if m:
            is_n = m.group(2).lower() == "n"
            is_protein = m.group(1).lower() == "protein"
            anch = (m.group(3) or "").strip() or aa_at_pos  # 残基后缀或 pos 残基
            if is_n:
                cands = []
                if anch:
                    cands.append("%s@%s^Any_N-term" % (base, anch))
                cands += ["%s@Any_N-term" % base, "%s@Protein_N-term" % base]
                site_val = pos if anch and pos.isdigit() else "0"
            else:
                cands = []
                if anch:
                    cands.append("%s@%s^Any_C-term" % (base, anch))
                cands += ["%s@Any_C-term" % base, "%s@Protein_C-term" % base]
                site_val = "-1"
            for key in cands:
                if vocab is None or key in vocab:
                    resolved = (key, site_val)
                    break
        else:
            aa = site.strip() or aa_at_pos or MOD_DEFAULT_SITE.get(base, "")
            if aa:
                key = "%s@%s" % (base, aa)
                if vocab is None or key in vocab:
                    resolved = (key, pos if pos.isdigit() else str(len(peptide)))
        if resolved is None:
            stats["unknown_tokens"].add(tok)
            return None
        out_m.append(resolved[0])
        out_s.append(resolved[1])
        stats["mapped_tokens"] += 1
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

    vocab = load_mod_vocab()
    stats = {"mapped_tokens": 0, "alias_names": 0, "unknown_tokens": set()}
    prec = {
        "spec_id": [], "sequence": [], "mods": [], "mod_sites": [],
        "charge": [], "nce": [], "instrument": [],
    }
    dropped = 0
    for r in rows:
        ms = pfind_mods_to_peptdeep(r["modinfo"] or "", r["peptide"], vocab, stats)
        if ms is None:
            dropped += 1
            continue
        m, s = ms
        prec["spec_id"].append(r["spec_id"])
        prec["sequence"].append(r["peptide"])
        prec["mods"].append(m or "")
        prec["mod_sites"].append(s or "")
        prec["charge"].append(int(r["charge"]))
        prec["nce"].append(float(args.ce))
        prec["instrument"].append(args.instrument)
    if stats["unknown_tokens"]:
        show = ", ".join(sorted(stats["unknown_tokens"])[:15])
        sys.stderr.write("[alphapeptdeep_driver] 未知名/位点跳过 %d 条肽"
                         "（词表直查+别名表均未命中，共 %d 个令牌）：%s\n"
                         % (dropped, len(stats["unknown_tokens"]), show))
    if not prec["spec_id"]:
        sys.stderr.write("[alphapeptdeep_driver] 全部输入肽的修饰无法解析（0/%d 可预测），"
                         "写空表退出\n" % len(rows))
        with open(args.out, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["spec_id", "b1", "b2", "y1", "y2"])
        return 0
    precursor_df = pd.DataFrame(prec)

    mgr = ModelManager(mask_modloss=False, device="cpu")
    mgr.load_installed_models()
    try:
        frag = mgr.predict_ms2(precursor_df, batch_size=int(args.batch))
    except Exception as exc:  # 兜底：词表键不应再漏，但仍可能撞位点组合
        keys = sorted(set(x for ms in precursor_df["mods"] if ms
                          for x in ms.split(";")))
        sys.stderr.write("[alphapeptdeep_driver] predict_ms2 异常：%r\n"
                         "已解析修饰键（前 15）：%s\n"
                         % (exc, ", ".join(keys[:15])))
        return 4

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
    sys.stderr.write("[alphapeptdeep_driver] 预测 %d/%d 条（修饰未知名跳过 %d；"
                     "碎片不足跳过 %d；词表=%s；instrument=%s ce=%s）\n"
                     % (len(out_rows), len(rows), dropped, skip,
                        "alphabase" if vocab else "内置别名表", args.instrument, args.ce))
    return 0


if __name__ == "__main__":
    sys.exit(main())
