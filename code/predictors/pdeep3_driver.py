# -*- coding: utf-8 -*-
"""pdeep3_driver.py — isoPSV 的 pDeep3 调用桥（v1.0，2026-10-04；同日实弹核验修订）。

用法：
  python pdeep3_driver.py --in <spec_id/peptide/modinfo/charge TSV> \
      --out <spec_id/b1/b2/y1/y2 TSV> [--instrument Lumos] [--ce 30]

输入列：spec_id, peptide, modinfo（pDeep 族 "2,Carbamidomethyl[C]" 格式，
可为空），charge。输出列：spec_id, b1, b2, y1, y2——强度逗号列表：
  b1/b2 ＝ b1..b_{n-1} 正序；y1/y2 ＝ y1..y_{n-1} 正序（R 侧组块时再翻回位
  点序）。由 R 侧 pred_pdeep3.R 组块。

pDeep3 侧走 API（2026-10-04 本机实弹核验，D:\\Bioinformatics\\pDeep3 ＋
pDeep2 自带 py37/TF1.15 环境）：
  from pDeep.cmd.tune_and_predict import get_prediction
  pred = get_prediction(peptides, tune_psm=None, raw=None, instrument=..., ce=...)
  m = pred.GetIntensities(seq, modinfo, charge)   # ndarray [n-1, 8]
  列序（官方文档）：[b+, b++, y+, y++, b-ModLoss+, b-ModLoss++, y-ModLoss+,
  y-ModLoss++]；行序＝酶切位点 1..n-1 ⇒ b 列即 b1..b_{n-1} 正序；y 列行 k
  ＝互补 y_{n-k}，即位点序 [y_{n-1}..y1]——与 pDeep2 predict.txt 的 y+1/y+2
  行同构（PEPTIDEK/Lumos/ce30 双引擎对比实证：pDeep3 y 列与 pDeep2 y 行
  形状逐位对应，故 y 正序＝行序倒置 m[::-1, col]）。
"""
import argparse
import csv
import sys


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--instrument", default="Lumos")
    ap.add_argument("--ce", default="30")
    return ap.parse_args()


def import_pdeep3():
    """按约定布局探测 pDeep3 根并导入：①环境变量 PDEEP3_ROOT；②PYTHONPATH/
    已安装环境；③python 所在 env 的上级兄弟目录（<env>/../../pDeep3，即
    D:\\Bioinformatics\\pDeep2\\py37 → D:\\Bioinformatics\\pDeep3 布局）；
    ④D:/Bioinformatics/pDeep3。"""
    import os
    cands = [os.environ.get("PDEEP3_ROOT", "")]
    py_dir = os.path.dirname(os.path.abspath(sys.executable))
    cands.append(os.path.normpath(os.path.join(py_dir, "..", "..", "pDeep3")))
    cands.append(r"D:\Bioinformatics\pDeep3")
    # v1.0：桌面克隆布局（客户端机 pDeep3 在 Desktop/pDeep3-master）
    cands.append(os.path.join(os.path.expanduser("~"), "Desktop", "pDeep3-master"))
    # v1.0：官方附件包布局 <包根>/predictor_env/{python37,pDeep3}——
    # python 在 predictor_env/python37 时其上一级即 pDeep3 兄弟目录
    cands.append(os.path.join(py_dir, "..", "pDeep3"))
    cands.append(os.path.join(py_dir, "..", "..", "pDeep3"))
    last_err = None
    for root in cands:
        try:
            if root and os.path.isdir(root) and root not in sys.path:
                sys.path.insert(0, root)
            from pDeep.cmd.tune_and_predict import get_prediction
            return get_prediction
        except ImportError as e:
            last_err = e
    raise last_err


def main():
    args = parse_args()
    try:
        get_prediction = import_pdeep3()
    except ImportError as e:
        sys.stderr.write(
            "[pdeep3_driver] 无法导入 pDeep3（%s）——请确认 --python 指向装有 "
            "pDeep 的环境（pFindStudio/pDeep3；或设环境变量 PDEEP3_ROOT 指向"
            "其仓库根；TF1.13-1.15/Py<=3.7 或 TF2）。\n" % e)
        return 2

    rows = []
    with open(args.inp, "r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    if not rows:
        sys.stderr.write("[pdeep3_driver] 输入为空：%s\n" % args.inp)
        return 3

    # v100-6（2026-10-06）：词表剔除——pDeep3 的 get_modification() 词表（与
    # featurize.py mod_feature 同源同键）不含的 ModName[site] 位点修饰，送预测
    # 前剔除（该位点按无修饰预测）；否则 get_mod_x 对词表外键 KeyError、整批
    # 中断（实弹：开放式搜索修饰 Carbamidomethyl[M]/HexNAc[T] 等命不中 1532
    # 条词表，2 条即崩掉 79,140 肽的全批）。
    # pDeep3 vocabulary filter: site-mods missing from get_modification() (the
    # exact key source of featurize.mod_feature) are dropped before prediction
    # (site predicted unmodified); otherwise get_mod_x raises KeyError and the
    # whole batch dies.
    try:
        from pDeep.config.modification import get_modification
        vocab = set(get_modification().keys())
    except Exception as e:
        vocab = None
        sys.stderr.write("[pdeep3_driver] 词表加载失败（%r），跳过剔除\n"
                         % (e,))
    drop_stat = {}
    n_aff = 0
    if vocab is not None:
        for r in rows:
            mi = r["modinfo"] or ""
            if not mi:
                continue
            items = [x for x in mi.strip(";").split(";") if x]
            keep = []
            hit = False
            for it in items:
                name = it.split(",", 1)[1] if "," in it else it
                if name in vocab:
                    keep.append(it)
                else:
                    drop_stat[name] = drop_stat.get(name, 0) + 1
                    hit = True
            if hit:
                n_aff += 1
                r["modinfo"] = ";".join(keep)
        if drop_stat:
            sys.stderr.write(
                "[pdeep3_driver] 词表剔除：%d/%d 肽段共移除 %d 个词表外修饰"
                "（%s），该位点按无修饰预测\n"
                % (n_aff, len(rows), sum(drop_stat.values()),
                   ",".join("%s*%d" % kv for kv in sorted(drop_stat.items()))))

    peptides = [(r["peptide"], r["modinfo"] or "", int(r["charge"])) for r in rows]
    ce = float(args.ce)
    pred = get_prediction(peptides, tune_psm=None, raw=None,
                          instrument=args.instrument, ce=ce)

    out_rows = []
    skip_none = 0
    for r in rows:
        n_frag = len(r["peptide"]) - 1
        if n_frag < 1:
            skip_none += 1
            continue
        try:
            m = pred.GetIntensities(r["peptide"], r["modinfo"] or "",
                                    int(r["charge"]))
        except Exception as e:  # 单肽异常不拖垮整批（下游按缺失处理）
            sys.stderr.write("[pdeep3_driver] spec %s 预测异常：%r\n"
                             % (r["spec_id"], e))
            skip_none += 1
            continue
        if m is None or len(m) < n_frag:
            skip_none += 1
            continue
        fmt = lambda v: ",".join("%.5f" % max(0.0, float(x)) for x in v[:n_frag])
        out_rows.append({
            "spec_id": r["spec_id"],
            "b1": fmt(m[:, 0]),            # b+1 正序 b1..b_{n-1}
            "b2": fmt(m[:, 1]),            # b+2 正序
            "y1": fmt(m[::-1, 2]),         # y+1：位点序倒置 → y1..y_{n-1} 正序
            "y2": fmt(m[::-1, 3]),         # y+2 同上
        })

    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["spec_id", "b1", "b2", "y1", "y2"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    sys.stderr.write("[pdeep3_driver] 预测 %d/%d 条（跳过 %d；instrument=%s "
                     "ce=%g）\n" % (len(out_rows), len(rows), skip_none,
                                    args.instrument, ce))
    return 0


if __name__ == "__main__":
    sys.exit(main())
