# -*- coding: utf-8 -*-
# MGF 重编号：pParse 标题 raw.scan.scan.charge.fraction.dta 的同 scan 多价态/拆分
# 会导致 MSFragger "duplicated scan numbers" 报错——改为唯一顺序 scan（1..N），
# 标题改写为 <stem>.<n>.<n>.<charge>.0.dta，并输出 n→原始标题 映射表
import sys

src, dst, map_path = sys.argv[1], sys.argv[2], sys.argv[3]
stem = dst.rsplit("/", 1)[-1].rsplit("\\", 1)[-1].replace(".mgf", "")

blocks = []
cur = {"charge": "2"}
for line in open(src, encoding="utf-8", errors="replace"):
    ls = line.strip()
    if ls == "BEGIN IONS":
        cur = {"charge": "2"}
    elif ls.startswith("TITLE="):
        cur["title"] = ls[6:]
    elif ls.startswith("CHARGE="):
        cur["charge"] = ls[7:].rstrip("+").split()[0]
    if "lines" not in cur:
        cur["lines"] = []
    cur["lines"].append(line)
    if ls == "END IONS":
        blocks.append(cur)

with open(dst, "w", encoding="utf-8", newline="\n") as out, \
     open(map_path, "w", encoding="utf-8", newline="\n") as m:
    m.write("new_scan\toriginal_title\n")
    for n, b in enumerate(blocks, 1):
        for line in b["lines"]:
            ls = line.strip()
            if ls.startswith("TITLE="):
                out.write("TITLE=%s.%d.%d.%s.0.dta\n" % (stem, n, n, b["charge"]))
            else:
                out.write(line)
        m.write("%d\t%s\n" % (n, b.get("title", "")))
print("rewrote", len(blocks), "spectra; mapping ->", map_path)
