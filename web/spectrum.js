/* spectrum.js - isoPSV 网页版谱图浏览器
 * spectrum.js - web spectrum viewer for the results report (canvas 2D, no deps).
 *
 * 双面板：上＝候选 1（原始），下＝下拉选择的竞争/校正候选。
 * 理论梯 a/b/c/x/y/z（1+/2+）＋M 系（母离子与 -H2O/-NH3）；
 * 20 ppm 窗最近峰匹配；颜色表与标注规则与 PDF v4 一致：
 *   电荷号不省略（y10+/y10++）；同峰多离子合并单标签（优先序 y>b>a>c>x>z>M，
 *   b/y 为准、b 与 y 同现以 y 为准）；相邻峰（≤50 Da）标注窜开（低丰度正常、
 *   高丰度抬到低丰度标签之上）。
 */
(function () {
  "use strict";
  var ION_COL = { a: "#004000", b: "#008000", c: "#CA1D52", x: "#650E29",
                  y: "#CA1D52", z: "#943AA4", M: "#808080", I: "#FF8C00", O: "#9ACD32" };
  var ION_PRI = { y: 1, b: 2, a: 3, c: 4, x: 5, z: 6, M: 7, I: 8, O: 9 };
  var CO = 27.994915, NH3 = 17.026549, H2O = 18.010565, H = 1.007276;
  var MOD = null;                       // /api/modmass 结果 {aa:{}, mods:{}}
  var CUR = null;                       // /api/spectrum 结果

  function aaMass(ch) {                 // 氨基酸质量（J=Xle=I/L）
    if (ch === "J" || ch === "I" || ch === "L") return (MOD.aa.I + MOD.aa.L) / 2;
    return MOD.aa[ch];
  }
  var AA3TO1 = { Ala:"A",Arg:"R",Asn:"N",Asp:"D",Cys:"C",Gln:"Q",Glu:"E",Gly:"G",
                 His:"H",Ile:"I",Leu:"L",Lys:"K",Met:"M",Phe:"F",Pro:"P",Ser:"S",
                 Thr:"T",Trp:"W",Tyr:"Y",Val:"V",Xle:"J",Sec:"U",Pyl:"O" };
  function parseMods(md) {              // "3,Oxidation[M];0,Gln->pyro-Glu[AnyN-termQ]" -> [{pos, mass, sap}]
    var out = [];
    if (!md) return out;
    md.split(";").forEach(function (t) {
      t = t.trim(); if (!t) return;
      var i = t.indexOf(","); if (i < 0) return;
      var pos = parseInt(t.slice(0, i), 10);
      var body = t.slice(i + 1);
      // SAP 型替换修饰 X->Y（残基名目标）——按目标残基质量替换该位（非增量）
      // SAP-type substitution X->Y (residue-name target): replace the residue
      // mass at that position (not an increment)
      var arrow = body.indexOf("->");
      if (arrow > 0) {
        var to3 = body.slice(arrow + 2).replace(/\[.*$/, "").trim();
        var to1 = AA3TO1[to3] || null;
        if (to1) {
          out.push({ pos: pos, mass: 0, sap: true, to: to1, from: body.slice(0, arrow).trim() });
          return;
        }
      }
      var name = body.replace(/\[.*$/, "");   // 去 [M] 定位后缀
      var mass = MOD.mods.hasOwnProperty(name) ? MOD.mods[name]:
                 (MOD.mods.hasOwnProperty(body) ? MOD.mods[body]: 0);
      out.push({ pos: pos, mass: mass || 0 });
    });
    return out;
  }
  /* 理论离子表：与 PDF v4 同口径（a=b-CO、c=b+NH3、x=y+CO、z=y-NH3、M 系）；
   * SAP 型 X->Y 替换修饰：该位残基质量替换为目标残基质量（不是加增量） */
  /* Theoretical-ion table (same convention as PDF v4; SAP-type X->Y
   * substitution replaces the residue mass at that position with the target
   * residue mass, not an increment). */
  function theoIons(seq, mod) {
    var n = seq.length, res = [], i;
    var pref = [0], mods = parseMods(mod);
    for (i = 0; i < n; i++) {
      var m = aaMass(seq[i]);
      mods.forEach(function (md) { if (md.sap && md.pos === i + 1) m = aaMass(md.to); });
      pref.push(pref[i] + m);
    }
    var total = pref[n];
    function bMass(k) {                  // b_k 中性质量（含位置≤k 的修饰，pos0=N 端计入全部 b）
      var m = 0;
      mods.forEach(function (md) { if (!md.sap && md.pos <= k) m += md.mass; });
      return pref[k] + m;
    }
    function yMass(k) {                  // y_k 中性质量（末 k 个残基＋H2O；含位置>len-k 修饰）
      var m = 0, lo = n - k + 1;
      mods.forEach(function (md) { if (!md.sap && (md.pos >= lo || md.pos >= n)) m += md.mass; });
      return pref[n] - pref[n - k] + H2O + m;
    }
    function push(mz, lab, ion, num, chg) { res.push({ mz: mz, lab: lab, ion: ion, num: num, chg: chg }); }
    for (var k = 1; k <= n; k++) {
      var bm = bMass(k), ym = yMass(k);
      push(bm + H, "b" + k + "+", "b", k, 1);   push((bm + 2 * H) / 2, "b" + k + "++", "b", k, 2);
      push(bm - CO + H, "a" + k + "+", "a", k, 1); push((bm - CO + 2 * H) / 2, "a" + k + "++", "a", k, 2);
      push(bm + NH3 + H, "c" + k + "+", "c", k, 1); push((bm + NH3 + 2 * H) / 2, "c" + k + "++", "c", k, 2);
      push(ym + H, "y" + k + "+", "y", k, 1);   push((ym + 2 * H) / 2, "y" + k + "++", "y", k, 2);
      push(ym + CO + H, "x" + k + "+", "x", k, 1); push((ym + CO + 2 * H) / 2, "x" + k + "++", "x", k, 2);
      push(ym - NH3 + H, "z" + k + "+", "z", k, 1); push((ym - NH3 + 2 * H) / 2, "z" + k + "++", "z", k, 2);
    }
    var m0 = 0; mods.forEach(function (md) { m0 += md.mass; });
    var pm = total + m0 + H2O;                        // 母离子区（M 系）
    push(pm + H, "M+", "M", null, 1);        push((pm + 2 * H) / 2, "M++", "M", null, 2);
    push(pm + H - H2O, "M-H2O+", "M", null, 1); push((pm + 2 * H - H2O) / 2, "M-H2O++", "M", null, 2);
    push(pm + H - NH3, "M-NH3+", "M", null, 1); push((pm + 2 * H - NH3) / 2, "M-NH3++", "M", null, 2);
    return res;
  }
  /* 20 ppm 窗最近最强峰匹配（窗内无峰不匹配——与 PDF v4 修 bug 后一致） */
  function matchPeaks(theo, peaks) {
    peaks.forEach(function (p) { p.mzSort = p.mz; });
    peaks.sort(function (a, b) { return a.mz - b.mz; });
    theo.forEach(function (t) {
      t.peak = null;
      var lo = t.mz * (1 - 20e-6), hi = t.mz * (1 + 20e-6);
      var best = null;
      peaks.forEach(function (p) {
        if (p.mz >= lo && p.mz <= hi && (!best || p.int > best.int)) best = p;
      });
      t.peak = best;
    });
  }
  /* 合并标签＋窜开（v4 规则），返回每峰 {x,y0,lab,col} */
  function annotate(theo, toX, toY) {
    var hit = theo.filter(function (t) { return t.peak; });
    hit.sort(function (a, b) {
      var d = ION_PRI[a.ion] - ION_PRI[b.ion];
      if (d) return d;
      d = (a.num || 0) - (b.num || 0); if (d) return d;
      return a.chg - b.chg;
    });
    var byPeak = {};
    hit.forEach(function (t) {
      (byPeak[t.peak.mz] = byPeak[t.peak.mz] || []).push(t);
    });
    var items = Object.keys(byPeak).map(Number).sort(function (a, b) { return a - b; })
      .map(function (mz) {
        var g = byPeak[mz];
        return { x: toX(mz), y: toY(g[0].peak.int), lab: g.map(function (t) { return t.lab; }).join(""),
                 col: ION_COL[g[0].ion], int: g[0].peak.int };
      });
    var CH_H = 9;                                     // 每字符纵向像素（11px 字体）
    for (var i = 1; i < items.length; i++) {          // 相邻 ≤50 Da 窜开（像素域等价换算）
      if (items[i].x - items[i - 1].x <= 50) {
        var j = i, k = i - 1;
        if (items[j].int < items[k].int) { var tt = j; j = k; k = tt; }  // 高丰度抬高
        var floorY = items[k].y0 !== undefined ? items[k].y0: items[k].y;
        floorY = floorY - (items[k].lab.length * CH_H + 3);
        if (items[j].y > floorY) items[j].y0 = floorY;
      }
    }
    return items;
  }
  function drawPanel(ctx, W, H, peaks, seq, mod, head) {
    var xmax = peaks.length ? Math.max.apply(null, peaks.map(function (p) { return p.mz; })): 1;
    xmax *= 1.03;
    var imax = Math.max.apply(null, peaks.map(function (p) { return p.int; }));
    var L = 52, R = 8, T = 30, B = 26;
    var toX = function (mz) { return L + (mz / xmax) * (W - L - R); };
    var toY = function (v) { return H - B - (v / imax) * (H - T - B); };
    ctx.clearRect(0, 0, W, H);
    ctx.strokeStyle = "#bbb"; ctx.lineWidth = 1;
    ctx.strokeRect(L, T, W - L - R, H - T - B);
    ctx.fillStyle = "#666"; ctx.font = "10px sans-serif"; ctx.textAlign = "right";
    for (var g = 0; g <= 1.0001; g += 0.25) {
      var yy = toY(g * imax <= imax ? g * imax: imax);
      ctx.fillText((100 * g).toFixed(0) + "%", L - 4, toY(g * imax) + 3);
    }
    ctx.textAlign = "center";
    for (var xm = 0; xm <= xmax; xm += Math.ceil(xmax / 10 / 100) * 100) {
      ctx.fillText(String(xm), toX(xm), H - 8);
    }
    peaks.forEach(function (p) {                      // 灰色基线峰
      ctx.strokeStyle = "#9aa3ab"; ctx.lineWidth = 0.7;
      ctx.beginPath(); ctx.moveTo(toX(p.mz), H - B); ctx.lineTo(toX(p.mz), toY(p.int)); ctx.stroke();
    });
    var theo = theoIons(seq, mod);
    matchPeaks(theo, peaks.map(function (p) { return { mz: p.mz, int: p.int }; }));
    var items = annotate(theo, toX, toY);
    ctx.font = "11px sans-serif";
    items.forEach(function (it) {                     // 命中峰着色＋标签（窜开）
      ctx.strokeStyle = it.col; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.moveTo(it.x, H - B); ctx.lineTo(it.x, it.y); ctx.stroke();
      ctx.fillStyle = it.col;
      var yStart = it.y0 !== undefined ? it.y0: it.y - 4;
      ctx.save(); ctx.translate(it.x, yStart); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "left"; ctx.fillText(it.lab, 0, -3); ctx.restore();
    });
    ctx.fillStyle = "#111"; ctx.font = "12px sans-serif"; ctx.textAlign = "left";
    ctx.fillText(head, L, T - 10);
  }
  function render() {
    if (!CUR || !MOD) return;
    var c1 = CUR.candidates[0] || ["", ""];
    var sel = document.getElementById("ispsv_cand_sel");
    var ci = sel ? parseInt(sel.value, 10): 1;
    if (!(ci >= 0 && ci < CUR.candidates.length)) ci = 0;
    var c2 = CUR.candidates[ci] || c1;
    var cv = document.getElementById("ispsv_canvas");
    var W = cv.width, H = cv.height;
    var ctx = cv.getContext("2d");
    ctx.clearRect(0, 0, W, H);
    var head2 = "bottom cand" + (ci + 1) + ": | " + (CUR.protein2 || "NA") + " | " + c2[0] +
                (c2[1] ? " (" + c2[1] + ")": "") + " |";
    var head1 = "top cand1: | " + (CUR.protein1 || "NA") + " | " + c1[0] +
                (c1[1] ? " (" + c1[1] + ")": "") + " |";
    drawPanel(ctx, W, H / 2, CUR.peaks, c1[0], c1[1], head1);
    ctx.save(); ctx.translate(0, H / 2);
    drawPanel(ctx, W, H / 2, CUR.peaks, c2[0], c2[1], head2);
    ctx.restore();
  }
  function loadSpectrum(fn) {
    fetch("/api/spectrum?file=" + encodeURIComponent(fn)).then(function (r) { return r.json(); })
      .then(function (j) {
        if (!j.ok) { document.getElementById("ispsv_msg").textContent = j.error || "failed"; return; }
        CUR = j;
        document.getElementById("ispsv_msg").textContent =
          j.title + " | " + j.conclusion + " | peaks " + j.peaks.length;
        var sel = document.getElementById("ispsv_cand_sel");
        sel.innerHTML = "";
        j.candidates.forEach(function (c, i) {
          var o = document.createElement("option");
          o.value = String(i); o.textContent = "cand" + (i + 1) + ": " + c[0] + (c[1] ? " (" + c[1] + ")": "");
          sel.appendChild(o);
        });
        if (j.winner_index != null && j.winner_index >= 0 && j.winner_index < j.candidates.length)
          sel.value = String(j.winner_index);      /* 胜选候选 */
        else if (j.candidates.length > 1)
          sel.value = "1";                          /* 无胜选标记时默认候选 2 */
        render();
      });
  }
  /* 谱图 PDF 卡片：扫描清单＋点击窗口内嵌（iframe 自带打印工具条） */
  function loadPdfList(dir) {
    var u = "/api/pdfs" + (dir ? "?dir=" + encodeURIComponent(dir): "");
    fetch(u).then(function (r) { return r.json(); }).then(function (j) {
      var box = document.getElementById("ispsv_pdf_list");
      if (!box) return;
      box.innerHTML = "";
      if (!j.ok || !j.pdfs.length) {
        box.innerHTML = "<p class='hint'>未找到 PDF no PDFs found</p>";
        return;
      }
      j.pdfs.forEach(function (p) {
        var a = document.createElement("a");
        a.href = "#"; a.textContent = p.name + "（" + p.size_mb + " MB）";
        a.style.marginRight = "6px";
        a.onclick = function (e) {
          e.preventDefault();
          document.getElementById("ispsv_pdf_frame").src =
            "/pdf?path=" + encodeURIComponent(p.path);
        };
        box.appendChild(a);
        var nw = document.createElement("a");
        nw.href = "/pdf?path=" + encodeURIComponent(p.path);
        nw.target = "_blank"; nw.textContent = "⧉新窗口";
        nw.style.marginRight = "14px";
        box.appendChild(nw);
      });
      var first = j.pdfs[0];                      /* 首个 PDF 直接嵌入窗口 */
      document.getElementById("ispsv_pdf_frame").src =
        "/pdf?path=" + encodeURIComponent(first.path);
    });
  }
  window.ispsvInit = function () {
    /* 谱图 PDF 卡片接线（与谱图查看器同页共存） */
    var pb = document.getElementById("ispsv_pdf_scan");
    var pi = document.getElementById("ispsv_pdf_dir");
    if (pb) {
      pb.onclick = function () { loadPdfList(pi ? pi.value.trim(): ""); };
      loadPdfList("");
    }
    fetch("/api/modmass").then(function (r) { return r.json(); }).then(function (m) {
      MOD = m;
      var input = document.getElementById("ispsv_file");
      var btn = document.getElementById("ispsv_load");
      var sel = document.getElementById("ispsv_cand_sel");
      btn.onclick = function () { if (input.value.trim()) loadSpectrum(input.value.trim()); };
      input.onkeydown = function (e) { if (e.key === "Enter") btn.onclick(); };
      sel.onchange = render;
      fetch("/api/spectra_list").then(function (r) { return r.json(); })
        .then(function (j) {
          var dl = document.getElementById("ispsv_datalist");
          if (dl && j.ok) j.files.slice(0, 500).forEach(function (f) {
            var o = document.createElement("option"); o.value = f; dl.appendChild(o);
          });
          if (j.ok && j.files.length) loadSpectrum(j.files[0]);
        });
    });
  };
  if (document.readyState !== "loading") window.ispsvInit();
  else document.addEventListener("DOMContentLoaded", window.ispsvInit);
})();
