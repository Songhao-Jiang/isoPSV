# tests/web_smoke.py — 网页安全冒烟（B-5，2026-09-26）
# 起服 → 断言 Host/Origin/白名单/端点行为 → 关服。非零退出即失败。
# 用法：python tests/web_smoke.py <client_root>
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import urllib.error

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(
    os.path.join(os.path.dirname(__file__), ".."))
PORT = 8780
BASE = "http://127.0.0.1:%d" % PORT


def js_raw_newline_in_strings(html):
    # v1.0（2026-09-28）：JS 单/双引号字符串内的真实换行 = 整块 <script> 语法
    # 错误（如 join("\n") 被写成 join(" + 换行 + ")），块内全部函数失效、按钮
    # 点了没反应。逐字符状态机扫描所有 <script> 块，返回问题行号列表；识别正则
    # 字面量（前一个有效字符是操作符/开括号或 return 等关键字时 / 起正则）以免
    # esc() 的 /"/g 之类被误判为字符串起止。
    _kw = {"return", "typeof", "case", "in", "of", "new", "throw", "else",
           "do", "delete", "void", "instanceof", "yield", "await"}
    issues = []
    i = 0
    n = len(html)
    while i < n:
        j = html.find("<script", i)
        if j < 0:
            break
        k = html.find("</script>", j)
        if k < 0:
            break
        block = html[j + html[j:j + 20].find(">") + 1:k]
        state = None
        line = 1
        p = 0
        m = len(block)
        last_sig = ""
        last_word = ""
        while p < m:
            c = block[p]
            if state is None:
                if c == "\n":
                    line += 1
                elif c in " \t\r":
                    pass
                elif c == "/" and block[p + 1:p + 2] == "/":
                    state = "lc"
                elif c == "/" and block[p + 1:p + 2] == "*":
                    state = "bc"
                elif c == "/":
                    if last_sig in "([{=,:;!&|?+-*%<>^~" or last_sig == "" \
                            or last_word in _kw:
                        q = p + 1
                        in_class = False
                        while q < m:
                            ch = block[q]
                            if ch == "\\":
                                q += 2
                                continue
                            if ch == "\n":
                                line += 1
                                q += 1
                                continue
                            if ch == "[":
                                in_class = True
                            elif ch == "]":
                                in_class = False
                            elif ch == "/" and not in_class:
                                break
                            q += 1
                        p = q
                    last_sig = "/"
                    last_word = ""
                elif c == '"':
                    state = "d"
                elif c == "'":
                    state = "s"
                elif c == "`":
                    state = "t"
                elif c.isalnum() or c in "_$":
                    last_word += c
                    last_sig = "w"
                else:
                    last_sig = c
                    last_word = ""
            elif state == "lc":
                if c == "\n":
                    state = None
                    line += 1
            elif state == "bc":
                if c == "\n":
                    line += 1
                elif block[p:p + 2] == "*/":
                    state = None
                    p += 1
            elif state in ("d", "s"):
                if c == "\\":
                    p += 1
                elif c == "\n":
                    issues.append(line)
                    state = None
                    line += 1
                elif (c == '"' and state == "d") or (c == "'" and state == "s"):
                    state = None
            elif state == "t":
                if c == "\\":
                    p += 1
                elif c == "`":
                    state = None
                elif c == "\n":
                    line += 1
            p += 1
        i = k + len("</script>")
    return issues


# N-14/N-20（2026-09-26，GLM）：用临时客户端根做冒烟——config 读写断言需要真实
# config_step1.txt；直接打仓库根会因文件缺失让 /api/fill 断言形同虚设。
CLIENT = tempfile.mkdtemp(prefix="isopsv_smoke_")
shutil.copy(os.path.join(ROOT, "config_step1_template.txt"),
            os.path.join(CLIENT, "config_step1.txt"))

srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "web", "server.py"), CLIENT],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
try:
    time.sleep(3)
    fails = []

    def req(path, method="GET", headers=None, data=None):
        r = urllib.request.Request(BASE + path, method=method, headers=headers or {},
                                   data=data)
        try:
            with urllib.request.urlopen(r, timeout=5) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()

    def ck(label, cond):
        print(("PASS" if cond else "FAIL"), label)
        if not cond:
            fails.append(label)

    # 首页正常（正常 Host）
    st, _ = req("/")
    ck("首页 200", st == 200)
    # v1.0（2026-09-28）：index.html 全部 <script> 块无"引号字符串内真实换行"
    # ——此类换行会让整块脚本语法失败、页面按钮全部失效（label_type 第二类
    # 添加按钮实弹暴露）
    with open(os.path.join(ROOT, "web", "index.html"), encoding="utf-8") as fh:
        _html = fh.read()
    _js_issues = js_raw_newline_in_strings(_html)
    ck("index.html JS 无断行字符串（行号 %s）" % _js_issues, not _js_issues)
    # 恶意 Host → 403（GET 与 POST）
    st, _ = req("/", headers={"Host": "evil.example.com"})
    ck("恶意 Host GET 403", st == 403)
    st, _ = req("/api/stop", method="POST", headers={"Host": "evil.example.com"},
                data=b"{}")
    ck("恶意 Host POST 403", st == 403)
    # 空 Host → 403（G-4）
    st, _ = req("/", headers={"Host": ""})
    ck("空 Host 403", st == 403)
    # 恶意 Origin → 403
    st, _ = req("/pdf?path=x.pdf", headers={"Origin": "http://evil.example.com"})
    ck("恶意 Origin /pdf 403", st == 403)
    # 越界 PDF → 404
    st, _ = req("/pdf?path=C:/Windows/System32/win.ini.pdf")
    ck("越界 PDF 404", st == 404)
    # config_keys 白名单
    st, body = req("/api/config_keys?file=..%2F..%2F..%2FWindows%2Fwin.ini")
    ck("config_keys 越界 ok:false", st == 200 and json.loads(body).get("ok") is False)
    # load_config 白名单（B-1）
    st, body = req("/api/load_config", method="POST",
                   data=json.dumps({"path": "C:/Windows/win.ini"}).encode())
    ck("load_config 越界 ok:false", st == 200 and json.loads(body).get("ok") is False)
    # N-7（2026-09-26，GLM）：合法文件名 + 任意目录残面——旧实现只锁 basename，
    # 会把 /foo/bar/config_step1.txt 里的内容回显出来；新实现 join 到固定 base，
    # 即读到的必须是客户端自己的 config_step1.txt（与直接请求同内容，而非外部文件）
    st1, body1 = req("/api/load_config", method="POST",
                     data=json.dumps({"path": "C:/Windows/config_step1.txt"}).encode())
    st2, body2 = req("/api/load_config", method="POST",
                     data=json.dumps({"path": "config_step1.txt"}).encode())
    fields1 = json.loads(body1).get("fields", {})
    fields2 = json.loads(body2).get("fields", {})
    ck("load_config 合法名+任意目录读的是客户端自身配置（N-7）",
       st1 == 200 and st2 == 200 and bool(fields1) and fields1 == fields2)
    # N-20（2026-09-26，GLM）：_sanitize_inline_field 须剔 tab/NUL——用 /api/fill
    # 写含 \t/\0 的键值，回读 config 不得出现拆列或 NUL 字节
    st, _ = req("/api/fill", method="POST",
                data=json.dumps({"file": "config_step1.txt", "key": "xkey\x00",
                                 "value": "v1\tv2"}).encode())
    ck("/api/fill 含 tab/NUL 200（N-20）", st == 200)
    st, body = req("/api/load_config", method="POST",
                   data=json.dumps({"path": "config_step1.txt"}).encode())
    fields = json.loads(body).get("fields", {})
    ck("config 无 NUL 键、无拆列歧义（N-20）",
       st == 200 and "xkey" in fields and "\x00" not in json.dumps(fields))
    # N-13 机制（2026-09-26）：/api/fill 空值 = 移除该键——前端清理陈旧
    # predictorN_* 键依赖此语义（删行后残留键会让 user_run 执行陈旧预测器）
    st, _ = req("/api/fill", method="POST",
                data=json.dumps({"file": "config_step1.txt", "key": "xkey",
                                 "value": ""}).encode())
    ck("/api/fill 空值 200（N-13）", st == 200)
    st, body = req("/api/load_config", method="POST",
                   data=json.dumps({"path": "config_step1.txt"}).encode())
    fields = json.loads(body).get("fields", {})
    ck("fill 空值后键被移除（N-13）", st == 200 and "xkey" not in fields)
    # pick_file 改 POST（排队③，2026-09-26）：GET 不再触发文件对话框（404），
    # 防跨站对话框轰炸；POST 弹框路径不做自动化断言（会挂起等待人工选择）
    st, _ = req("/api/pick_file?kind=file")
    ck("pick_file GET 404（不再支持 GET）", st == 404)
    # 结果浏览模式三端点（v1.0）：不存在的目录拒绝、合法目录激活、清除复位
    st, body = req("/api/set_browse", method="POST",
                   data=json.dumps({"path": os.path.join(CLIENT, "no_such_dir")}).encode())
    ck("set_browse 不存在的目录返回 ok=false", st == 200 and not json.loads(body).get("ok"))
    st, body = req("/api/set_browse", method="POST",
                   data=json.dumps({"path": CLIENT}).encode())
    ck("set_browse 合法目录激活", st == 200 and json.loads(body).get("ok"))
    st, body = req("/api/browse_state")
    ck("browse_state 报告激活状态", st == 200 and json.loads(body).get("active"))
    st, body = req("/api/clear_browse", method="POST", data=b"{}")
    ck("clear_browse 复位", st == 200)
    st, body = req("/api/browse_state")
    ck("清除后 browse_state 不激活", st == 200 and not json.loads(body).get("active"))
    # 结果浏览模式端到端（v1.0 本机自测）：合成"他机项目夹"（config 加
    # unlabeled MGF 加 result_* 报告），set_browse 后谱清单、单谱峰、报告页
    # 两分法摘要都从该夹读出，不依赖本机当前结果
    bproj = os.path.join(CLIENT, "other_machine_proj")
    os.makedirs(os.path.join(bproj, "unlabeled"), exist_ok=True)
    os.makedirs(os.path.join(bproj, "result_20260926120000", "1.unlabeled_result"),
                exist_ok=True)
    open(os.path.join(bproj, "config_step1.txt"), "w", encoding="utf-8").write(
        "label_type\tunlabeled\n")
    open(os.path.join(bproj, "unlabeled", "raw9_HCDFT.mgf"), "w", encoding="utf-8").write(
        "BEGIN IONS\nTITLE=raw9.1.1.2.0.dta\nPEPMASS=500.25 100\nCHARGE=2+\n"
        "100.1 500\n200.2 800\nEND IONS\n")
    hdr = "File_Name\tcandidates\tcandidate_PSM_score_by\tdecision_A\tvote_says_wrong\tchimeric\n"
    open(os.path.join(bproj, "result_20260926120000", "1.unlabeled_result", "report_A.tsv"),
         "w", encoding="utf-8").write(
        hdr + "raw9.1.1.2.0.dta\tPEPTIDEK()\t90\tcandidate_misidentification\t0\t0\n")
    open(os.path.join(bproj, "result_20260926120000", "1.unlabeled_result", "classification_A.tsv"),
         "w", encoding="utf-8").write(
        hdr + "raw9.1.1.2.0.dta\tPEPTIDEK()\t90\tcandidate_misidentification\t0\t0\n")
    open(os.path.join(bproj, "result_20260926120000", "1.unlabeled_result",
                      "two_part_report_summary.txt"), "w", encoding="utf-8").write(
        "part\trows\theadline\n"
        "1.annotation_correction_report\t1\tannotation corrections: 0\n"
        "2.novel_peptide_competition_report\t1\tnovel competition rows: 1 (credible 0, doubted 0)\n")
    st, body = req("/api/set_browse", method="POST",
                   data=json.dumps({"path": bproj}).encode())
    ck("浏览模式 E2E：set_browse 合成他机项目夹成功",
       st == 200 and json.loads(body).get("ok"))
    st, body = req("/api/spectra_list")
    ck("浏览模式 E2E：谱清单从该夹 report_A 读出（含 raw9 标题）",
       st == 200 and "raw9.1.1.2.0.dta" in body.decode("utf-8", "replace"))
    st, body = req("/api/spectrum?file=raw9.1.1.2.0.dta")
    ck("浏览模式 E2E：单谱峰从该夹 unlabeled MGF 读出（含 100.1）",
       st == 200 and b"100.1" in body)
    st, body = req("/api/report")
    ck("浏览模式 E2E：报告页两分法摘要从该夹读出（novel competition rows）",
       st == 200 and "novel competition rows" in body.decode("utf-8", "replace"))
    st, _ = req("/api/clear_browse", method="POST", data=b"{}")
    ck("浏览模式 E2E：clear 后复位", st == 200)
    # /api/stop 正常
    st, body = req("/api/stop", method="POST", data=b"{}")
    ck("/api/stop 200", st == 200)

    print("=== web_smoke:", len([f for f in fails]), "FAIL ===")
    if fails:
        sys.exit(1)
finally:
    srv.terminate()
    try:
        srv.wait(timeout=8)
    except Exception:
        srv.kill()
    # P3-3（2026-09-26，deepseek v1.0 审查）：临时客户端根随冒烟清理，
    # 不再每次泄漏一个目录
    shutil.rmtree(CLIENT, ignore_errors=True)
