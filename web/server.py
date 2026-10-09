# -*- coding: utf-8 -*-
"""web_shell/server.py - 网页版壳子的本地服务器（Python 标准库，无第三方依赖）。
The local server of the web shell (Python stdlib only, no third-party deps).

功能：① 生成 config_step1.txt（表单路径一律写绝对路径；写入用户指定的
结果保存文件夹）；② 一键运行总程序 user_run.R（等同命令行，日志实时可见）；
③ 停点续接（补写 config_step1/2/3 键）；④ 停止当前运行。
命令行格式完全保留：本壳子的"运行"就是执行
  Rscript code/user_run.R  （CONFIG=结果文件夹/config_step1.txt）
Endpoints: serve the form, generate the step-1 config into the user-chosen
result folder, run/stop the master program with a live log, and patch the
step-1/2/3 configs at stop points.

v1.0（2026-10-04）：网页界面语言切换（默认英文，侧栏按钮切中文，localStorage 记忆）。
前端三层机制：配对 span 显隐＋data-zh/data-en＋L(中文, English)；服务端消息保持
"中文
english" 双语对；R 管线日志不动。本文件仅随版升级，无逻辑改动。
v1.0 (2026-10-04): web UI language toggle (English by default, sidebar
button switches to Chinese, remembered in localStorage). Frontend uses three
layers (paired spans + data-zh/data-en + L(zh, en)); server messages stay
bilingual pairs; R pipeline logs unchanged. No logic changes in this file.
"""
import html
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.parse
from collections import Counter, defaultdict, deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HOST, PORT = "127.0.0.1", 8780
# v1.0（2026-10-03，实弹坐实）：不再接受 argv 传包根——bat 的 "%~dp0"
# 自带尾反斜杠，"G:\...\isoPSV_multiengine\" 的 \" 在 Windows 命令行解析
# 中＝转义引号，子进程实际收到尾带引号的假路径（本地复现：
# argv[1]='…isoPSV_multiengine"'）→ CLIENT_ROOT 全链失效（版本 unknown、
# 日志/配置路径错位）。脚本必然位于 <包根>\web\，由 __file__ 反推最稳。
# v1.0 (2026-10-03, confirmed live): the package root must NOT arrive via
# argv - the bat passes "%~dp0" whose trailing backslash escapes the closing
# quote, so the child receives a path ENDING IN A QUOTE (reproduced locally:
# argv[1]='...isoPSV_multiengine"') and every CLIENT_ROOT-derived path broke
# (version "unknown", misplaced logs/configs). The script always lives at
# <root>\web\, so derive the root from __file__ instead.
CLIENT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CODE_DIR = os.path.join(CLIENT_ROOT, "code")


def _probe_rscript():
    """R 解释器探测（P1-1，2026-09-25）：环境变量 → PATH → 注册表 →
    常见安装目录枚举；返回路径或 ""。旧实现硬编码 D 盘特定版本，是客户端
    部署最高频故障点。
    R interpreter probe (P1-1, 2026-09-25): environment variable -> PATH ->
    registry -> common install directories; returns a path or "". The old
    implementation hardcoded a D-drive version, the most frequent client
    deployment failure."""
    import shutil as _sh
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    w = _sh.which(exe)
    if w:
        return w
    if os.name == "nt":
        # F-2（2026-09-26）：按版本号降序取最新（旧实现枚举到第一个存在的
        # 注册表项即返回，且 glob 分支纯字典序——"R-4.10.0" < "R-4.6.1"，
        # 机器上多版本共存时会取到旧版）。
        # F-2 (2026-09-26): prefer the newest by version number (the old code
        # returned the first existing registry entry, and the glob branch was
        # plain lexicographic order - "R-4.10.0" < "R-4.6.1", so it could pick
        # an older version when several coexist).
        def _rver(p):
            m = re.search(r"R-(\d+)\.(\d+)(?:\.(\d+))?", p)
            if m:
                return tuple(int(x) for x in (m.group(1), m.group(2), m.group(3) or 0))
            return (0, 0, 0)
        cands = []
        try:
            import winreg as _reg
            for hive in (_reg.HKEY_LOCAL_MACHINE, _reg.HKEY_CURRENT_USER):
                for sub in (r"SOFTWARE\R-core\R", r"SOFTWARE\WOW6432Node\R-core\R"):
                    try:
                        k = _reg.OpenKey(hive, sub)
                        for i in range(_reg.QueryInfoKey(k)[0]):
                            ver = _reg.EnumKey(k, i)
                            try:
                                vk = _reg.OpenKey(k, ver)
                                p = _reg.QueryValueEx(vk, "InstallPath")[0]
                                cand = os.path.join(p, "bin", "x64", exe)
                                if os.path.exists(cand):
                                    cands.append(cand)
                                    continue
                                cand = os.path.join(p, "bin", exe)
                                if os.path.exists(cand):
                                    cands.append(cand)
                            except OSError:
                                pass
                    except OSError:
                        pass
        except ImportError:
            pass
        if cands:
            return sorted(cands, key=_rver, reverse=True)[0]
        import glob as _gb
        for root in ("C:/software", "D:/Bioinformatics", "C:/Program Files"):
            for pat in (os.path.join(root, "R-*", "bin", "x64", exe),
                        os.path.join(root, "R-*", "bin", exe),
                        os.path.join(root, "R", "R-*", "bin", "x64", exe),
                        os.path.join(root, "R", "R-*", "bin", exe)):
                hits = sorted(_gb.glob(pat), key=_rver, reverse=True)
                if hits:
                    return hits[0]
    return ""


RSCRIPT = os.environ.get("RSCRIPT") or _probe_rscript()
LOG_FILE = os.path.join(CLIENT_ROOT, "web_run.log")
INDEX = os.path.join(os.path.dirname(__file__), "index.html")


def _result_base(base):
    """项目目录 → 结果文件夹（pGlycoQuant 风格 result_YYYYMMDDHHMMSS，取最新）。
    浏览模式下改从 BROWSE_DIR 推导。
    The project folder -> result folder (pGlycoQuant-style
    result_YYYYMMDDHHMMSS; the newest wins). In browse mode it derives from
    BROWSE_DIR instead."""
    if BROWSE_DIR:
        base = BROWSE_DIR
    if not os.path.isdir(base):
        return base
    cands = [d for d in os.listdir(base)
             if os.path.isdir(os.path.join(base, d))
             and re.fullmatch(r"result_\d{14}", d)]
    return os.path.join(base, max(cands)) if cands else base


def _set_browse(path):
    """设置结果浏览目录：接受项目目录（含 result_*/unlabeled）或裸 result_* 目录。
    返回 (ok, msg)。空串清除浏览模式。"""
    global BROWSE_DIR, _BROWSE_MGF
    p = (path or "").strip()
    if not p:
        with _BROWSE_LOCK:
            BROWSE_DIR = None
            _BROWSE_MGF = None
        return True, ""
    if not os.path.isdir(p):
        return False, "目录不存在 not a directory: %s" % p
    rp = os.path.realpath(p)
    # 无标记谱目录：项目根（含 unlabeled/）就在选定目录或上一层
    mgf = ""
    for root in (rp, os.path.dirname(rp)):
        u = os.path.join(root, "unlabeled")
        if os.path.isdir(u):
            mgf = u
            break
    with _BROWSE_LOCK:
        BROWSE_DIR = rp
        _BROWSE_MGF = mgf or None
    return True, rp


def _browse_ctx():
    """浏览模式上下文 (pbase, rbase, label)；非浏览模式返回 None。
    label 优先取 rbase/2.labeled_result 内真实存在的 paired_report_<label>
    文件名（排序取定），项目根 config_step1.txt 的 label_type 作回退
    （UTF-8 优先、GBK 回退，与 _read_cfg 同口径），缺省 unlabeled。"""
    with _BROWSE_LOCK:
        bdir, bmgf = BROWSE_DIR, _BROWSE_MGF
    if not bdir:
        return None
    rbase = _result_base(bdir)
    proot = bdir
    if not os.path.isdir(os.path.join(proot, "unlabeled")):
        up = os.path.dirname(proot)
        if os.path.isdir(os.path.join(up, "unlabeled")):
            proot = up
    label = ""
    pr_dir = os.path.join(rbase, "2.labeled_result")
    if os.path.isdir(pr_dir):
        labs = sorted(re.match(r"paired_report_(.+)\.tsv$", f).group(1).lower()
                      for f in os.listdir(pr_dir)
                      if re.match(r"paired_report_.+\.tsv$", f))
        if labs:
            label = labs[0]
    if not label:
        cfg_path = os.path.join(proot, "config_step1.txt")
        if os.path.exists(cfg_path):
            try:
                text, _enc = read_config_text(cfg_path)
            except OSError:
                text = None
            if text:
                for ln in text.splitlines():
                    ln = ln.strip()
                    if ln.startswith("label_type"):
                        kv = ln.split("\t", 1)
                        if len(kv) == 2:
                            label = kv[1].strip().lower()
    return {"pbase": proot, "rbase": rbase, "label": label or "unlabeled",
            "mgf_dir": bmgf}

proc = None               # current running pipeline subprocess
_PROC_LOCK = threading.Lock()   # P1-6（2026-09-25）：串行化 run/stop 对 proc 的读写
RUNS = []                 # run history: list of dicts (started, ended, exit_code)
CURRENT_CFG = None        # path of the config_step1.txt the next run uses
CURRENT_RESULT = None     # result folder of the current config (its parent)
# 结果浏览模式（2026-09-26）：只读展示"另一台电脑跑完"的结果文件夹——
# 用户经 /api/set_browse 选定项目目录或 result_* 目录后，报告/谱图/PDF 全部
# 改从该目录读取，不再要求本机重跑。
BROWSE_DIR = None         # 选定的浏览目录（realpath）
_BROWSE_MGF = None        # 浏览模式的无标记谱目录（浏览目录同层 unlabeled/，可无）
_BROWSE_LOCK = threading.Lock()   # 自查整改（2026-09-26）：浏览状态整体换写，避免并发读到新旧混搭
_EXTRA_PDF_ROOTS = []     # /api/pdfs 手填目录的白名单登记（与 /pdf 服务共用，防"能扫不能看"）
_PS_PROBE = {"done": False, "result": ""}
def _ps_probe():
    """v1.0：一次性探测 PowerShell 可用性（页脚展示 picker: ok/failed）——
    弹窗问题的第一道远程诊断。v1.0: one-shot PowerShell liveness probe
    (footer shows picker: ok/failed) - the first line of remote diagnosis for
    the dialog issues."""
    if not _PS_PROBE["done"]:
        _PS_PROBE["done"] = True
        if os.name != "nt":
            _PS_PROBE["result"] = "n/a (non-windows)"
        else:
            import time as _t
            t0 = _t.time()
            try:
                r = subprocess.run(["powershell", "-NoProfile", "-Command",
                                    "Write-Output ok"], capture_output=True,
                                   text=True, timeout=20, encoding="utf-8",
                                   errors="replace")
                ms = int((_t.time() - t0) * 1000)
                _PS_PROBE["result"] = ("ok (%d ms)" % ms if
                                       (r.returncode == 0 and "ok" in (r.stdout or ""))
                                       else "failed: " + ((r.stderr or "").strip()[:120]
                                                          or ("rc=%s" % r.returncode)))
            except Exception as e:  # noqa
                _PS_PROBE["result"] = "failed: %s" % e
    return _PS_PROBE["result"]

_SRV_STARTED = time.strftime("%Y-%m-%d %H:%M:%S")
_TL_CACHE = {"size": -1, "mtime": -1, "timeline": None}  # timeline parse cache keyed by log size + mtime
_INDEX_BYTES = None       # index.html content cache
_INDEX_STAT = None        # (mtime, size) of the cached copy (reload on change)
_TK_LOCK = threading.Lock()  # serialize the native dialogs across HTTP threads


# v1.0（2026-10-03，v84 实弹定性：PS Forms 在客户端被无声悬挂——请求永
# 久 pending、无报错、探针 ok）——选择器改**进程内 ctypes 直调 Win32 对话
# 框 API**：comdlg32.GetOpenFileNameW（文件/多选）＋shell32.SHBrowseFor-
# FolderW（文件夹）。零子进程、零 .NET、零 COM——杀软/EDR 无拦截模式可盯。
# 作为首选路线；PS Forms/COM 降为后备（超时 60s，悬挂环境会快速越过）。
# v1.0 (2026-10-03, live v84: PS Forms hangs silently on the client -
# request pending forever, no error, probe ok) - the picker now calls the
# Win32 dialog APIs IN-PROCESS via ctypes: comdlg32.GetOpenFileNameW (files)
# and shell32.SHBrowseForFolderW (folders). No child process, no .NET, no
# COM - nothing for AV/EDR to pattern-match. This is the primary route;
# the PowerShell routes stay as fallbacks (60 s timeout each so hung
# environments skip them quickly).
# v1.0（2026-10-03，v84 实弹定性：PS Forms 在客户端被无声悬挂——请求永
# 久 pending、无报错、探针 ok）——选择器首选路线改**进程内 ctypes 直调
# Win32 对话框 API**：comdlg32.GetOpenFileNameW（文件/多选）＋shell32.
# SHBrowseForFolderW（文件夹）。零子进程、零 .NET、零 COM——杀软/EDR 无拦
# 截模式可盯。PS Forms/COM 降为后备（各 60s 超时，悬挂环境快速越过）。
# v1.0 (2026-10-03, live v84: PS Forms hangs silently on the client -
# request pending forever, no error, probe ok) - the PRIMARY picker route now
# calls the Win32 dialog APIs IN-PROCESS via ctypes: comdlg32.
# GetOpenFileNameW (files/multi) and shell32.SHBrowseForFolderW (folders).
# No child process, no .NET, no COM - nothing for AV/EDR to pattern-match.
# The PowerShell routes stay as fallbacks (60 s timeout each).
_PICK_ROUTE = [""]

def _pick_path_ctypes(kind):
    """进程内 Win32 对话框（ctypes）。返回 (ok, path_or_err)。
    In-process Win32 dialogs via ctypes. Returns (ok, path_or_err)."""
    import ctypes
    from ctypes import wintypes
    try:
        # v1.0（2026-10-03，v85 实弹终局）：用**控制台窗口**作对话框属主——
        # 服务器跑在 start_isoPSV.bat 的黑色 cmd 窗口里，这个窗口用户看得见、
        # 点得到，Windows 一定把对话框放在属主上面。此前 hwndOwner=None／离屏
        # 小窗／主屏居中窗体在该客户端全部不生效（窗口管理器忽略无前台属主的
        # 对话框），控制台窗口是这台机器上唯一确定可见的窗口句柄。
        # v1.0 (2026-10-03, live v85 final): use the CONSOLE WINDOW as the
        # dialog owner - the server runs in the black cmd window from
        # start_isoPSV.bat which the user can see and click; Windows places
        # owned dialogs on top of their owner. All previous owner strategies
        # (None / off-screen / primary-centered) failed on this client because
        # the window manager ignores dialogs without a foreground owner. The
        # console window is the one guaranteed-visible HWND on that machine.
        kernel32 = ctypes.windll.kernel32
        user32 = ctypes.windll.user32
        _hwnd = kernel32.GetConsoleWindow()
        if _hwnd:
            user32.SetForegroundWindow(_hwnd)

        if kind == "dir":
            shell32 = ctypes.windll.shell32
            ole32 = ctypes.windll.ole32

            class BROWSEINFOW(ctypes.Structure):
                _fields_ = [("hwndOwner", wintypes.HWND),
                            ("pidlRoot", ctypes.c_void_p),
                            ("pszDisplayName", wintypes.LPWSTR),
                            ("lpszTitle", wintypes.LPCWSTR),
                            ("ulFlags", wintypes.UINT),
                            ("lpfn", ctypes.c_void_p),
                            ("lParam", wintypes.LPARAM),
                            ("iImage", ctypes.c_int)]
            buf = ctypes.create_unicode_buffer(4096)
            bi = BROWSEINFOW()
            bi.hwndOwner = _hwnd
            bi.pszDisplayName = ctypes.cast(buf, wintypes.LPWSTR)
            bi.lpszTitle = "Select folder / 选择文件夹"
            bi.ulFlags = 0x1 | 0x40  # BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE
            shell32.SHBrowseForFolderW.restype = ctypes.c_void_p
            shell32.SHBrowseForFolderW.argtypes = [ctypes.POINTER(BROWSEINFOW)]
            shell32.SHGetPathFromIDListW.argtypes = [ctypes.c_void_p, wintypes.LPWSTR]
            shell32.SHGetPathFromIDListW.restype = wintypes.BOOL
            ole32.CoTaskMemFree.argtypes = [ctypes.c_void_p]
            pidl = shell32.SHBrowseForFolderW(ctypes.byref(bi))
            if not pidl:
                return False, "cancelled"
            path = ctypes.create_unicode_buffer(4096)
            ok = shell32.SHGetPathFromIDListW(pidl, path)
            ole32.CoTaskMemFree(pidl)
            if ok and path.value:
                return True, path.value
            return False, "no filesystem path"

        comdlg32 = ctypes.windll.comdlg32

        class OPENFILENAMEW(ctypes.Structure):
            _fields_ = [("lStructSize", wintypes.DWORD),
                        ("hwndOwner", wintypes.HWND),
                        ("hInstance", wintypes.HINSTANCE),
                        ("lpstrFilter", wintypes.LPCWSTR),
                        ("lpstrCustomFilter", wintypes.LPWSTR),
                        ("nMaxCustFilter", wintypes.DWORD),
                        ("nFilterIndex", wintypes.DWORD),
                        ("lpstrFile", wintypes.LPWSTR),
                        ("nMaxFile", wintypes.DWORD),
                        ("lpstrFileTitle", wintypes.LPWSTR),
                        ("nMaxFileTitle", wintypes.DWORD),
                        ("lpstrInitialDir", wintypes.LPCWSTR),
                        ("lpstrTitle", wintypes.LPCWSTR),
                        ("Flags", wintypes.DWORD),
                        ("nFileOffset", wintypes.WORD),
                        ("nFileExtension", wintypes.WORD),
                        ("lpstrDefExt", wintypes.LPCWSTR),
                        ("lCustData", wintypes.LPARAM),
                        ("lpfnHook", ctypes.c_void_p),
                        ("lpTemplateName", wintypes.LPCWSTR),
                        ("pvReserved", ctypes.c_void_p),
                        ("dwReserved", wintypes.DWORD),
                        ("FlagsEx", wintypes.DWORD)]
        nul = chr(0)
        flt = ctypes.create_unicode_buffer(
            "".join(["All files", nul, "*.*", nul, nul, nul]))
        fbuf = ctypes.create_unicode_buffer(65536)
        ofn = OPENFILENAMEW()
        ofn.lStructSize = ctypes.sizeof(OPENFILENAMEW)
        ofn.hwndOwner = _hwnd
        ofn.lpstrFilter = ctypes.cast(flt, wintypes.LPCWSTR)
        ofn.lpstrFile = ctypes.cast(fbuf, wintypes.LPWSTR)
        ofn.nMaxFile = len(fbuf)
        ofn.lpstrTitle = "Select file(s) / 选择文件"
        ofn.Flags = 0x80000 | 0x1000 | 0x800 | 0x8  # EXPLORER|FILEMUSTEXIST|PATHMUSTEXIST|NOCHANGEDIR
        if kind == "files":
            ofn.Flags |= 0x200  # ALLOWMULTISELECT
        if not comdlg32.GetOpenFileNameW(ctypes.byref(ofn)):
            return False, "cancelled"
        parts = [s for s in fbuf.raw.split(nul) if s]
        if kind == "files" and len(parts) > 1:
            d = parts[0]
            return True, chr(10).join(os.path.join(d, x) for x in parts[1:])
        return True, parts[0]
    except Exception as e:  # noqa
        return False, "ctypes: %s" % e


def _pick_path_ps(kind):
    """v1.0（2026-10-03，v81 实弹：进程已新版仍无弹窗）：三层强化——
    (1) 属主窗体改"离屏小窗"（Minimized 在部分环境拖累对话框前置）；
    (2) 文件夹加 Shell.Application COM 后备（不走 .NET Forms，策略受限机器
        也能弹）；(3) 失败返回 stderr 诊断（此前前端静默吞掉，"没弹窗也没
        报错"）。返回 (ok, path_or_err)。
    v1.0 (2026-10-03, live v81 - new process confirmed, still no dialog):
    (1) owner switched to an off-screen tiny window (Minimized hurt
    foregrounding on some setups); (2) folders fall back to the
    Shell.Application COM browser (no .NET Forms - works under restrictive
    policies); (3) failures carry stderr diagnostics (the frontend used to
    swallow them). Returns (ok, path_or_err)."""
    if os.name != "nt":
        return False, "non-windows"
    # 离屏小窗属主：TopMost＋屏幕外 1×1（Minimized 的 Activate 在部分环境
    # 无效，对话框跟着沉底）
    # off-screen tiny owner: TopMost + 1x1 at (-32000,-32000) (a Minimized
    # owner fails to foreground the dialog on some setups)
    # v1.0（2026-10-03，v83 实弹）：-32000 离屏定位在"副屏位于左侧负坐
    # 标"的多显示器布局上会落在真实屏幕——对话框跟着跑到不可见屏。改为定
    # 位在主屏工作区中心（TopMost 属主），任何布局下对话框都出现在主屏。
    # v1.0 (2026-10-03, live v83): the -32000 off-screen position lands on
    # a REAL monitor when a secondary display sits at negative coordinates -
    # the dialog followed it off-screen. Center the owner on the primary
    # screen's working area instead so the dialog always shows there.
    owner = ("$b=[System.Windows.Forms.Screen]::PrimaryScreen.WorkingArea;"
             "$w = New-Object Windows.Forms.Form;"
             "$w.TopMost=$true; $w.ShowInTaskbar=$false;"
             "$w.StartPosition='Manual';"
             "$w.Location=New-Object System.Drawing.Point(([int]($b.Width/2)),([int]($b.Height/2)));"
             "$w.Size=New-Object System.Drawing.Size(1,1);"
             "$w.Show();")
    if kind == "dir":
        ps_forms = ("[Console]::OutputEncoding=[Text.Encoding]::UTF8;"
                    "Add-Type -AssemblyName System.Windows.Forms;" + owner +
                    "$f = New-Object Windows.Forms.FolderBrowserDialog;"
                    "$f.Description = 'Select folder / 选择文件夹';"
                    "$r = $f.ShowDialog($w); $w.Close();"
                    "if($r -eq 'OK'){$f.SelectedPath}")
        # COM 后备：Shell.Application 不依赖 .NET Forms（Forms 被策略禁用时可用）
        # COM fallback: Shell.Application needs no .NET Forms (works where
        # Forms is blocked by policy)
        ps_com = ("[Console]::OutputEncoding=[Text.Encoding]::UTF8;"
                  "$sh = New-Object -ComObject Shell.Application;"
                  "$f = $sh.BrowseForFolder(0, 'Select folder / 选择文件夹', 0x41, 0);"
                  "if($f){$f.Self.Path}")
        attempts = [ps_forms, ps_com]
    elif kind == "files":
        attempts = [("[Console]::OutputEncoding=[Text.Encoding]::UTF8;"
                     "Add-Type -AssemblyName System.Windows.Forms;" + owner +
                     "$f = New-Object Windows.Forms.OpenFileDialog;"
                     "$f.Multiselect = $true;"
                     "$f.Title = 'Select files / 选择文件（可多选）';"
                     "$r = $f.ShowDialog($w); $w.Close();"
                     "if($r -eq 'OK'){$f.FileNames -join \"`n\"}")]
    else:
        attempts = [("[Console]::OutputEncoding=[Text.Encoding]::UTF8;"
                     "Add-Type -AssemblyName System.Windows.Forms;" + owner +
                     "$f = New-Object Windows.Forms.OpenFileDialog;"
                     "$f.Title = 'Select file / 选择文件';"
                     "$r = $f.ShowDialog($w); $w.Close();"
                     "if($r -eq 'OK'){$f.FileName}")]
    last_err = ""
    for ps in attempts:
        try:
            out = subprocess.run(
                ["powershell", "-NoProfile", "-STA", "-Command", ps],
                capture_output=True, text=True, timeout=300,
                encoding="utf-8", errors="replace")
            p = (out.stdout or "").strip()
            if out.returncode == 0 and p:
                return True, p
            last_err = ((out.stderr or "").strip() or
                        ("rc=%s empty-output" % out.returncode))
        except Exception as e:  # noqa
            last_err = "%s: %s" % (type(e).__name__, e)
    return False, last_err[:400]


def esc(s):
    """HTML 转义（把路径/配置值/日志文本插入报告页前使用）。
    HTML-escape a string before inserting it into report HTML."""
    return html.escape(str(s), quote=True)


def read_config_text(path):
    """读取参数文件文本与原编码：UTF-8（含 BOM）优先，失败退 GBK。
    返回 (文本, 编码) 或 (None, None)。
    Read a parameter file as text with its encoding: UTF-8 (BOM-tolerant) first,
    GBK as the fallback. Returns (text, encoding) or (None, None)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    for enc in ("utf-8-sig", "gbk"):
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return None, None


def get_index_bytes():
    """返回 index.html 字节内容；文件修改后自动重读（无需重启服务器）。
    Return the index.html bytes; reload when the file changes (no server
    restart needed)."""
    global _INDEX_BYTES, _INDEX_STAT
    try:
        st = os.stat(INDEX)
        key = (st.st_mtime_ns, st.st_size)
    except OSError:
        key = None
    if _INDEX_BYTES is None or key != _INDEX_STAT:
        with open(INDEX, "rb") as fh:
            _INDEX_BYTES = fh.read()
        _INDEX_STAT = key
    return _INDEX_BYTES


def parse_timeline(text):
    """解析运行日志为模块时间线：
    模块开始/完成横幅（带 @ HH:MM:SS 时间戳，中文行与英文行成对）、停点、失败，
    以及每个模块内的阶段计时行（[计时] tag：已用 Xs / CPU Ys 与 [timing] 英文行）。
    模块可嵌套（如 pUniFind 自动 rescore/denovo 嵌在用户第二步内）：用栈维护
    打开中的模块，结束/停点/失败横幅按标签回填到对应模块，计时行归当前栈顶。
    Parse the run log into a module timeline: module start/done banners with
    @ HH:MM:SS timestamps (the Chinese and English lines are paired), stop
    points, failures, and per-module phase timing lines ([计时] tag: elapsed Xs
    / CPU Ys plus the [timing] English line). Modules may nest (the pUniFind
    auto rescore/denovo runs nest inside user step 2): a stack tracks the open
    modules, done/stop/failed banners are matched back to their module by tag,
    and timing lines attach to the stack top."""
    mods = []       # 所有模块，按首次出现顺序排列
    stack = []      # 当前打开中的模块（栈顶＝接收阶段计时行）
    last_phase = None

    def close_by(tag, status, end_ts):
        """按中文标签关闭对应模块（自栈顶向下找；找不到则关栈顶以兼容旧日志）。
        Close the module whose Chinese tag matches, searching from the top of the
        stack down; fall back to the top module for older logs without the tag."""
        for i in range(len(stack) - 1, -1, -1):
            if stack[i]["name_zh"] == tag:
                m = stack.pop(i)
                m["status"] = status
                m["end"] = end_ts
                return
        if stack:
            stack[-1]["status"] = status
            stack[-1]["end"] = end_ts

    for line in text.splitlines():
        m = re.search(r"==== \[总程序 / Master\] (.+?) 开始 @ (\d\d:\d\d:\d\d) ====", line)
        if m:
            mod = {"tag": m.group(1), "name_zh": m.group(1), "name_en": m.group(1),
                   "start": m.group(2), "end": "", "status": "running", "phases": []}
            mods.append(mod)
            stack.append(mod)
            last_phase = None
            continue
        m = re.search(r"==== \[Master\] (.+?) start @ (\d\d:\d\d:\d\d) ====", line)
        if m and stack:
            stack[-1]["name_en"] = m.group(1)
            continue
        m = re.search(r"==== \[总程序 / Master\] (.+?) 完成 @ (\d\d:\d\d:\d\d) ====", line)
        if m:
            close_by(m.group(1), "done", m.group(2))
            last_phase = None
            continue
        m = re.search(r"==== \[总程序 / Master\] (.+?) 停点（等待用户输入）@ (\d\d:\d\d:\d\d) ====", line)
        if m:
            close_by(m.group(1), "stop", m.group(2))
            last_phase = None
            continue
        m = re.search(r"\[总程序 / Master\] (.+?) 失败 \(exit (\d+)\)", line)
        if m:
            close_by(m.group(1), "failed", time.strftime("%H:%M:%S"))
            last_phase = None
            continue
        m = re.search(r"\[计时\] (.+?)：已用 ([\d.]+)s / CPU ([\d.]+)s", line)
        if m and stack:
            ph = {"tag": m.group(1), "name_zh": m.group(1), "name_en": m.group(1),
                  "elapsed": float(m.group(2)), "cpu": float(m.group(3))}
            stack[-1]["phases"].append(ph)
            last_phase = ph
            continue
        m = re.search(r"\[timing\] (.+?): elapsed ([\d.]+)s / CPU ([\d.]+)s", line)
        if m and last_phase is not None:
            last_phase["name_en"] = m.group(1)
    return mods


def current_timeline():
    """解析当前日志（按文件大小缓存），并回填最近一次运行的时间线。
    Parse the current log (cached by file size) and back-fill the timeline of
    the latest run."""
    if not os.path.exists(LOG_FILE):
        return None
    st = os.stat(LOG_FILE)
    if (_TL_CACHE["size"] == st.st_size and _TL_CACHE["mtime"] == st.st_mtime_ns
            and _TL_CACHE["timeline"] is not None):
        return _TL_CACHE["timeline"]
    with open(LOG_FILE, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    tl = parse_timeline(text)
    _TL_CACHE["size"] = st.st_size
    _TL_CACHE["mtime"] = st.st_mtime_ns
    _TL_CACHE["timeline"] = tl
    return tl


def build_step1_config(fields):
    # all_mgf 自动推导自 mgf（阶段 B 的原始 MGF 文件夹＝输入 mgf 所在文件夹）：
    # mgf 填单个文件夹 → 直接用它；填多个文件且同在一个文件夹 → 用该文件夹；
    # 多个文件分散在不同文件夹 → 留空（阶段 B 停点手填）
    # all_mgf is derived from mgf (the stage-B original MGF folder = the folder
    # of the input mgfs): a single folder -> used as-is; multiple files in ONE
    # folder -> that folder; files spread over folders -> empty (stage B stops
    # and asks the user to fill it)
    mgf_raw = fields.get("mgf", "")
    entries = [ln.strip() for ln in str(mgf_raw).replace(";", "\n").splitlines() if ln.strip()]
    all_mgf = ""
    if len(entries) == 1 and os.path.isdir(entries[0]):
        all_mgf = entries[0]
    elif entries:
        dirs = [os.path.dirname(os.path.abspath(p)) for p in entries]
        if len(set(dirs)) == 1 and os.path.isdir(dirs[0]):
            all_mgf = dirs[0]
    # P0-8（2026-09-26，GLM 审查）：单值字段拒绝换行——否则字段值可注入
    # 任意键行（"mgf\t<值>" 拼接被换行拆成新键），构成任意内容写原语
    # v1.0 修复：label_type 自 v1.0 起是多行键（每行一个标记），此前
    # 被本名单误剥换行、多行值粘连成单行——移出名单，改在下方逐行净化
    #（换行是结构、制表符才是键注入载体，逐行剔除）。
    single = ("search_engine", "psm_file", "fasta", "enzyme",
              "auto_search", "pfind_exe", "use_punifind", "punifind_exe",
              "punifind_rescore", "punifind_denovo", "punifind_python",
              "use_pdeep2", "pdeep2_exe", "instrument", "pdf_on",
              "denovo_lookback", "mem_guard_mb",
              "unlabeled_channel",
              "candidate_scope", "proteogenomics_mode")
    for k in single:
        v = str(fields.get(k, "") or "")
        if "\n" in v or "\r" in v:
            fields[k] = v.replace("\r", "").replace("\n", "").strip()
    # v1.0：label_type 逐行净化——每行去空白与制表符（P0-8 同口径防键注入）
    lt_rows = [x.strip().replace("\t", " ")
               for x in str(fields.get("label_type", "") or "").replace("\r", "\n").split("\n")]
    lt_rows = [x for x in lt_rows if x] or ["unlabeled"]
    lines = [
        "# 网页版壳子生成的用户第一步参数文件",
        "# web-shell generated step-1 parameter file",
        "result_dir	" + fields.get("result_dir", ""),
        "search_engine\t" + fields.get("search_engine", "pfind"),
        "psm_file\t" + fields.get("psm_file", ""),
        "fasta\t" + fields.get("fasta", ""),
        # v1.0（2026-09-28 用户口径）：label_type 为多行键（每行一个标记：
        # N15/C13、旧预设名或 v1.0 自定义行 名称|位置|氨基酸|组成），
        # 空=无标记；网页表单按类 1+行合成多行值
        "label_type\t" + "\n".join(lt_rows),
        # v1.0（2026-09-26）：酶切类型下拉，默认 trypsin
        "enzyme\t" + (fields.get("enzyme", "") or "trypsin"),
        "mgf\t" + fields.get("mgf", ""),
        "all_mgf\t" + all_mgf,
        "auto_search\t" + fields.get("auto_search", "1"),
        "pfind_exe\t" + fields.get("pfind_exe", ""),
        # v1.0（2026-10-09 用户口径）：python_exe 键不再生成——重编号解释器
        # 由 R 侧探测（PATH → py -3 → python3 → 软件自带 runtime/python）；
        # 旧 config 手写的 python_exe 键仍被兼容读取
        # v1.0: the python_exe key is no longer written - the renumbering
        # interpreter is probed R-side (PATH -> py -3 -> python3 -> the bundled
        # runtime/python); legacy python_exe keys in old configs still work
        "use_punifind\t" + fields.get("use_punifind", "1"),
        "punifind_exe\t" + fields.get("punifind_exe", ""),
        "punifind_rescore\t" + fields.get("punifind_rescore", "0"),
        "punifind_denovo\t" + fields.get("punifind_denovo", "0"),
        "punifind_python\t" + fields.get("punifind_python", ""),
        "use_pdeep2\t" + fields.get("use_pdeep2", "1"),
        "pdeep2_exe\t" + fields.get("pdeep2_exe", ""),
        "instrument\t" + (fields.get("instrument", "") or "Lumos"),
        "pdf_on\t" + (fields.get("pdf_on", "") or "0"),
        "denovo_lookback\t" + (fields.get("denovo_lookback", "") or "keep"),
        # v1.0（2026-09-27 用户口径）：全候选竞争模式（忽略引擎 1% FDR）
        # 与工作流程门控（全量/纯蛋白质基因组学/已注释矫正）
        "candidate_scope\t" + (fields.get("candidate_scope", "") or "filtered"),
        "proteogenomics_mode\t" + (fields.get("proteogenomics_mode", "") or "0"),
        "",
    ]
    # F-9（2026-09-26）：内存守卫上限（MB）——留空不写键，user_run 沿用默认
    # （物理内存 80% 自适应）；非法值由 user_run 忽略。
    # F-9 (2026-09-26): memory guard ceiling in MB - empty keeps the default
    # (80% of physical RAM) in user_run; invalid values are ignored there.
    mem_guard = str(fields.get("mem_guard_mb", "") or "").strip()
    if mem_guard:
        try:
            v = float(mem_guard)
            if v > 0:
                lines.insert(-1, "ISOMEM_LIMIT_MB\t" + mem_guard)
        except ValueError:
            pass
    return "\n".join(lines)


def run_pipeline():
    global proc, CURRENT_CFG

    def _watch_and_copy_log(result_dir, start_time):
        """李治反馈①：进程退出（完成/报错）后把 web_run.log 拷贝到 result_dir
        作为带时间戳的运行日志。watcher 线程 daemon，server 退出时自动结束。
        tester feedback: copy web_run.log to result_dir as a timestamped
        run log after the pipeline exits (success or error)."""
        import threading, shutil as _sh
        try:
            _p = proc
            if _p is None:
                return
            _p.wait()
            if result_dir and os.path.isdir(result_dir):
                _ts = time.strftime("%Y%m%d_%H%M%S")
                _dst = os.path.join(result_dir,
                                    "run_log_%s.txt" % _ts)
                _sh.copy2(LOG_FILE, _dst)
        except Exception:
            pass

    # P1-6（2026-09-25）：锁内完成 检查→启动 全程，消除两个并发 /api/run
    # 同时通过 poll() 检查、各自 Popen 的 check-then-act 竞态（前一个句柄
    # 被覆盖后无法 stop，两个流程写同一结果目录互相污染）。
    with _PROC_LOCK:
        if proc is not None and proc.poll() is None:
            return False, "流程已在运行\nalready running"
        if proc is not None and proc.poll() is not None:
            try:
                proc.wait()          # 回收已退出进程，避免句柄泄漏
            except Exception:
                pass
            proc = None
        # 李治反馈②：CURRENT_CFG 为空时先找 CURRENT_RESULT（上次 generate 写入的
        # result_dir）下的 config_step1.txt，再退回包根
        # tester feedback: when CURRENT_CFG is unset, look in the last generated
        # result_dir before falling back to the package root
        cfg = CURRENT_CFG
        if not cfg or not os.path.exists(cfg):
            _rd_try = (CURRENT_RESULT or "").strip()
            if _rd_try:
                _cfg_try = os.path.join(_rd_try, "config_step1.txt")
                if os.path.isfile(_cfg_try):
                    cfg = _cfg_try
                    CURRENT_CFG = cfg
        if not cfg or not os.path.exists(cfg):
            cfg = os.path.join(CLIENT_ROOT, "config_step1.txt")
        if not os.path.exists(cfg):
            return False, ("参数文件不存在——请先到'参数汇总'页点'生成参数文件'\n"
                           "config file not found - go to Summary page and click Generate config\n"
                           "(looked in: " + (CURRENT_RESULT or CLIENT_ROOT) + ")")
        # 施佳辉反馈（2026-10-08）：result_dir 空或缺失时第一步报 not found
        # tester feedback: step 1 errors when result_dir is empty or missing
        _rd_for_log = ""
        try:
            for _ln in open(cfg, encoding="utf-8", errors="replace"):
                if _ln.startswith("result_dir\t"):
                    _rd = _ln.split("\t", 1)[1].strip()
                    _rd_for_log = _rd
                    if not _rd:
                        return False, ("result_dir 为空——请在 STEP 1 填写结果保存文件夹后点"
                                       "生成参数文件再运行\nresult_dir is empty - fill it in STEP 1, "
                                       "click Generate config, then Run")
                    if not os.path.isdir(_rd):
                        try:
                            os.makedirs(_rd, exist_ok=True)
                        except OSError as e:
                            return False, ("result_dir 无法创建：" + _rd + "\n" + str(e) +
                                           "\nresult_dir cannot be created: " + _rd)
                    break
            else:
                return False, ("config_step1.txt 缺少 result_dir 键——请点生成参数文件\n"
                               "result_dir key missing from config_step1.txt - click Generate config")
        except Exception as e:
            return False, "读取配置文件失败\nfailed to read config: %s" % e
        if not os.path.exists(RSCRIPT):
            # P2-8（2026-09-25 代码审查）：报错文案与自动探测同步
            return False, "未探测到 Rscript，请安装 R（≥4.1，推荐 4.5.x）或设置 RSCRIPT 环境变量指向 Rscript.exe\nRscript not found; install R (>=4.1, 4.5.x recommended) or set the RSCRIPT environment variable"
        env = dict(os.environ)
        env.pop("LC_CTYPE", None)
        env.pop("LC_ALL", None)
        env["CONFIG"] = cfg
        _TL_CACHE["size"] = -1          # log is truncated: reset the timeline cache
        _TL_CACHE["mtime"] = -1
        _TL_CACHE["timeline"] = None
        with open(LOG_FILE, "w", encoding="utf-8") as fh:
            fh.write("==== 网页版壳子 web shell 启动时间 started: %s ====\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            fh.flush()
            try:
                proc = subprocess.Popen(
                    [RSCRIPT, os.path.join(CODE_DIR, "user_run.R")],
                    cwd=CODE_DIR, env=env, stdout=fh, stderr=subprocess.STDOUT)
            except OSError as e:
                fh.write("启动失败 launch failed: %s\n" % e)
                proc = None
                return False, "启动失败\nlaunch failed: %s" % e
        RUNS.append({"started": time.strftime("%Y-%m-%d %H:%M:%S"),
                     "ended": "", "exit_code": None})
        RUNS[:] = RUNS[-20:]           # keep the run history bounded
        # 李治反馈①：watcher 线程等进程退出后拷贝日志到 result_dir
        import threading as _th
        _th.Thread(target=_watch_and_copy_log,
                   args=(_rd_for_log, RUNS[-1]["started"]), daemon=True).start()
        return True, "已启动\nstarted (pid %d)" % proc.pid


def pick_path(kind):
    """弹系统文件/文件夹选择框；Windows 走 PowerShell 原生对话框（_pick_path
    内 v1.0 路线），非 Windows 走 tkinter；失败时用户可直接在文本框输入。
    Open a native file/folder dialog; PowerShell native on Windows (the
    v1.0 route inside _pick_path), tkinter elsewhere; on failure the user
    can simply type the path into the box instead."""
    with _TK_LOCK:
        return _pick_path(kind)


# ---------------------------------------------------------------------------
# 一键环境检测（用户口径 2026-09-21：网页按钮 → 自动扫描 → 回填路径）
# one-click environment detection: web button -> scan installs -> autofill paths
# 与 code\sysinfo_probe.R 同一套检测口径（R 版为命令行独立脚本，此为网页版）
# ---------------------------------------------------------------------------
_DETECT_ROOTS = ["D:/Bioinformatics", "C:/Program Files", "C:/Program Files (x86)",
                 "C:/Bioinformatics", "C:/software", "E:/Bioinformatics", "F:/Bioinformatics"]

def _walk_bounded(root, max_depth=3, max_nodes=30000, deadline=None):
    """限深限量的目录遍历（跳过 reparse point，避免死循环与超时）。
    Depth- and node-bounded directory walk (reparse points skipped)."""
    n = 0
    stack = [(root, 0)]
    while stack and n < max_nodes:
        if deadline is not None and (n & 15) == 0 and time.time() > deadline:
            return            # v1.0：总时限到即停（AV 慢机防拖死）
        cur, d = stack.pop()
        n += 1
        try:
            it = sorted(os.scandir(cur), key=lambda e: e.name)
        except OSError:
            continue
        for e in it:
            try:
                if e.is_dir(follow_symlinks=False):
                    if d < max_depth:
                        stack.append((e.path, d + 1))
                else:
                    yield e.name.lower(), os.path.dirname(e.path), e.path
            except OSError:
                continue

def detect_environment(extra_roots=None):
    """扫描常见安装根目录定位管线依赖；返回 engines/system/verdict。
    extra_roots：用户在界面上追加的扫描根目录（v1.0，最多 4 个、须存在、
    去重后插队在最前）；None/空＝仅内置根目录（行为与旧版逐字节一致）。
    Scan the common install roots for pipeline dependencies.
    extra_roots: user-added scan roots from the UI (v1.0; at most 4, must
    exist, deduplicated and placed first); None/empty = built-in roots only."""
    eng = dict.fromkeys(
        ["pfind_exe", "pdeep2_exe", "punifind_exe", "punifind_python",
         "msfragger", "comet", "sage", "maxquant", "metamorpheus", "msconvert",
         "xtandem", "msbooster", "pdeep3", "alphapeptdeep", "predictor_python"], "")
    # v1.0：用户显式指定的扫描根（extra_roots）单独跟踪——这些目录
    # 走全量深搜（12 层/12 万节点/120 s 预算），即"输入路径的所有文件搜索"；
    # 内置根维持限深限量（深 Program Files 更紧）。
    # v1.0: user-specified roots (extra_roots) are tracked separately and
    # walked exhaustively (12 levels / 120k nodes / 120 s) - a full file search
    # under the given path; built-in roots keep their bounded walks.
    extras = []
    if extra_roots:
        for r in extra_roots[:4]:
            r = str(r).strip().strip('"')
            if r and os.path.isdir(r) and r.lower() not in [x.lower() for x in extras]:
                extras.append(r)
    # v1.0：用户填了扫描范围 → 只扫用户目录，跳过全部内置根
    # （客户端 AV 可挂死 os.scandir，内置大根在 AV 机上不可控——用户指哪扫哪）
    if extras:
        roots = list(extras)
    else:
        roots = []
        for r in _DETECT_ROOTS:
            if r.lower() not in [x.lower() for x in roots] and os.path.isdir(r):
                roots.append(r)
    t0 = time.time()
    for root in roots:
        rl = root.lower()
        if rl in [x.lower() for x in extras]:
            walk_kw = dict(max_depth=12, max_nodes=120000, deadline=t0 + 120)
        elif rl.startswith("c:/program files"):
            walk_kw = dict(max_depth=4, max_nodes=8000, deadline=t0 + 45)
        else:
            walk_kw = dict(max_depth=6, max_nodes=30000, deadline=t0 + 45)
        for nm, dname, fpath in _walk_bounded(root, **walk_kw):
            if not eng["pfind_exe"] and nm == "pfind.exe":
                eng["pfind_exe"] = dname
            elif not eng["pdeep2_exe"] and nm == "predict.py" and dname.lower().endswith("pdeep2"):
                # v99k: go up 2 levels to the pDeep2 ROOT (predict.py sits at
                # .../pDeep2/pDeep/pDeep2/predict.py; the parent dir is the
                # INNER pDeep2 — the root with py37/python37 is 2 levels up)
                _pd2root = os.path.dirname(os.path.dirname(dname))
                eng["pdeep2_exe"] = (_pd2root
                    if os.path.basename(_pd2root).lower().startswith("pdeep2")
                    else dname)
            elif not eng["punifind_exe"] and nm == "official_score_workflow.exe":
                eng["punifind_exe"] = dname
            elif not eng["punifind_python"] and nm == "torch" and dname.lower().endswith("site-packages"):
                eng["punifind_python"] = os.path.dirname(os.path.dirname(dname))
            elif not eng["msfragger"] and nm.startswith("msfragger") and nm.endswith(".jar"):
                eng["msfragger"] = fpath
            elif not eng["comet"] and nm.startswith("comet") and nm.endswith(".exe"):
                eng["comet"] = fpath
            elif not eng["sage"] and nm == "sage.exe":
                eng["sage"] = fpath
            elif not eng["maxquant"] and nm == "maxquantcmd.exe":
                eng["maxquant"] = fpath
            elif not eng["metamorpheus"] and nm == "cmd.exe" and "metamorpheus" in fpath.lower():
                eng["metamorpheus"] = fpath
            elif not eng["metamorpheus"] and nm == "metamorpheusgui.exe":
                eng["metamorpheus"] = fpath
            elif not eng["msconvert"] and nm == "msconvert.exe":
                eng["msconvert"] = fpath
            elif not eng["xtandem"] and nm == "tandem.exe":
                eng["xtandem"] = fpath
            elif not eng["msbooster"] and nm.endswith(".jar") and "msbooster" in nm:
                eng["msbooster"] = fpath
            # v1.0：pDeep3（仓库特征＝pDeep/cmd/tune_and_predict.py）与
            # AlphaPeptDeep（site-packages/peptdeep 包）——取仓库根
            elif not eng["pdeep3"] and nm == "tune_and_predict.py" and "pdeep" in dname.lower():
                d2 = os.path.dirname(dname)
                if os.path.basename(d2).lower() == "pdeep":
                    eng["pdeep3"] = os.path.dirname(d2)
            elif not eng["alphapeptdeep"] and nm == "peptdeep" and dname.lower().endswith("site-packages"):
                eng["alphapeptdeep"] = os.path.dirname(os.path.dirname(dname))
            # v1.0：附件包 python37 → predictor_python
            # 修：查父目录名（basename 含 python37）而非父目录内容——嵌套扫描
            # 时 dname=…/python37，basename 就是 python37，直接命中
            elif not eng["predictor_python"] and nm == "python.exe" and                     "python37" in os.path.basename(dname).lower():
                eng["predictor_python"] = os.path.join(dname, nm)
    # 定向兜底：torch 虚拟环境（walk 深度不够时直接探测 *venv* 目录）
    # targeted fallback: torch venv probed directly (below the walk depth)
    if not eng["punifind_python"]:
        import glob as _gb
        for root in roots:
            for cand in _gb.glob(os.path.join(root, "*venv*")) + \
                         _gb.glob(os.path.join(root, "*python*")):
                if os.path.isdir(os.path.join(cand, "Lib", "site-packages", "torch")):
                    eng["punifind_python"] = cand
                    break
            if eng["punifind_python"]:
                break
    # pfind_exe 约定为安装根（含 bin\pFind.exe 时回退一级）
    # pfind_exe convention: the install ROOT (step up from ...\bin)
    if eng["pfind_exe"].lower().endswith("\\bin") or eng["pfind_exe"].lower().endswith("/bin"):
        eng["pfind_exe"] = os.path.dirname(eng["pfind_exe"])
    # 兜底：PATH 上直接找
    import shutil as _sh
    for key, exe in (("sage", "sage.exe"), ("comet", "comet.win64.exe"),
                     ("msconvert", "msconvert.exe"), ("maxquant", "MaxQuantCmd.exe")):
        if not eng[key]:
            w = _sh.which(exe)
            if w:
                eng[key] = w
    # 系统信息
    sysinfo = {"cpu": "", "cores": os.cpu_count() or 0, "ram_gb": 0.0,
               "disks": [], "gpu": ""}
    # v1.0：弃 PowerShell（客户端 AV 悬挂 PS 子进程的实弹前科）——纯 Python
    try:
        sysinfo["cpu"] = os.environ.get("PROCESSOR_IDENTIFIER", "")
    except Exception:
        pass
    try:
        import ctypes
        class _MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        ms = _MS(); ms.dwLength = ctypes.sizeof(_MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
        sysinfo["ram_gb"] = round(ms.ullTotalPhys / 1024**3, 1)
    except Exception:
        pass
    for drv in "CDEFG":
        p = drv + ":\\"
        if os.path.isdir(p):
            try:
                import shutil as _sh2
                t, used, free = _sh2.disk_usage(p)
                sysinfo["disks"].append("%s: 剩余 %.0f/%.0f GB" % (drv, free/1024**3, t/1024**3))
            except Exception:
                pass
    # v1.0：GPU 恢复显示——改用注册表 winreg 直读显示适配器类键（零子进程，
    # 与 v98 弃 PS 原则兼容：Display Adapter 类 {4d36e968-…} 下各 000N 子键的
    # DriverDesc 即显卡名）。失败留空不阻断。
    # v1.0: GPU restored via an in-process registry read (no subprocess,
    # keeping the v98 no-PS principle): DriverDesc under the display-adapter
    # class {4d36e968-...} subkeys 000N gives the adapter names. Empty on
    # failure, never blocking.
    try:
        import winreg
        _cls = r"SYSTEM\CurrentControlSet\Control\Class\{4d36e968-e325-11ce-bfc1-08002be10318}"
        _gpus = []
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, _cls) as _k:
            _i = 0
            while True:
                try:
                    _sub = winreg.EnumKey(_k, _i); _i += 1
                except OSError:
                    break
                if not _sub.isdigit():
                    continue
                try:
                    with winreg.OpenKey(_k, _sub) as _sk:
                        _desc, _ = winreg.QueryValueEx(_sk, "DriverDesc")
                    # 剔除远程/虚拟显示适配器（向日葵 Oray、ToDesk、Parsec、
                    # VMware、Basic Display）——实弹本机曾混入 OrayIddDriver
                    # drop remote/virtual display adapters (Oray, ToDesk,
                    # Parsec, VMware, Basic Display)
                    _d = _desc.lower() if _desc else ""
                    if (_desc and _desc not in _gpus and not any(
                            x in _d for x in ("oray", "todesk", "parsec",
                                              "vmware", "virtual", "basic display"))):
                        _gpus.append(_desc)
                except OSError:
                    continue
        sysinfo["gpu"] = " / ".join(_gpus[:3])
    except Exception:
        sysinfo["gpu"] = ""
    core4 = all(eng[k] for k in ("pfind_exe", "pdeep2_exe", "punifind_exe", "punifind_python"))
    verdict = ("核心依赖齐备（pFind/pDeep2/pUniFind/Python 环境）\ncore dependencies ready"
               if core4 else "部分核心依赖未检出，请按下方清单手动填写\nsome core dependencies missing - fill manually")
    return {"ok": True, "engines": eng, "system": sysinfo, "verdict": verdict}
# ---------------------------------------------------------------------------

_DETECT_CACHE = None
_DETECT_TS = 0.0


def _pick_path(kind):
    # v1.0：Windows 首选**进程内 ctypes** Win32 对话框（零子进程——PS
    # Forms 在实弹客户端被无声悬挂，请求永久挂起）；ctypes 非取消性失败才
    # 落 PowerShell 后备；两路都取消→"未选择"；都失败→带两路诊断。
    # v1.0: on Windows the IN-PROCESS ctypes dialogs come first (no child
    # process - PS Forms hung silently on the live client, leaving the
    # request pending forever); PowerShell is the fallback only when ctypes
    # fails for a non-cancel reason; both-cancel -> "cancelled";
    # both-fail -> diagnostics from both routes.
    with _TK_LOCK:
        if os.name == "nt":
            ok, r = _pick_path_ctypes(kind)
            if ok:
                _PICK_ROUTE[0] = "ctypes"
                return True, r
            if r == "cancelled":
                _PICK_ROUTE[0] = "ctypes"
                return False, "未选择\ncancelled"
            _PICK_ROUTE[0] = "ctypes-failed->ps"
            ok2, r2 = _pick_path_ps(kind)
            if ok2:
                _PICK_ROUTE[0] = "ps"
                return True, r2
            if r2 in ("non-windows", "cancelled"):
                return False, "未选择\ncancelled"
            return False, ("选择器执行失败，请在文本框手动输入路径\npicker "
                           "failed - type the path manually:\nctypes: %s\nps: %s"
                           % (r, r2))
        if os.name != "nt":
            try:
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk()
                root.withdraw()
                root.attributes("-topmost", True)
                if kind == "dir":
                    p = filedialog.askdirectory(title="选择文件夹 Select folder")
                elif kind == "files":
                    ps_ = filedialog.askopenfilenames(title="选择多个文件 Select files")
                    p = "\n".join(ps_) if ps_ else ""
                else:
                    p = filedialog.askopenfilename(title="选择文件 Select file")
                root.destroy()
                if p:
                    return True, p
                return False, "未选择\ncancelled"
            except Exception as e:  # noqa
                return False, "无法打开系统对话框，请直接在文本框手动输入路径\ndialog unavailable - type the path manually: %s" % e
        # Windows 上 PowerShell 也没弹出（用户取消或被策略禁用）
        # PowerShell returned nothing (user cancelled or blocked by policy)
        return False, "未选择\ncancelled"


def read_config_fields(path):
    """读取一个参数文件（config_step1.txt 格式：键<TAB>值＋# 注释），
    返回其中全部键值对。编码先试 UTF-8（网页壳子生成的），失败退 GBK
    （记事本在中文 Windows 上保存的）。
    Read a parameter file (config_step1.txt layout: key<TAB>value + # comments)
    and return all key/value pairs. Encoding: UTF-8 first (web-shell generated),
    then GBK (Notepad saves ANSI/GBK on Chinese Windows)."""
    try:
        text, _enc = read_config_text(path)
    except OSError as e:
        return False, "无法读取文件\ncannot read the file: %s" % e, {}
    if text is None:
        return False, "文件编码无法识别\nfile encoding not recognized", {}
    fields = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        kv = line.split("\t", 1)
        key = kv[0].strip()
        if not key:
            continue
        fields[key] = kv[1].strip() if len(kv) > 1 else ""
    return True, "", fields


def count_tsv_col(path, col):
    """流式统计一个 TSV 文件某列的取值计数（返回 (Counter, 总行数)；
    文件缺失或列不存在返回 (None, None)）。
    Stream-count the values of one TSV column (returns (Counter, total rows);
    (None, None) when the file or column is missing)."""
    if not os.path.exists(path):
        return None, None
    c = Counter()
    total = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if col not in hdr:
                return None, None
            j = hdr.index(col)
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if j < len(f):
                    c[f[j]] += 1
                total += 1
    except OSError:
        return None, None
    return c, total


def norm_state(v):
    """归一历史成对状态口径：Paired/paired/Matched/matched/Double/double 等历史
    产物取值统一为 Single/Paired（大小写与旧术语归一；未知值原样返回）。
    Normalize historical paired-state spellings (Paired/paired/Matched/matched/
    Double/double) to Single/Paired; unknown values pass through unchanged."""
    v = (v or "").strip().lower()
    return {"single": "Single", "paired": "Paired",
            "matched": "Paired", "double": "Paired"}.get(v, (v or "").strip())


def count_tsv_col_any(path, cols):
    """按列名优先级读第一个存在的列；状态类列的键经 norm_state 归一后合并计数
    （大小写/历史口径统一显示）。
    Read the first existing column among cols; state columns get their keys
    normalized and merged via norm_state (case/legacy spellings unified)."""
    for c in cols:
        dist, total = count_tsv_col(path, c)
        if dist is not None:
            if "state" in c:
                merged = {}
                for k, v in dist.items():
                    nk = norm_state(k)
                    merged[nk] = merged.get(nk, 0) + v
                return merged, total
            return dist, total
    return None, None


def _read_cfg(base):
    """读取结果文件夹 config_step1.txt 的键值（UTF-8 优先 GBK 兜底）。
    Read the key/value pairs of the result folder's config_step1.txt."""
    fields = {}
    p = os.path.join(base, "config_step1.txt")
    if not os.path.exists(p):
        return fields
    try:
        text, _enc = read_config_text(p)
    except OSError:
        return fields
    if not text:
        return fields
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        kv = ln.split("\t", 1)
        if len(kv) == 2:
            fields[kv[0].strip()] = kv[1].strip()
    return fields


_MGF_CACHE = {}     # mgf 路径 -> {TITLE: (块起始字节, 块结束字节)}
_MODMASS = None     # AA/修饰质量表缓存


def _config_mgf_files():
    """从 config_step1 的 all_mgf/mgf 键展开 MGF 文件清单（文件夹展开 *.mgf）。
    浏览模式下只采用浏览项目同层 unlabeled/（另一台电脑的结果，本机 config
    的 mgf 路径不适用，避免同标题旧 MGF 抢谱）。"""
    bc = _browse_ctx()
    if bc:
        files = []
        if bc.get("mgf_dir"):
            for f in sorted(os.listdir(bc["mgf_dir"])):
                if f.lower().endswith(".mgf"):
                    files.append(os.path.join(bc["mgf_dir"], f))
        seen, uniq = set(), []
        for f in files:
            k = os.path.normcase(os.path.realpath(f))
            if k not in seen:
                seen.add(k)
                uniq.append(f)
        return uniq
    fields = _read_cfg(CURRENT_RESULT or CLIENT_ROOT)
    files = []
    for key in ("all_mgf", "mgf"):
        for p in fields.get(key, "").replace("\r", "\n").split("\n"):
            p = p.strip()
            if not p:
                continue
            if os.path.isdir(p):
                for f in sorted(os.listdir(p)):
                    if f.lower().endswith(".mgf"):
                        files.append(os.path.join(p, f))
            elif os.path.exists(p):
                files.append(p)
    # 浏览模式：谱图来自浏览项目同层的 unlabeled/（另一台电脑的结果，本机
    # config 的 mgf 路径不存在）
    seen, uniq = set(), []
    for f in files:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def _index_mgf(path):
    """建立 MGF 的 TITLE 字节偏移索引（懒建缓存；标题->块字节区间）。
    Build a TITLE byte-offset index of an MGF (title -> block byte range),
    cached lazily on first use."""
    if path in _MGF_CACHE:
        return _MGF_CACHE[path]
    idx = {}
    try:
        with open(path, "rb") as fh:
            off = 0
            start = None
            title = None
            for raw in fh:
                n = len(raw)
                line = raw.decode("utf-8", "replace").strip()
                if line == "BEGIN IONS":
                    start = off
                    title = None
                elif line.startswith("TITLE=") and start is not None:
                    title = line[6:].strip()
                elif line == "END IONS" and start is not None and title is not None:
                    idx[title] = (start, off + n)
                    start = None
                off += n
    except OSError:
        idx = {}
    _MGF_CACHE[path] = idx
    return idx


def find_spectrum(file_name):
    """按 File_Name（.dta 标题）在配置 MGF 中定位谱图：精确标题优先，退而
    endswith/包含核心 token。返回 (peaks, title) 或 (None, None)。
    Locate a spectrum by File_Name (the .dta title) across the configured MGFs:
    exact title first, then endswith/contains fallback on the core token.
    Returns (peaks, title) or (None, None)."""
    core = file_name[:-4] if file_name.lower().endswith(".dta") else file_name
    for p in _config_mgf_files():
        idx = _index_mgf(p)
        hit = None
        used = None
        if file_name in idx:
            hit, used = idx[file_name], file_name
        else:
            for t, rng in idx.items():
                if t.endswith(file_name) or t.endswith(core) or core in t:
                    hit, used = rng, t
                    break
        if hit is None:
            continue
        peaks = []
        try:
            with open(p, "rb") as fh:
                fh.seek(hit[0])
                block = fh.read(hit[1] - hit[0]).decode("utf-8", "replace")
        except OSError:
            return None, None
        for line in block.splitlines():
            t = line.strip()
            if not t or "=" in t or t in ("BEGIN IONS", "END IONS"):
                continue
            parts = t.split()
            if len(parts) >= 2:
                try:
                    peaks.append({"mz": float(parts[0]), "int": float(parts[1])})
                except ValueError:
                    pass
        return peaks, used
    return None, None


def parse_candidates(s):
    """候选串拆解：'SEQ1(mod1)|SEQ2(mod2)' -> [[seq, mod], ...]
    （rfind 取最外层括号，兼容修饰内含 [] 的 SAP 记法）。
    Split a candidates string into [seq, mod] pairs (rfind takes the outermost
    parentheses, tolerant of [] inside SAP-style modifications)."""
    out = []
    for part in (s or "").split("|"):
        part = part.strip()
        if not part:
            continue
        i = part.rfind("(")
        if i > 0 and part.endswith(")"):
            out.append([part[:i].strip(), part[i + 1:-1]])
        else:
            out.append([part, ""])
    return out


def _label_of(base):
    """当前结果文件夹的 label_type（小写；缺省 unlabeled）。浏览模式下取浏览
    上下文的 label。"""
    bc = _browse_ctx()
    if bc:
        return bc["label"]
    return (_read_cfg(base).get("label_type", "unlabeled") or "unlabeled").lower()


def spectrum_info(file_name):
    """取一张谱的候选串与结论。候选全列表读 paired_report（conflict_report 的
    candidates 列仅含胜选单条）；结论与胜选肽段读 conflict_report（回退
    stageA report_A）。返回 {"cands": [[seq,mod]...], "winner_index": int,
    "conclusion": str}。Fetch one spectrum's candidates and conclusion: the
    full candidate list from paired_report (the conflict_report candidates
    column holds only the winner), the conclusion and winner peptide from
    conflict_report (stageA report_A as fallback)."""
    pbase = CURRENT_RESULT or CLIENT_ROOT
    base = _result_base(pbase)
    label = _label_of(pbase)
    cands = None
    pr_path = os.path.join(base, "2.labeled_result", "paired_report_%s.tsv" % label)
    if os.path.exists(pr_path):
        try:
            with open(pr_path, "r", encoding="utf-8", errors="replace") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                jf, jc = hdr.index("File_Name"), hdr.index("candidates")
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) > jf and f[jf] == file_name:
                        cands = parse_candidates(f[jc] if jc < len(f) else "")
                        break
        except OSError:
            pass
    if cands is None:
        cands = []
    for rel in ("2.labeled_result/conflict_report_%s.tsv" % label, "1.unlabeled_result/report_A.tsv"):
        p = os.path.join(base, rel)
        if not os.path.exists(p):
            continue
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                if "File_Name" not in hdr or "candidates" not in hdr:
                    continue
                jf = hdr.index("File_Name")
                jc = hdr.index("candidates")
                jw = hdr.index("winner_sequence") if "winner_sequence" in hdr else -1
                jwm = hdr.index("winner_modification") if "winner_modification" in hdr else -1
                jcs = hdr.index("conclusion_status") if "conclusion_status" in hdr else -1
                jde = hdr.index("decision_A") if "decision_A" in hdr else -1
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) > jf and f[jf] == file_name:
                        wseq = f[jw] if 0 <= jw < len(f) else ""
                        wmod = f[jwm] if 0 <= jwm < len(f) else ""
                        wi = -1
                        for i, c in enumerate(cands):
                            if c[0] == wseq and (not wmod or not c[1] or c[1] == wmod):
                                wi = i
                                break
                        concl = ""
                        if 0 <= jcs < len(f) and f[jcs]:
                            concl = f[jcs]
                        elif 0 <= jde < len(f):
                            concl = f[jde]
                        return {"cands": cands, "winner_index": wi,
                                "conclusion": concl}
        except OSError:
            continue
    return None


def spectra_list(limit=2000):
    """可浏览谱清单＝最后确定的"原始结果错误"谱图（用户口径 2026-09-11）：
    冲突报告 conclusion_status == "original result wrong"（标记样本）；
    无冲突报告时回退 report_A 的 decision_A == "candidate_misidentification"。
    The browsable list = the FINALIZED original-result-wrong spectra (user
    ruling 2026-09-11): conflict_report conclusion_status == "original result
    wrong" (labeled samples); fallback report_A decision_A ==
    "candidate_misidentification"."""
    pbase = CURRENT_RESULT or CLIENT_ROOT
    base = _result_base(pbase)
    label = _label_of(pbase)
    for rel, flag_col, flag_val in (
            ("2.labeled_result/conflict_report_%s.tsv" % label,
             "conclusion_status", "original result wrong"),
            ("1.unlabeled_result/report_A.tsv",
             "decision_A", "candidate_misidentification")):
        p = os.path.join(base, rel)
        if not os.path.exists(p):
            continue
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                if "File_Name" not in hdr:
                    continue
                jf = hdr.index("File_Name")
                jf_flag = hdr.index(flag_col) if flag_col in hdr else -1
                out = []
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if not (len(f) > jf and f[jf]):
                        continue
                    if 0 <= jf_flag < len(f) and f[jf_flag] != flag_val:
                        continue
                    out.append([f[jf], f[jf_flag] if 0 <= jf_flag < len(f) else ""])
                    if len(out) >= limit:
                        break
                if out:
                    return True, out
        except OSError:
            continue
    return False, []


def load_modmass():
    """解析 code/AA_mass_dataframe.txt -> {"aa": {字母: 单同位素质量},
    "mods": {修饰名: 质量增量}}（懒缓存；单字母行入 aa，其余入 mods）。
    Parse the AA/mod mass table into {"aa": {letter: mono mass}, "mods": {name:
    delta}} (cached; single-letter rows go to aa, the rest to mods)."""
    global _MODMASS
    if _MODMASS is not None:
        return _MODMASS
    aa, mods = {}, {}
    p = os.path.join(CODE_DIR, "AA_mass_dataframe.txt")
    try:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ja = hdr.index("AA") if "AA" in hdr else 0
            jm = hdr.index("Mono_weight") if "Mono_weight" in hdr else 7
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) <= ja or not f[ja]:
                    continue
                m = None
                if len(f) > jm:
                    try:
                        m = float(f[jm])
                    except ValueError:
                        m = None
                if m is None:
                    for cell in f[ja + 1:]:
                        try:
                            m = float(cell)
                            break
                        except ValueError:
                            continue
                if m is None:
                    continue
                k = f[ja]
                if len(k) == 1 and k.isalpha():
                    aa[k] = m
                else:
                    mods[k] = m
    except OSError:
        pass
    _MODMASS = {"aa": aa, "mods": mods}
    return _MODMASS


def find_pdfs(dirs, depth=3, cap=500):
    """在给定目录（限深度）内扫描 PDF 文件，返回 [{name, path, size_mb}]。
    Scan for PDF files under the given directories (depth-limited)."""
    out, seen = [], set()
    for d in dirs:
        if not d or not os.path.isdir(d):
            continue
        base_depth = d.rstrip(os.sep).count(os.sep)
        for root, sub, files in os.walk(d):
            if root.count(os.sep) - base_depth >= depth:
                sub[:] = []
            for f in sorted(files):
                if f.lower().endswith(".pdf"):
                    p = os.path.join(root, f)
                    if p not in seen:
                        seen.add(p)
                        try:
                            sz = os.path.getsize(p) / 1e6
                        except OSError:
                            continue
                        out.append({"name": f, "path": p,
                                    "size_mb": round(sz, 1)})
                    if len(out) >= cap:
                        return out
    return out


def pdf_scan_dirs(extra=None):
    """默认扫描根：当前结果文件夹（用户口径 2026-09-11：只展示当前运行的
    PDF，不再扫软件根目录以免归档旧 PDF 混入）；无当前结果时退回客户端根；
    另加用户手填的文件夹。
    Default scan roots: the CURRENT result folder only (user ruling 2026-09-11:
    show only the current run's PDFs; the client-root scan is removed so
    archived PDFs no longer leak in); falls back to the client root when no
    result exists; plus an optional user-supplied folder."""
    dirs = []
    if extra and os.path.isdir(extra):
        dirs.append(extra)
    if BROWSE_DIR:
        # 浏览模式只展示浏览目录的 PDF（本机结果不混入）
        dirs.append(_result_base(BROWSE_DIR))
    elif CURRENT_RESULT:
        dirs.append(CURRENT_RESULT)
    if not dirs:
        dirs.append(CLIENT_ROOT)
    return dirs


def count_conflicts(path):
    """流式统计 report_A.tsv 的总行数与冲突谱行数（Peptide_state 或
    Modification_state != Same）。
    Stream-count the total rows and conflict rows of report_A.tsv (Peptide_state
    or Modification_state != Same)."""
    if not os.path.exists(path):
        return None, None
    n_total = 0
    n_conf = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            jp = hdr.index("Peptide_state") if "Peptide_state" in hdr else -1
            jm = hdr.index("Modification_state") if "Modification_state" in hdr else -1
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                n_total += 1
                if jp >= 0 and jp < len(f) and f[jp] != "Same":
                    n_conf += 1
                elif jm >= 0 and jm < len(f) and f[jm] != "Same":
                    n_conf += 1
    except OSError:
        return None, None
    return n_total, n_conf


def count_error_reasons(path):
    """流式统计错误原因（仅 vote_says_wrong=1 且原因非空的行）。
    Stream-count the error reasons (only rows with vote_says_wrong=1 and a
    non-empty reason)."""
    if not os.path.exists(path):
        return Counter()
    c = Counter()
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            jw = hdr.index("vote_says_wrong") if "vote_says_wrong" in hdr else -1
            je = hdr.index("error_reason") if "error_reason" in hdr else -1
            if jw < 0 or je < 0:
                return Counter()
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if jw < len(f) and f[jw] == "1" and je < len(f) and f[je]:
                    c[f[je]] += 1
    except OSError:
        return Counter()
    return c


def read_tsv_cols(path, cols):
    """读取 TSV 的指定几列（返回行元组列表；文件或任一列缺失返回 None）。
    Read the given columns of a TSV (returns a list of row tuples; None when
    the file or any column is missing)."""
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if any(c not in hdr for c in cols):
                return None
            idx = [hdr.index(c) for c in cols]
            out = []
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                out.append(tuple(f[i] if i < len(f) else "" for i in idx))
            return out
    except OSError:
        return None


def _rows(d):
    out = []
    for k, v in d.items():
        if v:
            # 键/值均转义：error_reason 等可能携带引擎输出片段
            # escape both: error_reason may carry third-party engine output
            out.append("<tr><td>%s</td><td style='text-align:right'>%s</td></tr>"
                       % (html.escape(str(k)), html.escape(str(v))))
    return "".join(out) or "<tr><td colspan='2'>—</td></tr>"


def build_report_html():
    """按当前结果文件夹生成双语 HTML 结果报告（运行时间线＋阶段 A/B 统计＋
    文件清单）。所有数字直接流式读自结果文件夹。
    Build the bilingual HTML results report from the current result folder
    (run timeline + stage A/B statistics + output file list); all numbers are
    stream-read directly from the result folder."""
    base = CURRENT_RESULT or CLIENT_ROOT
    bc = _browse_ctx()
    if bc:
        base = bc["pbase"]
    rbase = _result_base(base)         # 结果文件在 result_时间戳 子文件夹内
    cfg = {}
    cfg_path = os.path.join(base, "config_step1.txt")
    if os.path.exists(cfg_path):
        try:
            text, _enc = read_config_text(cfg_path)
        except OSError:
            text = None
        if text:
            for ln in text.splitlines():
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                kv = ln.split("\t", 1)
                if len(kv) == 2:
                    cfg[kv[0].strip()] = kv[1].strip()
    label = (bc["label"] if bc else (cfg.get("label_type", "unlabeled") or "unlabeled").lower())

    tl = None
    tl_hint = ""
    if bc:
        # 自查整改（2026-09-26）：浏览模式下不展示本机时间线；浏览项目内有
        # web_run.log 时解析之，否则明确提示
        log_cand = next((c for c in (os.path.join(bc["pbase"], "web_run.log"),
                                     os.path.join(bc["rbase"], "web_run.log"))
                         if os.path.isfile(c)), None)
        if log_cand:
            try:
                with open(log_cand, "r", encoding="utf-8", errors="replace") as fh:
                    tl = parse_timeline(fh.read())
            except OSError:
                tl = None
        if not tl:
            tl_hint = "浏览目录内无运行日志（本机时间线不展示）"
    else:
        tl = current_timeline() or []
    if tl_hint and not tl:
        tl_rows = "<tr><td colspan='4' class='hint'>%s</td></tr>" % esc(tl_hint)
    else:
        tl_rows = "".join(
        "<tr><td><b>%s</b><br><span class='en'>%s</span></td>"
        "<td>%s</td><td>%s</td><td>%s</td></tr>"
        % (esc(m.get("name_zh", "")), esc(m.get("name_en", "")),
           esc(m.get("start", "—")), esc(m.get("end", "—") or "进行中 running"),
           esc({"done": "完成 done", "stop": "停点 stop", "failed": "失败 failed",
                "running": "进行中 running"}.get(m.get("status", ""), "—")))
        for m in tl)

    # 阶段 A
    rep_a = os.path.join(rbase, "1.unlabeled_result", "report_A.tsv")
    n_total, n_conf = count_conflicts(rep_a)
    dec, _ = count_tsv_col(rep_a, "decision_A")
    cla = os.path.join(rbase, "1.unlabeled_result", "classification_A.tsv")
    cls, _ = count_tsv_col(cla, "conclusion_status")
    cls_wrong = count_error_reasons(cla)
    sum_txt = ""
    sum_a = os.path.join(rbase, "1.unlabeled_result", "summary_A.txt")
    if os.path.exists(sum_a):
        try:
            with open(sum_a, "r", encoding="utf-8", errors="replace") as fh:
                lines = fh.read().splitlines()
                if lines and "\t" in lines[0]:
                    lines = lines[1:]          # 去掉表头行
                sum_txt = "<br>".join(esc(l.replace("\t", " : ")) for l in lines)
        except OSError:
            sum_txt = ""

    # 阶段 B（标记样本）
    stb_html = "<p class='hint'>无标记样本——无阶段 B 输出。unlabeled sample - no stage B outputs.</p>"
    if label != "unlabeled":
        pr = os.path.join(rbase, "2.labeled_result", "paired_report_%s.tsv" % label)
        ds, n_pr = count_tsv_col_any(pr, ["paired_state", "double_state"])
        pm, _ = count_tsv_col(pr, "pair_missing")
        cr = os.path.join(rbase, "2.labeled_result", "conflict_report_%s.tsv" % label)
        cs, n_cr = count_tsv_col(cr, "conclusion_status")
        # 分组与肽段层面统计：Single / Paired 分开 ＋ 肽段三 panel
        # grouped and peptide-level stats: Single vs Paired + three panels
        # 列名兼容：状态列 paired_state_B（旧版 double_state_B）；肽段序列列
        # a.Sequence（多引擎版）或 X1.Sequence（经典版）；error_reason/
        # vote_says_wrong 可选——先读表头一次性探测可用列，缺列自动降级
        # column-name tolerance: the state column is paired_state_B
        # (double_state_B in legacy); the peptide-sequence column is a.Sequence
        # (multi-engine) or X1.Sequence (classic); error_reason/vote_says_wrong
        # are optional - probe the header once and degrade gracefully
        gr = None
        try:
            with open(cr, "r", encoding="utf-8", errors="replace") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
        except OSError:
            hdr = []
        if hdr:
            st_col = "paired_state_B" if "paired_state_B" in hdr else (
                "double_state_B" if "double_state_B" in hdr else "")
            seq_col = "a.Sequence" if "a.Sequence" in hdr else (
                "X1.Sequence" if "X1.Sequence" in hdr else "")
            vote_col = "vote_says_wrong" if "vote_says_wrong" in hdr else ""
            err_col = "error_reason" if "error_reason" in hdr else ""
            need = [st_col, "conclusion_status", "winner_sequence"] + \
                   ([vote_col] if vote_col else []) + \
                   ([err_col] if err_col else []) + \
                   ([seq_col] if seq_col else [])
            if all(c in hdr for c in need) and st_col and seq_col:
                gr = read_tsv_cols(cr, need)
                gr = [(g[0], g[1], g[need.index(err_col)] if err_col else "",
                       g[need.index(vote_col)] if vote_col else "",
                       g[need.index("winner_sequence")],
                       g[need.index(seq_col)]) for g in gr]
        con_g = {"Single": Counter(), "Paired": Counter()}
        er_g = {"Single": Counter(), "Paired": Counter()}
        pep_g = {"Single": defaultdict(set), "Paired": defaultdict(set),
                 "Combined": defaultdict(set)}
        n_npc = 0
        if gr:
            for st_raw, cst, ers, wr, wseq, x1 in gr:
                st = norm_state(st_raw)
                pep = wseq or x1
                if st in con_g:
                    con_g[st][cst] += 1
                    if wr == "1" and ers:
                        er_g[st][ers] += 1
                    pep_g[st][cst].add(pep)
                else:
                    n_npc += 1
                pep_g["Combined"][cst].add(pep)

        def _gtable(cnt_s, cnt_p):
            keys = sorted(set(cnt_s) | set(cnt_p))
            rows = []
            for k in keys:
                a, b = cnt_s.get(k, 0), cnt_p.get(k, 0)
                if a or b:
                    rows.append("<tr><td>%s</td><td style='text-align:right'>%d</td>"
                                "<td style='text-align:right'>%d</td>"
                                "<td style='text-align:right'>%d</td></tr>"
                                % (esc(k), a, b, a + b))
            if not rows:
                return "<tr><td colspan='4'>—</td></tr>"
            tot = "<tr><th>合计 total</th><th style='text-align:right'>%d</th>" \
                  "<th style='text-align:right'>%d</th><th style='text-align:right'>%d</th></tr>" \
                  % (sum(cnt_s.values()), sum(cnt_p.values()),
                     sum(cnt_s.values()) + sum(cnt_p.values()))
            return "".join(rows) + tot

        def _pep_panel(name, sets):
            total = len(set().union(*sets.values())) if sets else 0
            rows = "".join("<tr><td>%s</td><td style='text-align:right'>%d</td></tr>"
                           % (esc(k), len(v))
                           for k, v in sorted(sets.items(), key=lambda kv: -len(kv[1])))
            return ("<p><b>%s</b>：去重肽段 unique peptides <b>%d</b></p><table>%s</table>"
                    % (esc(name), total, rows or "<tr><td colspan='2'>—</td></tr>"))

        stb_html = ("<p>成对匹配报告 paired report：<b>%s</b> 行 rows；pair_missing 标记 %s 行。</p>"
                    % (n_pr if n_pr is not None else "—",
                       pm.get("1", 0) if pm else "—") +
                    "<h4>paired_state 分布 distribution</h4><table>" + _rows(dict(ds or {})) +
                    "</table>" +
                    "<h4>冲突谱结论状态（按 Single / Paired 分组）conflict conclusion status by group</h4>"
                    "<table><tr><th>结论 conclusion</th><th>Single</th><th>Paired</th><th>合计 total</th></tr>"
                    + _gtable(con_g["Single"], con_g["Paired"]) + "</table>"
                    + ("<p class='hint'>not pair-checked（未核查）：%d 行 rows。</p>" % n_npc if n_npc else "") +
                    "<h4>错误原因（原始结果错误判定行，按组）error reasons by group</h4>"
                    "<table><tr><th>原因 reason</th><th>Single</th><th>Paired</th><th>合计 total</th></tr>"
                    + _gtable(er_g["Single"], er_g["Paired"]) + "</table>" +
                    "<h4>肽段层面统计（三个 panel）peptide-level stats (three panels)</h4>"
                    + _pep_panel("Single 谱贡献肽段 peptides from Single spectra", pep_g["Single"])
                    + _pep_panel("Paired 谱贡献肽段 peptides from Paired spectra", pep_g["Paired"])
                    + _pep_panel("合并（全部冲突谱）combined", pep_g["Combined"]))

    # 输出文件清单（config 在项目根；结果在 result_时间戳 子文件夹内）
    # output file list (config at the project root; results inside the
    # result_<timestamp> subfolder)
    files = []
    file_roots = [("config_step1.txt", base)] + \
                 [("1.unlabeled_result/report_A.tsv", rbase),
                  ("1.unlabeled_result/classification_A.tsv", rbase),
                  ("1.unlabeled_result/summary_A.txt", rbase),
                  ("merge/peptide.txt", rbase), ("merge/predict.txt", rbase)]
    if label != "unlabeled":
        file_roots += [("2.labeled_result/paired_report_%s.tsv" % label, rbase),
                       ("2.labeled_result/conflict_report_%s.tsv" % label, rbase)]
    for rel, root in file_roots:
        p = os.path.join(root, rel)
        if os.path.exists(p):
            files.append(rel + "  (" + ("%.1f MB" % (os.path.getsize(p) / 1e6)) + ")")
    files_html = "".join("<li>%s</li>" % esc(f) for f in files) or "<li>—</li>"

    # 两分法结果呈现（用户口径 2026-09-12）：软件固定输出两个部分——
    # ① 已注释库错误鉴定报告（原始鉴定蛋白在用户 FASTA 内，判校正＝错误匹配清单）
    # ② 可信新肽段报告（原始鉴定蛋白不在用户 FASTA 内，credible_novel＝竞争保留＋非混合谱）
    # two-part result presentation (user ruling 2026-09-12): the software always
    # reports both parts - (1) annotated-library misidentification report (the
    # original protein is in the user FASTA; judged-corrected rows are the wrong
    # matches), (2) credible novel-peptide report (original protein not in the
    # user FASTA; credible_novel = competition-kept and non-chimeric)
    twopart_rows = read_two_part_summary(rbase)
    if twopart_rows:
        tpb = []
        for sub, part, rows, head in twopart_rows:
            fp = os.path.join(rbase, sub, part)
            tpb.append(
                "<div class='tp-block'><h4>第 %s 部分 — %s <span class='en'>Part %s — %s</span></h4>"
                "<p class='hint'>%s</p>"
                "<ul><li>%s  <span class='sz'>(%.1f MB)</span></li></ul></div>" % (
                    esc(part[0]), esc(twopart_title(part, "zh")),
                    esc(part[0]), esc(twopart_title(part, "en")),
                    esc(head), esc(part),
                    os.path.getsize(fp) / 1e6
                    if os.path.exists(fp) else 0.0))
        twopart_html = (
            "<div class='card'><h2>两分法结果 <span class='en'>two-part results"
            "（① 已注释错误鉴定 ② 可信新肽段）</span></h2>"
            "<p class='hint'>软件固定同时输出两个部分，无需指定用途："
            "① 用户 FASTA 内肽段的错误鉴定报告（判校正＝错误匹配清单）；"
            "② 用户 FASTA 外肽段的可信新肽段报告（可信＝竞争后保留且非混合谱）。</p>"
            + "".join(tpb) + "</div>")
    else:
        twopart_html = ""

    # 谱图 PDF 卡片（确信原始结果错误一本，50 页口径；spectrum.js 驱动）
    # the spectrum-PDF card (confirmed-original-wrong volume, 50-page ruling)
    pdfon = read_pdf_on_flag(rbase)
    pdf_status = ("<span class='st on'>PDF_ON 已开启 enabled</span>" if pdfon == "1"
                  else "<span class='st off'>PDF_ON 默认关闭 disabled (default)</span>")
    pdfcard_html = (
        "<div class='card'><h2>谱图 PDF <span class='en'>spectrum PDFs"
        "（确信原始结果错误 confirmed-original-wrong，50 页口径）</span></h2>"
        "<p class='hint'>" + pdf_status +
        "；自动扫描结果文件夹与软件根目录；也可输入任一文件夹路径后点查找。"
        "点击文件名在窗口内查看，iframe 顶部为浏览器自带 PDF 工具条（含打印）。</p>"
        "<input id='ispsv_pdf_dir' size='52' placeholder='文件夹路径（可选）'>"
        " <button id='ispsv_pdf_scan'>查找 scan</button>"
        "<div id='ispsv_pdf_list'></div>"
        "<iframe id='ispsv_pdf_frame'"
        " style='width:100%;height:640px;border:1px solid #dce5ee;background:#fff'></iframe>"
        "</div>")

    # 谱图浏览卡片（spectrum.js 交互呈现；数据走 /api/spectrum 等端点）
    # spectrum-viewer card (interactive via spectrum.js; data from the APIs)
    spec_html = (
        "<div class='card'><h2>谱图浏览 <span class='en'>spectrum viewer</span></h2>"
        "<p class='hint'>双面板：上＝候选 1（原始），下＝下拉选择的竞争/校正候选。"
        "理论梯 a/b/c/x/y/z（1+/2+）＋M 系，20 ppm 窗最近峰匹配；标注规则与 PDF 一致"
        "（电荷号不省略、同峰合并 by 为准 y 优先、相邻峰窜开）。</p>"
        "<input id='ispsv_file' list='ispsv_datalist' size='52'"
        " placeholder='File_Name（.dta 标题，可输入筛选）'>"
        "<datalist id='ispsv_datalist'></datalist>"
        " <button id='ispsv_load'>加载 load</button>"
        " <select id='ispsv_cand_sel'></select>"
        "<p class='hint' id='ispsv_msg'></p>"
        "<canvas id='ispsv_canvas' width='900' height='560'"
        " style='max-width:100%;background:#fff;border:1px solid #dce5ee'></canvas>"
        "</div><script src='/spectrum.js'></script>")

    pct_str = ("%.1f%%" % (100.0 * n_conf / n_total)) if n_total and n_conf is not None else "—"

    # 判校正行的成对标记证据（correction_paired_summary.txt）
    # paired-label evidence among corrected rows
    paired_card = ""
    if label != "unlabeled":
        cps = read_correction_paired_summary(rbase)
        if cps:
            lines = cps.strip().splitlines()
            # 结构：前两行 = 汇总键/值（如 corrections_total=4839）；
            # 第三行起 = 五类分列表（category n paired_labeled paired_cos_gt_0.7）
            # layout: lines 0-1 are the summary key/value row
            # (e.g. corrections_total=4839); line 2 onwards is the per-category
            # table (category / n / paired_labeled / paired_cos_gt_0.7)
            html_bits = []
            if len(lines) >= 2:
                keys = lines[0].split("\t")
                vals = lines[1].split("\t")
                kvs = "".join(
                    "<span class='kpi'>%s <b>%s</b></span>" % (esc(k.replace("_", " ")), esc(vals[i]) if i < len(vals) else "—")
                    for i, k in enumerate(keys))
                html_bits.append("<p>%s</p>" % kvs)
            if len(lines) >= 4:
                rows_html = []
                for ln in lines[3:]:
                    cells = ln.split("\t")
                    rows_html.append("<tr>" + "".join(
                        "<td>%s</td>" % esc(c) for c in cells) + "</tr>")
                hdr_cells = "".join("<th>%s</th>" % esc(c) for c in lines[2].split("\t"))
                html_bits.append("<table><tr>%s</tr>%s</table>"
                                 % (hdr_cells, "".join(rows_html)))
            paired_card = ("<div class='card hl-g'><h2>判校正行的成对标记证据"
                           " <span class='en'>paired-label evidence of corrected rows</span></h2>"
                           + "".join(html_bits) + "</div>")
    # 五类合并错误分布（主表 diff_type → MERGE5 映射）
    # five merged error categories
    merge5_card = ""
    m5 = count_merge5_errors(rbase)
    if m5 and m5["total"]:
        rows = "".join(
            "<tr><td>%s</td><td style='text-align:right'><b>%d</b></td>"
            "<td style='width:42%%;min-width:180px'>%.1f%%"
            "<div class='bar'><div style='width:%.1f%%'></div></div></td></tr>"
            % (esc(k), v, 100.0 * v / m5["total"], 100.0 * v / m5["total"])
            for k, v in m5["categories"].items())
        merge5_card = ("<div class='card hl'><h2>错误类型分布（五类合并）"
                       " <span class='en'>error-type distribution (5 merged categories)</span></h2>"
                       "<p>判校正 <b>" + str(m5["total"]) + "</b> 行</p>"
                       "<table><tr><th>类型 type</th><th style='text-align:right'>谱图 spectra</th>"
                       "<th>占比 pct</th></tr>" + rows + "</table></div>")
    # 策略归因（哪个引擎贡献了被替换的候选）
    # strategy attribution of corrected candidates
    attribution_card = ""
    attr = count_strategy_attribution(rbase)
    if attr:
        rows = "".join("<tr><td>%s</td><td><b>%d</b></td></tr>" % (esc(k), v)
                       for k, v in attr.items())
        attribution_card = ("<div class='card hl-o'><h2>校正来源策略归因"
                            " <span class='en'>strategy attribution of corrections</span></h2>"
                            "<table><tr><th>策略槽 slot</th><th>谱图 spectra</th></tr>%s</table>"
                            "<p class='hint'>被替换候选的序列精确匹配合并表策略列；"
                            "未归属行多为仲裁自身 SAP 机制产生的同序列蛋白归属调整。</p></div>" % rows)
    html = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8"><title>isoPSV 结果报告</title>
<style>
 :root { --ac:#0b6b8a; --ac2:#0e8a5f; --ink:#1e2d3d; --sub:#5b7284; --line:#d5dfe8; --soft:#eaf3f8; }
 * { box-sizing:border-box; }
 body { font-family:"Segoe UI","Microsoft YaHei","PingFang SC",sans-serif; margin:0; padding:22px 28px 60px; color:var(--ink); background:#f4f7fa; line-height:1.6; }
 .en { color:#8ca3b5; font-weight:normal; font-size:.88em; }
 h1 { font-size:1.35em; margin:0 0 2px; color:var(--ac); display:flex; align-items:center; gap:8px; }
 h1::before { content:"◈"; font-size:1.05em; }
 h2 { font-size:1.08em; margin:0 0 10px; color:var(--ac); display:flex; align-items:center; gap:8px;
      border-left:4px solid var(--ac); padding-left:10px; }
 h4 { margin:12px 0 4px; color:var(--sub); font-size:.95em; }
 .card { background:#fff; border:1px solid var(--line); border-radius:10px; padding:16px 20px;
         margin:0 0 14px; box-shadow:0 1px 3px rgba(20,50,80,.07); }
 table { border-collapse:collapse; margin:6px 0; font-size:.9em; width:100%%; }
 td, th { border:1px solid #dde5ee; padding:5px 10px; text-align:left; }
 th { background:var(--soft); color:var(--ac); font-weight:600; }
 tr:nth-child(even) td { background:#f8fafc; }
 tr:hover td { background:#f0f6fa; }
 td b { color:var(--ac); }
 .hint { color:#8ca3b5; font-size:.88em; }
 li { margin:2px 0; }
 .kpi { display:inline-flex; align-items:baseline; gap:6px; background:var(--soft); color:var(--ac);
        border-radius:8px; padding:3px 12px; margin:0 8px 4px 0; font-size:.92em; }
 .kpi b { font-size:1.12em; }
 .bar { height:7px; background:#e3ecf2; border-radius:4px; overflow:hidden; margin-top:3px; min-width:80px; }
 .bar > div { height:100%%; background:linear-gradient(90deg,#4a90c9,var(--ac)); border-radius:4px; }
 .card.hl { border-left:4px solid var(--ac); }
 .card.hl-g { border-left:4px solid var(--ac2); }
 .card.hl-o { border-left:4px solid #c9973f; }
 .small { font-size:.85em; color:var(--sub); }
</style></head><body>
<h1>isoPSV 结果报告 <span class="en">results report</span></h1>
<p class="hint">结果文件夹 result folder：%s<br>生成时间 generated：%s</p>
<div class="card"><h2>运行时间线 <span class="en">run timeline</span></h2>
<table><tr><th>模块 module</th><th>开始 start</th><th>结束 end</th><th>状态 status</th></tr>%s</table>
</div>
<div class="card"><h2>阶段 A <span class="en">stage A</span></h2>
<p>无标记谱图 unlabeled spectra：<b>%s</b>　冲突谱 conflict spectra：<b>%s</b>（%s）</p>
<p>判定 decision：%s</p>
<p>%s</p>
<h4>逐谱结论状态 conclusion status（全部谱图 all spectra）</h4>
<table>%s</table>
<h4>错误原因（原始结果错误判定行）error reasons</h4>
<table>%s</table>
</div>
<div class="card"><h2>阶段 B <span class="en">stage B</span>（label_type = %s）</h2>
%s</div>
<div class="card"><h2>谱图浏览 <span class="en">spectrum viewer</span></h2>
%s</div>
<div class="card"><h2>谱图 PDF <span class="en">spectrum PDFs</span></h2>
%s</div>
%s%s%s%s
<div class="card"><h2>输出文件 <span class="en">output files</span></h2>
<ul>%s</ul></div>
</body></html>""" % (
        esc(base), time.strftime("%Y-%m-%d %H:%M:%S"), tl_rows,
        esc(n_total if n_total is not None else "—"),
        esc(n_conf if n_conf is not None else "—"),
        pct_str,
        esc("　".join("%s: %d" % kv for kv in sorted((dec or {}).items()))),
        sum_txt,
        _rows(dict(cls or {})), _rows(dict(cls_wrong or {})),
        esc(label), stb_html, spec_html, pdfcard_html, twopart_html,
        paired_card, merge5_card, attribution_card, files_html)
    return html


def _sanitize_inline_field(v):
    """剥离换行/制表符/NUL（P0-8 通用化，2026-09-26 MiniMax P0-2；N-20 扩展，
    2026-09-26 GLM）：单值字段含换行会被拼接成新的键值行（任意键注入）；含 \t
    会在 R 端按制表符解析出歧义键，含 \0 会让 R readLines 报错——一律剔除。
    多行值字段（raw 等）不在此处理。"""
    return str(v or "").replace("\r", "").replace("\n", "").replace("\t", "").replace("\x00", "").strip()


def build_step2_config(strategies, msfragger_xmx=""):
    """从表单策略行生成 config_step2.txt 文本（search_file + search_file_type 逐行；
    MaxQuant 行附 search_file_raw 记录 raw 输入——MaxQuant 不能解析 MGF，
    用户需先对 raw 跑 MaxQuant，再把 msms.txt 填为结果文件）。具备自动搜库能力
    的引擎行（msfragger/comet/sage/metamorpheus/maxquant）可附 exe 软件路径——
    写入 <engine>_exe 与 <engine>_mode=1，管线按内置模板回退参数并自动运行
    （2026-09-24 网页策略行一键填充路径配套）；msfragger_xmx 非空时写入
    全局键（N-10，2026-09-25）。
    Build the config_step2.txt text from the form strategy rows (search_file +
    search_file_type one pair per line; MaxQuant rows also record search_file_raw
    - MaxQuant cannot parse MGF, so the user runs it on the raws first and the
    msms.txt goes into the result-file box). Rows of auto-capable engines
    (msfragger/comet/sage/metamorpheus/maxquant) may carry an exe software path -
    it is written as <engine>_exe with <engine>_mode=1 and the pipeline falls
    back to the bundled parameter templates and runs the engine automatically
    (companion to the 2026-09-24 web one-click path fill)."""
    lines = ["# config_step2.txt — 多引擎策略配置（网页版生成）",
             "# config_step2.txt - multi-engine strategy configuration (generated by the web shell)"]
    auto_capable = {"msfragger", "comet", "sage", "metamorpheus", "maxquant"}
    wrote_msconvert = False
    for s in strategies:
        sf = _sanitize_inline_field(s.get("file"))
        st = _sanitize_inline_field(s.get("type"))
        exe = _sanitize_inline_field(s.get("exe"))
        # v1.0：行内显式 mode（auto/manual）；缺省时沿用旧推断
        # （file 空＋exe 非空＝auto）以兼容旧表单。manual 行忽略 exe。
        # explicit per-row mode (v1.0); rows without mode keep the legacy
        # inference (empty file + exe = auto). exe is ignored on manual rows.
        mode = _sanitize_inline_field(s.get("mode") or "")
        is_auto = ((mode == "auto") if mode else (not sf and bool(exe))) and st in auto_capable
        if not sf and not is_auto:
            continue
        # 杨翰反馈（2026-10-08）：auto 模式下 search_file 为空但 search_file_type
        # 仍写入 → 管线看到 type 无 file 挂在"等待手动输入"。写入占位标记，
        # 管线检测到后走适配器自动搜库而非停等。
        # tester feedback: auto mode had empty search_file causing the pipeline
        # to hang at "awaiting input"; write a marker so the adapter runs.
        if is_auto:
            if not exe:
                continue
            if not sf:
                sf = "(auto:" + st + ")"
        if sf:
            lines.append("search_file\t" + sf)
        lines.append("search_file_type\t" + st)
        if st == "maxquant":
            # raw 输入每行一个文件；写入单行（";" 连接）以保持与 search_file 行序一一对应
            # one raw per line in the form; joined with ";" into ONE config line so
            # the parallel key stays row-aligned with search_file
            raw = (s.get("raw") or "").strip()
            if raw:
                raw1 = ";".join(x.strip() for x in raw.replace("\r", "").split("\n") if x.strip())
                lines.append("search_file_raw\t" + raw1)
        if is_auto:
            lines.append(st + "_exe\t" + exe)
            lines.append(st + "_mode\t1")
            mc = _sanitize_inline_field(s.get("msconvert"))
            if mc and not wrote_msconvert:
                lines.append("msconvert_exe\t" + mc)
                wrote_msconvert = True
    if msfragger_xmx:
        lines.append("msfragger_xmx\t" + _sanitize_inline_field(msfragger_xmx))
    return "\n".join(lines) + "\n"


def get_active_strategies():
    """读取当前结果文件夹的合并表策略槽列表（列名前缀 a–z）。
    Read the strategy-slot letters from the current result folder's merge table."""
    base = CURRENT_RESULT or CLIENT_ROOT
    rbase = _result_base(base)
    merge = os.path.join(rbase, "merge", "1.combine_result.spectra")
    if not os.path.exists(merge):
        return []
    slots = []
    try:
        with open(merge, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
            for col in header.rstrip("\n").split("\t"):
                m = re.match(r"^([a-z])\.Sequence$", col)
                if m:
                    slots.append(m.group(1))
    except OSError:
        pass
    return sorted(slots)


def read_correction_paired_summary(rbase):
    """读取 correction_paired_summary.txt（主表导出时自动生成的成对标记统计）。
    Read correction_paired_summary.txt (auto-generated with the master table)."""
    path = os.path.join(rbase, "2.labeled_result", "correction_paired_summary.txt")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


def read_two_part_summary(rbase):
    """读取 two_part_report_summary.txt → [(子目录, 文件名, 行数, 摘要行), ...]。
    先 2.labeled_result（标记样本），回退 1.unlabeled_result（无标记样本；
    v1.0 用户口径：无标记两分法报告落在 report 路径）。
    Read two_part_report_summary.txt -> [(subdir, file, rows, headline), ...].
    Checks 2.labeled_result first (labeled samples), then 1.unlabeled_result
    (unlabeled samples; v1.0 ruling: the unlabeled two-part reports live in
    the report folder)."""
    for sub in ("2.labeled_result", "1.unlabeled_result"):
        path = os.path.join(rbase, sub, "two_part_report_summary.txt")
        if not os.path.exists(path):
            continue
        out = []
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                next(fh, None)
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3 and parts[0]:
                        out.append((sub, parts[0],
                                    int(parts[1]) if parts[1].isdigit() else 0,
                                    parts[2]))
        except OSError:
            return []
        return out
    return []


def twopart_title(part, lang):
    """两分法面板标题（part 文件名 → 中/英标题）。
    Title of a two-part panel (file name -> zh/en title)."""
    if part.startswith("1."):
        return "注释校正报告 annotation correction" if lang == "zh" \
            else "annotation-correction report"
    return "新肽段竞争报告 novel peptide competition" if lang == "zh" \
        else "novel-peptide competition report"


def read_pdf_on_flag(rbase):
    """读取 result_*/pdf_on_flag.txt（"PDF_ON=1"/"PDF_ON=0"），无文件返回 ""。
    Read result_*/pdf_on_flag.txt ("PDF_ON=1"/"PDF_ON=0"); "" when absent."""
    path = os.path.join(rbase, "pdf_on_flag.txt")
    if not os.path.exists(path):
        return ""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            line = fh.readline().strip()
            return "1" if line.endswith("1") else "0"
    except OSError:
        return ""


MERGE5 = {
    "chimera_same_charge": "混合谱(同价) chimera same charge",
    "chimera_different_charge": "异价嵌合 chimera different charge",
    "n_terminal_non_specific_cleavage": "N端非特异 N-term",
    "different_modification": "修饰改变 Modification",
    "modification_location_changed": "修饰改变 Modification",
    "same_weight_aa_pairs": "序列替换 Substitution",
    "peptide_length_less_than_7aa": "序列替换 Substitution",
    "aa_order_change": "序列替换 Substitution",
    "sap": "sap",
}


def count_merge5_errors(rbase):
    """从主表 diff_type 统计五类合并错误分布（谱图数）。
    Count the five merged error categories (spectra) from the master-table diff_type."""
    master = os.path.join(rbase, "2.labeled_result", "master_table.tsv")
    if not os.path.exists(master):
        return None
    cnt = Counter()
    n = 0
    try:
        with open(master, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            iDt = header.index("diff_type") if "diff_type" in header else -1
            if iDt < 0:
                return None
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if iDt < len(parts):
                    dt = parts[iDt].strip()
                    if dt:
                        n += 1
                        cnt[MERGE5.get(dt, "其他 Other")] += 1
    except (OSError, StopIteration):
        pass
    return {"total": n, "categories": dict(cnt.most_common())}


def count_strategy_attribution(rbase):
    """从主表逐行判定哪一策略报告了被替换的候选（按序列精确匹配合并表策略列）。
    Attribute corrected candidates to strategies (exact sequence match against the
    merge-table strategy columns). Returns a Counter of strategy letters."""
    master = os.path.join(rbase, "2.labeled_result", "master_table.tsv")
    merge = os.path.join(rbase, "merge", "1.combine_result.spectra")
    if not os.path.exists(master) or not os.path.exists(merge):
        return None
    # 读合并表策略列（letter -> row_idx -> sequence）
    slot_seqs = {}
    slot_order = []
    try:
        with open(merge, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            seq_idx = {}
            for col in hdr:
                m = re.match(r"^([a-z])\.Sequence$", col)
                if m:
                    seq_idx[m.group(1)] = hdr.index(col)
                    slot_order.append(m.group(1))
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                row = {}
                for l, idx in seq_idx.items():
                    row[l] = parts[idx] if idx < len(parts) else ""
                slot_seqs[parts[0]] = row
    except OSError:
        return None
    # 读主表判校正行
    cnt = Counter()
    try:
        with open(master, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            iFn = hdr.index("File_Name")
            iCr = hdr.index("corrected_sequence")
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                dt_idx = hdr.index("diff_type")
                if dt_idx >= len(parts) or not parts[dt_idx].strip():
                    continue
                fn = parts[iFn]
                cs = parts[iCr] if iCr < len(parts) else ""
                row = slot_seqs.get(fn, {})
                attributed = False
                for l in slot_order:
                    if row.get(l, "") == cs and cs:
                        cnt[l] += 1
                        attributed = True
                        break
                if not attributed:
                    cnt["(unattributed)"] += 1
    except (OSError, ValueError):
        pass
    return dict(cnt.most_common())


def patch_config(fname, key, value, base_override=""):
    """补写 config_step1/2/3 的一个键：删除已有活动键行后按值行数写入
    （多行值＝列表键，如 search_file/mgf，每行一个）。已有活动键则替换，否则追加。
    v1.0：base_override 允许前端指定绝对路径写入位置（config_predictor
    写到 result_dir 而非包根）。
    Patch one key of config_step1/2/3: drop the existing active key lines, then
    write one line per value. v1.0: base_override lets the frontend pin
    the write location (config_predictor goes to result_dir)."""
    base = base_override if base_override and os.path.isdir(base_override) else (CURRENT_RESULT or CLIENT_ROOT)
    if fname not in ("config_step1.txt", "config_step2.txt", "config_step3.txt", "config_predictor.txt"):
        return False, "非法的参数文件名\ninvalid config file name: " + fname
    path = os.path.join(base, fname)
    if not os.path.exists(path):
        # config_predictor.txt / config_step3.txt 可能首次写入——创建空文件。
        # v1.0：config_step3 加入创建白名单——它由管线在第 2 步后才生成，
        # W8 预写（weights_B/threshold_B）发生在任何运行之前；管线重建该文件
        # 时会保留既有用户键（v1.0 old_* 回填语义），预写即生效。
        # config_predictor.txt / config_step3.txt may be first-write - create
        # an empty file. v1.0: config_step3 joins the create list - it is
        # generated by the pipeline only after step 2, while the W8 pre-write
        # (weights_B/threshold_B) happens before any run; regeneration carries
        # existing user keys forward (v1.0 old_* semantics), so the
        # pre-write takes effect.
        if fname in ("config_predictor.txt", "config_step3.txt"):
            try:
                open(path, "w", encoding="utf-8").close()
            except OSError:
                return False, "无法创建文件\ncannot create: " + path
        else:
            return False, "文件不存在\nfile not found: " + fname
    try:
        text, enc = read_config_text(path)
    except OSError as e:
        return False, "无法读取文件\ncannot read the file: %s" % e
    if text is None:
        return False, "文件编码无法识别\nfile encoding not recognized"
    lines = text.splitlines()
    # MiniMax P0-2（2026-09-26）：key 脱换行；多行值按行拆分后每行值也脱换行
    key = _sanitize_inline_field(key)
    values = [_sanitize_inline_field(v) for v in value.replace("\r", "\n").split("\n") if v.strip()]
    kept = []
    for l in lines:
        s = l.strip()
        if s and not s.startswith("#") and s.split("\t", 1)[0].strip() == key:
            continue
        kept.append(l)
    kept.extend(key + "\t" + v for v in values)
    # 杨翰反馈②：Windows 文件句柄残留导致覆写失败——先删再写保干净
    # tester feedback: delete before write to avoid Windows file handle issues
    try:
        os.remove(path)
    except OSError:
        pass
    with open(path, "w", encoding="utf-8" if enc == "utf-8-sig" else enc, newline="") as fh:
        fh.write("\n".join(kept) + "\n")
    return True, "已写入 %s = %s\nwritten: %s = %s" % (key, " / ".join(values), key, " / ".join(values))


def list_config_keys(fname):
    """列出参数文件中的键（含被注释的占位键），供停点续接页面自动生成键列表。
    List the keys of a parameter file (including commented placeholder keys) so
    the stop-point resume page can generate its key list automatically."""
    base = CURRENT_RESULT or CLIENT_ROOT
    # P0-7（2026-09-26，GLM 审查）：文件名白名单——os.path.join 对绝对路径/
    # ..\..\ 直接透传，读端此前可任意读本地文本（写端 patch_config 有同名
    # 白名单，读端缺失）
    if fname not in ("config_step1.txt", "config_step2.txt", "config_step3.txt", "config_predictor.txt"):
        return False, "非法的参数文件名\ninvalid config file name: " + fname, []
    path = os.path.join(base, fname)
    if not os.path.exists(path):
        return False, "文件不存在\nfile not found: " + fname, []
    try:
        text, _enc = read_config_text(path)
    except OSError as e:
        return False, "无法读取文件\ncannot read the file: %s" % e, []
    if text is None:
        return False, "文件编码无法识别\nfile encoding not recognized", []
    keys = {}
    order = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        active = not s.startswith("#")
        body = s[1:].strip() if not active else s
        kv = body.split("\t", 1)
        if len(kv) != 2:
            continue
        key = kv[0].strip()
        if not key:
            continue
        val = kv[1].strip()
        if key not in keys:
            keys[key] = {"key": key, "value": val, "active": active}
            order.append(key)
        elif active and not keys[key]["active"]:
            keys[key] = {"key": key, "value": val, "active": active}
    return True, "", [keys[k] for k in order]


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json; charset=utf-8", extra_headers=None):
        data = body if isinstance(body, bytes) else body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        # v1.0（2026-10-03，客户端实弹）：此前不发缓存头，浏览器启发式缓存
        # 旧 index.html——升级包后用户不强制刷新就仍是旧页面（新按钮"点了没反
        # 应"假象）。本地服务全响应 no-store，升级即生效。
        # v1.0 (2026-10-03, live client feedback): no cache headers were sent,
        # so browsers heuristically cached the old index.html - after an upgrade
        # the stale page kept serving until a hard refresh. Serve everything
        # no-store; upgrades now take effect on a plain refresh.
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra_headers or {}).items():
            self.send_header(k, v)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        # P0-9（2026-09-26，GLM 审查）/G-4（空 Host 也拒绝）：Host 校验封死
        # DNS rebinding 远程利用面——只收 127.0.0.1:8780 / localhost:8780
        host = self.headers.get("Host", "")
        # MiniMax P1-8（2026-09-26）：无端口 Host 也归一为 :8780 判定（只收本机端口）
        if host in ("127.0.0.1", "localhost"):
            host = "%s:%d" % (host, PORT)
        if host not in ("127.0.0.1:%d" % PORT, "localhost:%d" % PORT):
            self._send(403, json.dumps({"error": "forbidden host"}))
            return
        path = urllib.parse.urlparse(self.path).path
        if path in ("/", "/index.html"):
            self._send(200, get_index_bytes(), "text/html; charset=utf-8")
        elif path == "/api/version":
            # v1.0：服务端自报家门（磁盘 VERSION＋运行解释器），页脚显示；
            # "磁盘已升级但进程是旧的"从此一眼可辨
            # v1.0: server self-identification (disk VERSION + running
            # interpreter), shown in the page footer
            self._send(200, json.dumps(
                {"version": _pkg_version(),
                 "python": sys.executable,
                 "started": _SRV_STARTED,
                 "picker": _ps_probe()}, ensure_ascii=False))
        elif path == "/api/browse_state":
            self._send(200, json.dumps({"path": BROWSE_DIR or "",
                                        "active": bool(BROWSE_DIR)},
                                       ensure_ascii=False))
        elif path == "/api/report":
            # 结果报告：按当前结果文件夹实时生成的双语 HTML（可整页新窗口打开）
            # results report: bilingual HTML generated live from the current result folder
            try:
                rep_html = build_report_html()
            except Exception as e:  # noqa
                self._send(500, json.dumps({"error": "report generation failed: %s" % e}))
                return
            self._send(200, rep_html, "text/html; charset=utf-8")
        elif path == "/api/status":
            running = proc is not None and proc.poll() is None
            # 运行结束时回填历史记录的结束时间与退出码，并快照模块时间线
            # when the run has ended, back-fill its history entry and snapshot the module timeline
            if not running and RUNS and RUNS[-1]["ended"] == "" and proc is not None:
                RUNS[-1]["ended"] = time.strftime("%Y-%m-%d %H:%M:%S")
                RUNS[-1]["exit_code"] = proc.poll()
            tl_now = current_timeline()
            cur_phase = ""
            if running and tl_now:
                last = tl_now[-1]
                cur_phase = last["name_zh"] + " / " + last["name_en"]
            self._send(200, json.dumps({
                "running": running,
                "pid": proc.pid if running else None,
                "exit_code": None if running else (proc.poll() if proc is not None else None),
                "client_root": CLIENT_ROOT,
                "config": CURRENT_CFG or "",
                "result_dir": CURRENT_RESULT or "",
                "rscript": RSCRIPT,
                "runs": RUNS[-10:][::-1],
                "timeline": tl_now or [],
                "current_phase": cur_phase,
            }, ensure_ascii=False))
        elif path == "/spectrum.js":
            # 谱图浏览器脚本（静态；改文件后重启服务生效）
            # the spectrum-viewer script (static; restart to reload)
            p = os.path.join(os.path.dirname(__file__), "spectrum.js")
            try:
                with open(p, "rb") as fh:
                    self._send(200, fh.read(), "application/javascript; charset=utf-8")
            except OSError:
                self._send(404, "spectrum.js missing")
        elif path == "/api/strategies":
            # 当前结果文件夹的策略槽列表（a–z 字母，来自合并表列名）
            # strategy slots of the current result folder (letters from the merge table)
            slots = get_active_strategies()
            self._send(200, json.dumps({"ok": True, "slots": slots,
                                        "n": len(slots)}, ensure_ascii=False))
        elif path == "/api/modmass":
            # AA/修饰质量表（谱图浏览器理论离子计算用）
            # the AA/modification mass table (for theoretical-ion computation)
            mm = load_modmass()
            self._send(200, json.dumps({"ok": bool(mm["aa"]),
                                        "aa": mm["aa"], "mods": mm["mods"]}))
        elif path == "/api/spectra_list":
            ok, files = spectra_list()
            self._send(200, json.dumps({"ok": ok, "files": [f[0] for f in files],
                                        "n": len(files)}, ensure_ascii=False))
        elif path == "/api/spectrum":
            # 单张谱：MGF 峰列表＋候选串＋结论（conflict 优先回退 report_A）
            # one spectrum: MGF peaks + candidates + conclusion
            # （urlparse().path 不含查询串，参数在 query 里）
            # (urlparse().path excludes the query string; params live in query)
            qs = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
            fn = (qs.get("file", [""])[0]).strip()
            if not fn:
                self._send(400, json.dumps({"ok": False, "error": "file missing"}))
                return
            peaks, title = find_spectrum(fn)
            if peaks is None:
                self._send(404, json.dumps(
                    {"ok": False, "error": "spectrum not found in configured MGFs"},
                    ensure_ascii=False))
                return
            info = spectrum_info(fn) or {}
            self._send(200, json.dumps({
                "ok": True, "title": title, "peaks": peaks,
                "candidates": info.get("cands", []),
                "winner_index": info.get("winner_index", -1),
                "conclusion": info.get("conclusion", "")}, ensure_ascii=False))
        elif path == "/api/pdfs":
            # 谱图 PDF 清单（结果文件夹/软件根/web_shell[\pdf]/指定文件夹，深度 3）
            # the spectrum-PDF list (result folder / client root / web_shell[\pdf]
            # / a user-supplied folder, depth-limited)
            qs = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
            extra = (qs.get("dir", [""])[0]).strip()
            if extra and os.path.isdir(extra):
                rp = os.path.realpath(extra)
                if rp not in _EXTRA_PDF_ROOTS:
                    _EXTRA_PDF_ROOTS.append(rp)
                    del _EXTRA_PDF_ROOTS[:-10]     # 只记最近 10 个手填目录
            pdfs = find_pdfs(pdf_scan_dirs(extra or None))
            self._send(200, json.dumps({"ok": True, "pdfs": pdfs}, ensure_ascii=False))
        elif path == "/pdf":
            # 服务一个 PDF（本机工具，供结果页 iframe 窗口内查看/打印）
            # serve one PDF (local tool; embedded in the report page iframe)
            # v1.0（2026-09-26，审查 P0-1/H2）：路径白名单——realpath 必须落在
            # 客户端根目录（或结果目录）之下，禁止读取任意本地 PDF；响应头加
            # nosniff。B-13（2026-09-26）：报告页 iframe 需内联查看，故不加
            # Content-Disposition: attachment（注释与实现对齐）。
            # v1.0 (2026-09-26, review P0-1/H2): path whitelist - the realpath
            # must stay under the client root (or a result folder); no arbitrary
            # local PDF reads. nosniff added. B-13 (2026-09-26): the report-page
            # iframe views the PDF inline, so no Content-Disposition: attachment.
            origin = self.headers.get("Origin", "")
            if origin and urllib.parse.urlparse(origin).netloc not in (
                    "127.0.0.1:%d" % PORT, "localhost:%d" % PORT):
                self._send(403, json.dumps({"error": "forbidden origin"}))
                return
            # 原始 UTF-8 查询串（未百分号编码）经 http.server 的 latin-1 解码会
            # 变成乱码——先按 latin-1 字节还原再按 UTF-8 解码（已编码路径不受影响）
            qs = urllib.parse.parse_qs(
                urllib.parse.urlparse(self.path.encode("latin-1").decode("utf-8", "replace")).query)
            p = (qs.get("path", [""])[0]).strip()
            allowed_roots = [CLIENT_ROOT, _result_base(CLIENT_ROOT)]
            for r in ([CURRENT_RESULT] if CURRENT_RESULT else []) + \
                     ([BROWSE_DIR] if BROWSE_DIR else []) + \
                     list(_EXTRA_PDF_ROOTS):
                if r and r not in allowed_roots:
                    allowed_roots.append(r)
            # 自查整改（2026-09-26）：白名单根与请求路径同口径 realpath+normcase，
            # 防 junction/符号链接前缀失配；并入 CURRENT_RESULT/浏览目录/手填目录，
            # 使"能扫描就能打开"（清单与白名单共用同一份根集合）
            allowed_roots = [os.path.normcase(os.path.realpath(r))
                             for r in allowed_roots]
            try:
                rp = os.path.realpath(p)
                rn = os.path.normcase(os.path.normpath(rp))
                inside = any(rn == os.path.normcase(os.path.normpath(r)) or
                             rn.startswith(os.path.normcase(os.path.normpath(r)) + os.sep)
                             for r in allowed_roots)
            except OSError:
                inside = False
            if (not p.lower().endswith(".pdf") or not os.path.isfile(rp) or not inside):
                self._send(404, "pdf not found")
                return
            try:
                with open(rp, "rb") as fh:
                    data = fh.read()
            except OSError as e:
                self._send(500, "cannot read pdf: %s" % e)
                return
            self._send(200, data, "application/pdf",
                       {"X-Content-Type-Options": "nosniff"})
        elif path == "/api/log":
            tail = ""
            if os.path.exists(LOG_FILE):
                try:
                    with open(LOG_FILE, "r", encoding="utf-8", errors="replace") as fh:
                        tail = "".join(deque(fh, maxlen=300))
                except OSError:
                    tail = ""
            self._send(200, json.dumps({"log": tail}, ensure_ascii=False))
        elif path == "/api/config_keys":
            q = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
            ok, msg, keys = list_config_keys((q.get("file") or [""])[0])
            self._send(200, json.dumps({"ok": ok, "msg": msg, "keys": keys}, ensure_ascii=False))
        elif path == "/api/detect":
            try:
                global _DETECT_CACHE, _DETECT_TS
                now = time.time()
                if _DETECT_CACHE is None or now - _DETECT_TS > 60:
                    _DETECT_CACHE = detect_environment()
                    _DETECT_TS = now
                self._send(200, json.dumps(_DETECT_CACHE, ensure_ascii=False))
            except Exception as e:  # 探测失败返回 500 JSON 而非掐断连接
                self._send(500, json.dumps({"ok": False, "error": str(e)}, ensure_ascii=False))
            except Exception as e:
                self._send(500, json.dumps({"ok": False, "error": str(e)}, ensure_ascii=False))
        else:
            self._send(404, json.dumps({"error": "not found"}))

    def do_POST(self):
        # P0-9（2026-09-26，GLM 审查）/G-4（空 Host 也拒绝）：与 do_GET 同口径
        host = self.headers.get("Host", "")
        # MiniMax P1-8（2026-09-26）：无端口 Host 也归一为 :8780 判定（只收本机端口）
        if host in ("127.0.0.1", "localhost"):
            host = "%s:%d" % (host, PORT)
        if host not in ("127.0.0.1:%d" % PORT, "localhost:%d" % PORT):
            self._send(403, json.dumps({"error": "forbidden host"}))
            return
        # 2026-09-26 冒烟修复（v1.0）：proc 必须一并声明 global——P0-A 在
        # /api/stop 分支内新增 `proc = None` 赋值后，proc 在 do_POST 内变为
        # 局部名，`p = proc` 读取先于赋值触发 UnboundLocalError（停止按钮
        # 全路径崩溃）。冒烟测试 POST /api/stop 捕获。
        # 2026-09-26 smoke fix (v1.0): proc must also be declared global -
        # after P0-A added `proc = None` inside the /api/stop branch, proc
        # became function-local and the `p = proc` read raised
        # UnboundLocalError (the stop button crashed on every call). Caught by
        # the smoke test POSTing /api/stop.
        global CURRENT_CFG, CURRENT_RESULT, proc
        path = urllib.parse.urlparse(self.path).path
        origin = self.headers.get("Origin", "")
        if origin and urllib.parse.urlparse(origin).netloc not in (
                "127.0.0.1:%d" % PORT, "localhost:%d" % PORT):
            self._send(403, json.dumps({"error": "forbidden origin"}))
            return
        try:
            n = int(self.headers.get("Content-Length", 0) or 0)
            if n > 10 * 1024 * 1024:
                self._send(413, json.dumps({"error": "payload too large"}))
                return
            body = json.loads(self.rfile.read(n).decode("utf-8") if n else "{}")
        except (ValueError, UnicodeDecodeError):
            self._send(400, json.dumps({"error": "invalid json"}))
            return
        if path == "/api/pick_file":
            # 排队③（2026-09-26）：改 POST——GET 会允许任意站点触发文件对话框
            # （跨站对话框轰炸）；POST 经 Host/Origin 校验后同进程锁内弹框
            kind = body.get("kind") or "file"
            ok, msg = pick_path(kind)
            # v1.0：响应带 route（ctypes/ps/ctypes-failed->ps）——远程
            # 排障一眼分清走的哪条路线
            self._send(200, json.dumps({"ok": ok, "path": msg if ok else "",
                                        "msg": msg, "route": _PICK_ROUTE[0]},
                                       ensure_ascii=False))
        elif path == "/api/detect":
            # v1.0：POST body {extra_roots:[...]}＝用户在界面追加的扫描根
            # 目录（最多 4 个、须存在；detect_environment 内去重并置于内置
            # 根之前）。带 extras 绕过 60s 缓存（用户主动扫描，多等一轮可接
            # 受）；不带的 POST 等价 GET（走缓存）。
            # v1.0: POST body {extra_roots:[...]} = user-added scan roots
            # from the UI (<=4, must exist; detect_environment dedupes and
            # puts them before the built-in roots). With extras the 60 s cache
            # is bypassed; without extras the POST equals the GET (cached).
            try:
                extra = body.get("extra_roots") or None
                if extra:
                    self._send(200, json.dumps(
                        {"ok": True, "extra": True,
                         **detect_environment(extra_roots=extra)},
                        ensure_ascii=False))
                else:
                    global _DETECT_CACHE, _DETECT_TS
                    now = time.time()
                    if _DETECT_CACHE is None or now - _DETECT_TS > 60:
                        _DETECT_CACHE = detect_environment()
                        _DETECT_TS = now
                    self._send(200, json.dumps(_DETECT_CACHE, ensure_ascii=False))
            except Exception as e:
                self._send(500, json.dumps({"ok": False, "error": str(e)},
                                           ensure_ascii=False))
        elif path == "/api/resolve_path":
            # v1.0：浏览器原生选框只给文件/文件夹名，此端点在常见安装
            # 根目录搜索同名条目返回完整路径（找不到则前端填名让用户补全）
            # v1.0: the browser-native picker only yields the file/folder
            # NAME; this endpoint searches the common install roots for a
            # matching entry and returns the full path (if not found the
            # frontend fills the name and the user completes the prefix)
            name = (body.get("name") or "").strip()
            if not name or len(name) > 260:
                self._send(200, json.dumps({"ok": False, "path": "",
                                            "msg": "invalid name"}, ensure_ascii=False))
                return
            hits = []
            roots = [r for r in _DETECT_ROOTS if os.path.isdir(r)]
            for root in roots:
                try:
                    for entry in os.scandir(root):
                        if entry.name == name:
                            hits.append(entry.path)
                except OSError:
                    continue
                if len(hits) >= 5:
                    break
            if not hits:
                for root in roots:
                    try:
                        for sub in os.scandir(root):
                            if sub.is_dir(follow_symlinks=False):
                                cand = os.path.join(sub.path, name)
                                if os.path.exists(cand):
                                    hits.append(cand)
                                    if len(hits) >= 5:
                                        break
                    except OSError:
                        continue
            if len(hits) == 1:
                self._send(200, json.dumps({"ok": True, "path": hits[0],
                                            "msg": ""}, ensure_ascii=False))
            elif len(hits) > 1:
                self._send(200, json.dumps({"ok": True, "path": hits[0],
                                            "alternatives": hits,
                                            "msg": "multiple"}, ensure_ascii=False))
            else:
                self._send(200, json.dumps({"ok": False, "path": "",
                                            "msg": "not found in scan roots"}, ensure_ascii=False))
        elif path == "/api/list_dir":
            # v1.0：内嵌文件浏览器——列出指定目录的子目录（与可选文件）
            # v1.0: in-page file browser - list subdirectories (and optional files)
            dpath = (body.get("path") or "").strip()
            show_files = bool(body.get("files"))
            if not dpath:
                # 返回可用盘符
                import string
                drives = []
                for letter in string.ascii_uppercase:
                    drv = letter + ":\\"
                    if os.path.isdir(drv):
                        drives.append(drv)
                self._send(200, json.dumps({"ok": True, "drives": drives,
                                            "dirs": [], "files": [], "path": ""},
                                           ensure_ascii=False))
                return
            dpath = os.path.normpath(dpath)
            if not os.path.isdir(dpath):
                self._send(200, json.dumps({"ok": False, "msg": "not a dir: %s" % dpath,
                                            "dirs": [], "files": [], "path": dpath},
                                           ensure_ascii=False))
                return
            dirs = []
            files = []
            try:
                for entry in os.scandir(dpath):
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            dirs.append(entry.name)
                        elif show_files and entry.is_file():
                            files.append(entry.name)
                    except OSError:
                        continue
            except OSError as e:
                self._send(200, json.dumps({"ok": False, "msg": str(e),
                                            "dirs": [], "files": [], "path": dpath},
                                           ensure_ascii=False))
                return
            parent = os.path.dirname(dpath)
            self._send(200, json.dumps({"ok": True, "dirs": sorted(dirs),
                                        "files": sorted(files), "path": dpath,
                                        "parent": parent if parent != dpath else ""},
                                       ensure_ascii=False))
        elif path == "/api/set_browse":
            # 结果浏览模式（2026-09-26）：选定另一台电脑的结果项目目录/结果目录，
            # 报告、谱图、PDF 全部只读改从该目录读取
            ok, msg = _set_browse(body.get("path", ""))
            self._send(200, json.dumps({"ok": ok,
                                        "path": msg if ok else "",
                                        "msg": ("" if ok else msg)},
                                       ensure_ascii=False))
        elif path == "/api/clear_browse":
            _set_browse("")
            self._send(200, json.dumps({"ok": True}, ensure_ascii=False))
        elif path == "/api/generate_config":
            fields = body.get("fields", {})
            if not isinstance(fields, dict):
                self._send(400, json.dumps({"error": "fields must be an object"}))
                return
            txt = build_step1_config(fields)
            result_dir = (fields.get("result_dir") or "").strip()
            # 结果保存文件夹为空时退回客户端目录
            # an empty result folder falls back to the client root
            if not result_dir:
                result_dir = CLIENT_ROOT
            try:
                os.makedirs(result_dir, exist_ok=True)
                p = os.path.join(result_dir, "config_step1.txt")
                with open(p, "w", encoding="utf-8", newline="") as fh:
                    fh.write(txt)
                # 写后验证：文件确实存在才返回 ok
                if not os.path.isfile(p):
                    self._send(200, json.dumps({"ok": False,
                                                "msg": ("写入后验证失败——文件未出现在磁盘上："
                                                        + p + "\n请检查磁盘是否可写（U 盘写保护？）"
                                                        + "\npost-write verification failed: " + p
                                                        + " (check disk write access)")},
                                               ensure_ascii=False))
                    return
            except OSError as e:
                self._send(200, json.dumps({"ok": False,
                                            "msg": "无法写入结果文件夹\ncannot write result folder: %s" % e},
                                           ensure_ascii=False))
                return
            CURRENT_CFG = p
            CURRENT_RESULT = result_dir
            self._send(200, json.dumps({"ok": True, "path": p, "result_dir": result_dir,
                                        "content": txt}, ensure_ascii=False))
        elif path == "/api/generate_step2":
            strategies = body.get("strategies", [])
            if not isinstance(strategies, list) or not strategies:
                self._send(200, json.dumps({"ok": False,
                                            "msg": "策略列表为空\nstrategy list is empty"},
                                           ensure_ascii=False))
                return
            txt = build_step2_config(strategies, (body.get("msfragger_xmx") or "").strip())
            base = (body.get("result_dir") or "").strip() or CURRENT_RESULT or CLIENT_ROOT
            # config_step2 必须带 base_config 键（user_step2 靠它定位步骤 1 参数；
            # 网页版此前只写 search_file 行，缺键使 read_param_file(NA) 崩溃）
            # the step-2 file must carry base_config (user_step2 resolves the
            # step-1 parameters through it; the old web output only wrote the
            # search_file lines and the missing key crashed read_param_file(NA))
            cfg1 = os.path.join(base, "config_step1.txt")
            txt = ("base_config\t" + cfg1.replace("\\", "/") + "\n\n" + txt)
            try:
                os.makedirs(base, exist_ok=True)
                p = os.path.join(base, "config_step2.txt")
                with open(p, "w", encoding="utf-8", newline="") as fh:
                    fh.write(txt)
            except OSError as e:
                self._send(200, json.dumps({"ok": False,
                                            "msg": "无法写入\ncannot write: %s" % e},
                                           ensure_ascii=False))
                return
            self._send(200, json.dumps({"ok": True, "path": p, "content": txt},
                                       ensure_ascii=False))
        elif path == "/api/run":
            # 李治反馈：全局变量在多线程 HTTP 服务中不可靠——前端直接传 config 路径
            # do_POST 顶部已有 global CURRENT_CFG，此处直接赋值
            _cfg_from_fe = (body.get("config_path") or "").strip()
            if _cfg_from_fe and os.path.isfile(_cfg_from_fe):
                CURRENT_CFG = _cfg_from_fe
            ok, msg = run_pipeline()
            self._send(200, json.dumps({"ok": ok, "msg": msg}, ensure_ascii=False))
        elif path == "/api/stop":
            # P1-2/P0-A（2026-09-25）：锁内取句柄快照。Windows 上 terminate()
            # 是硬杀父进程——先 terminate 则 taskkill /T（进程树）条件恒假、
            # 子进程（Searcher/java/PSOCK worker）全部残留；故 Windows 先
            # taskkill 杀树，POSIX 保持温和→强制阶梯。成功后置 proc=None。
            with _PROC_LOCK:
                p = proc
            tk_msg = ""
            if p is not None and p.poll() is None:
                if os.name == "nt":
                    # §3.1（2026-09-26）：检查 taskkill 返回码——权限/跨用户
                    # 失败时不再静默空等，把失败信息回给前端。
                    tk = subprocess.run(["taskkill", "/T", "/F", "/PID", str(p.pid)],
                                        capture_output=True)
                    if tk.returncode != 0:
                        tk_msg = ("（taskkill 失败 rc=%d：%s）"
                                  % (tk.returncode,
                                     (tk.stderr or tk.stdout or b"").decode("mbcs", "replace").strip()))
                    try:
                        p.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        try:
                            p.kill()
                            p.wait(timeout=5)
                        except Exception:
                            pass
                else:
                    try:
                        p.terminate()
                    except OSError:
                        pass
                    try:
                        p.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        try:
                            p.kill()
                            p.wait(timeout=5)
                        except Exception:
                            pass
                leftover = p.poll() is None
                if not leftover:
                    with _PROC_LOCK:
                        proc = None
                self._send(200, json.dumps(
                    {"ok": not leftover,
                     "msg": ("已停止\nstopped" if not leftover
                             else "已请求强制结束，若残留请手动结束 pid %d\nforce-kill requested; if it remains, end pid %d manually" % p.pid) + tk_msg},
                    ensure_ascii=False))
            else:
                self._send(200, json.dumps({"ok": True, "msg": "无运行中的流程\nno running process"}, ensure_ascii=False))
        elif path == "/api/fill":
            ok, msg = patch_config(body.get("file", ""), body.get("key", ""), body.get("value", ""),
                                   body.get("base", ""))
            self._send(200, json.dumps({"ok": ok, "msg": msg}, ensure_ascii=False))
        elif path == "/api/load_config":
            # v99l：本机工具（localhost + Host 校验），允许从任意路径加载
            # config——保留 basename 白名单（只收 config_step1/2/3.txt），
            # 弃"join 固定 base"限制（用户从文件浏览器选到哪就从哪读）
            p = (body.get("path") or "").strip()
            fname = os.path.basename(p.replace("\\", "/"))
            if fname not in ("config_step1.txt", "config_step2.txt", "config_step3.txt", "config_predictor.txt"):
                self._send(200, json.dumps({"ok": False,
                                            "msg": "invalid config file name: " + p,
                                            "fields": {}}, ensure_ascii=False))
                return
            if not os.path.isfile(p):
                self._send(200, json.dumps({"ok": False,
                                            "msg": "file not found: " + p,
                                            "fields": {}}, ensure_ascii=False))
                return
            ok, msg, fields = read_config_fields(p)
            self._send(200, json.dumps({"ok": ok, "msg": msg, "fields": fields}, ensure_ascii=False))
        elif path == "/api/mkdir":
            # 文件浏览器：新建文件夹
            p = (body.get("path") or "").strip()
            name = (body.get("name") or "").strip()
            if not p or not name:
                self._send(200, json.dumps({"ok": False, "msg": "path/name required"}, ensure_ascii=False))
                return
            target = os.path.join(p, name)
            try:
                os.makedirs(target, exist_ok=True)
                self._send(200, json.dumps({"ok": True, "path": target}, ensure_ascii=False))
            except OSError as e:
                self._send(200, json.dumps({"ok": False, "msg": str(e)}, ensure_ascii=False))
        elif path == "/api/rename":
            # 文件浏览器：重命名文件夹
            old_p = (body.get("old") or "").strip()
            new_p = (body.get("new") or "").strip()
            if not old_p or not new_p:
                self._send(200, json.dumps({"ok": False, "msg": "old/new required"}, ensure_ascii=False))
                return
            try:
                os.rename(old_p, new_p)
                self._send(200, json.dumps({"ok": True, "path": new_p}, ensure_ascii=False))
            except OSError as e:
                self._send(200, json.dumps({"ok": False, "msg": str(e)}, ensure_ascii=False))
        elif path == "/api/check_config":
            # 李治反馈②：Run 前检查 result_dir 内是否有 config_step1.txt
            rd = (body.get("result_dir") or "").strip()
            cfg1 = os.path.join(rd, "config_step1.txt") if rd else ""
            ok = bool(cfg1 and os.path.isfile(cfg1))
            self._send(200, json.dumps({"ok": ok,
                                        "path": cfg1}, ensure_ascii=False))
        else:
            self._send(404, json.dumps({"error": "not found"}))

    def log_message(self, fmt, *args):
        pass  # quiet


def _pkg_version():
    """v1.0：读包根 VERSION 文件（网页页脚与 /api/version 用——远程排障
    一眼分清'磁盘包已升级但进程还是旧'的情形）。
    v1.0: read the package VERSION file (shown in the page footer and
    served by /api/version - makes 'disk upgraded but process stale'
    instantly visible during remote debugging)."""
    try:
        with open(os.path.join(CLIENT_ROOT, "VERSION"), "rb") as fh:
            return fh.read().decode("utf-8", "replace").strip()
    except OSError:
        return "unknown"


def _kill_stale_listener():
    """v1.0（2026-10-03，客户端实弹）：升级包后旧 server 进程仍占 8780，
    新进程 bind 失败而浏览器一直连旧进程（页面新、API 旧、"按钮没反应"）。
    按端口精准结束旧监听（Get-NetTCPConnection → OwningProcess），只动
    8780 监听者，不碰其他 python。
    v1.0 (2026-10-03, live): after an upgrade the OLD server keeps port
    8780; the new one fails to bind and the browser keeps talking to the
    stale process (new page, old API, "dead buttons"). Kill the listener on
    8780 precisely via Get-NetTCPConnection; other python processes are
    untouched."""
    if os.name != "nt":
        return False
    ps = ("Get-NetTCPConnection -LocalPort %d -State Listen -ErrorAction "
          "SilentlyContinue | ForEach-Object { Write-Host ('stale PID ' + "
          "$_.OwningProcess); Stop-Process -Id $_.OwningProcess -Force -"
          "ErrorAction SilentlyContinue }" % PORT)
    try:
        subprocess.run(["powershell", "-NoProfile", "-Command", ps],
                       capture_output=True, timeout=30)
        return True
    except Exception:  # noqa
        return False


def main():
    try:
        srv = ThreadingHTTPServer((HOST, PORT), Handler)
    except OSError:
        # v1.0：端口被旧进程占用 → 自动清场重试一次；再失败给双语指引
        # v1.0: port held by a stale process - clear it and retry once;
        # if it still fails, print bilingual guidance
        print("!! 端口 %d 被占用（多半是升级前的旧服务进程）——尝试自动结束旧进程…"
              % PORT)
        print("!! port %d busy (likely the pre-upgrade server) - killing the "
              "stale process..." % PORT)
        _kill_stale_listener()
        time.sleep(1)
        try:
            srv = ThreadingHTTPServer((HOST, PORT), Handler)
        except OSError as e:
            print("!! 仍无法绑定端口 / still cannot bind: %s" % e)
            print("!! 请手动结束旧进程后重开本窗口：任务管理器搜 python.exe（命令行含 web\\server.py）→ 结束任务")
            print("!! manually kill the old server (Task Manager: python.exe running web\\server.py), then reopen")
            raise SystemExit(1)
    print("=== 网页版壳子 web shell（%s）===" % _pkg_version())
    print("客户端目录 client root:", CLIENT_ROOT)
    # F-2（2026-09-26）：把探测选中的 Rscript 路径与版本打出来（远程排障时
    # 一眼可见多版本共存时选了哪个）。
    # F-2 (2026-09-26): log the probed Rscript path and version, so it is
    # immediately visible which R was picked when several coexist.
    print("Rscript 路径 Rscript:", RSCRIPT)
    try:
        _rv = subprocess.run([RSCRIPT, "--version"], capture_output=True, timeout=10)
        print("R 版本 R version:", (_rv.stdout or _rv.stderr or b"").decode("utf-8", "replace").strip())
    except Exception as _e:
        print("R 版本探测失败 R version probe failed:", _e)
    print("请在浏览器打开 open: http://%s:%d" % (HOST, PORT))
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
