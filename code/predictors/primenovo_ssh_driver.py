# -*- coding: utf-8 -*-
"""primenovo_ssh_driver.py — isoPSV 的 PrimeNovo SSH/VM 调用桥

场景：PrimeNovo 部署在 Linux 主机/VM（ctcdecode 无 Windows 原生构建），
isoPSV 在 Windows 宿主上运行。本驱动把一个 MGF 推到远端、远端跑
PrimeNovo denovo、拉回 denovo.tsv，全程 SSH 免密（公钥认证）。

用法（由 user_step2.R 的 denovo_engine_auto_run 自动调用）：
  python primenovo_ssh_driver.py --host jsh@192.168.216.128 \
      --mgf <local.mgf> --model <remote-path-or-local> --out <local.tsv> \
      [--repo ~/pi-PrimeNovo] [--python ~/miniconda3/envs/primenovo/bin/python] \
      [--remote-mgf /tmp/isopsv_pn/<basename>] [--keep-remote]

本地 model 路径自动推到远端固定位（~/isopsv_pn_model.ckpt，仅首次）；
remote-mgf 默认 /tmp/isopsv_pn/ 前缀＋原文件名；远端 tsv 拉回写到 --out。
退出码：0 成功；2 参数错误；3 远端运行失败；4 传输失败。
"""
import argparse
import os
import posixpath
import sys

import paramiko

REMOTE_DIR = "/tmp/isopsv_pn"
REMOTE_MODEL = "isopsv_pn_model.ckpt"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True, help="user@host")
    ap.add_argument("--mgf", required=True)
    ap.add_argument("--model", required=True,
                    help="本地或远端 ckpt 路径（Windows 路径＝本地推，/开头＝远端已有")
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", default="~/pi-PrimeNovo")
    ap.add_argument("--python", default="~/miniconda3/envs/primenovo/bin/python")
    ap.add_argument("--remote-mgf", default="")
    ap.add_argument("--keep-remote", action="store_true")
    return ap.parse_args()


def main():
    a = parse_args()
    if not os.path.isfile(a.mgf):
        print("[pn-ssh] local MGF not found:", a.mgf, file=sys.stderr)
        return 2
    user, _, host = a.host.partition("@")
    if not host:
        print("[pn-ssh] --host must be user@host", file=sys.stderr)
        return 2
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        c.connect(host, username=user, timeout=30)
    except Exception as exc:
        print("[pn-ssh] SSH connect failed:", exc, file=sys.stderr)
        return 4
    sftp = c.open_sftp()
    try:
        # 远端目录
        def run(cmd, timeout=3600):
            _, o, e = c.exec_command(cmd, timeout=timeout)
            rc = o.channel.recv_exit_status()
            out = o.read().decode(errors="replace").strip()
            err = e.read().decode(errors="replace").strip()
            return rc, out, err

        rc, _, err = run("mkdir -p %s" % REMOTE_DIR)
        if rc != 0:
            print("[pn-ssh] mkdir failed:", err, file=sys.stderr)
            return 4

        # ckpt：本地推（仅首次，按远端 mtime 指纹跳过）或远端路径直用
        model_remote = a.model if a.model.startswith("/") else (
            posixpath.join(REMOTE_DIR, REMOTE_MODEL))
        if not a.model.startswith("/"):
            need = True
            try:
                local_size = os.path.getsize(a.model)
                rc, out, _ = run("stat -c %%s %s 2>/dev/null" % model_remote)
                if rc == 0 and out.strip() == str(local_size):
                    need = False
            except OSError:
                pass
            if need:
                print("[pn-ssh] pushing model (one-time, %.0f MB)..." % (
                    os.path.getsize(a.model) / 1048576), file=sys.stderr)
                sftp.put(a.model, model_remote)

        # MGF 推
        mgf_remote = a.remote_mgf or posixpath.join(
            REMOTE_DIR, os.path.basename(a.mgf))
        sftp.put(a.mgf, mgf_remote)

        # 远端运行（repo 根目录 cwd；denovo.tsv 落 cwd）
        tsv_remote_name = os.path.basename(a.out)
        cmd = "cd {repo} && rm -f denovo.tsv && " \
              "timeout 7200 {py} -m PrimeNovo.PrimeNovo --mode=denovo " \
              "--peak_path={mgf} --model={model} > {dir}/run.log 2>&1; " \
              "rc=$?; [ $rc -eq 0 ] && [ -f denovo.tsv ] && " \
              "mv denovo.tsv {dir}/{tsv}; exit $rc".format(
                  repo=a.repo, py=a.python, mgf=mgf_remote,
                  model=model_remote, dir=REMOTE_DIR, tsv=tsv_remote_name)
        rc, out, err = run(cmd, timeout=7200)
        if rc != 0:
            # 回贴远端日志尾部（ 纪律）
            _, log_tail, _ = run("tail -15 %s/run.log" % REMOTE_DIR)
            print("[pn-ssh] remote run failed (exit %s):\n%s" % (
                rc, log_tail), file=sys.stderr)
            return 3

        # 拉回
        sftp.get(posixpath.join(REMOTE_DIR, tsv_remote_name), a.out)
        if not a.keep_remote:
            run("rm -f %s %s/%s" % (mgf_remote, REMOTE_DIR, tsv_remote_name))
        print("[pn-ssh] done:", a.out, file=sys.stderr)
        return 0
    except Exception as exc:
        print("[pn-ssh] transfer error:", exc, file=sys.stderr)
        return 4
    finally:
        sftp.close()
        c.close()


if __name__ == "__main__":
    sys.exit(main())
