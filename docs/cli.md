# CLI Guide / 命令行指南

---

# 中文版

```bat
cd /d "G:\isoPSV_multiengine"
set CONFIG=G:\项目文件夹\config_step1.txt
"%RSCRIPT%" code\user_run.R
```

- 只与总入口 `code/user_run.R` 交互；`CONFIG` 环境变量指向 config_step1.txt
  绝对路径（中文路径一律走环境变量，不进命令行参数）。
- **三轮运行**：第一轮＝第一步提取无标记谱＋第二步搜库策略（手动搜库时停点
  提示填 config_step2 的 search_file）；第二轮＝合并候选表＋生成 config_step3
  （手动 pDeep2 时停点提示填 predict_file）；第三轮＝阶段 A/B 仲裁＋全部报告。
- 断点续跑：中断后重跑同一条命令即从断点续接（完成标记/指纹自动跳过）。
- 全部键值见 [parameters.md](parameters.md)；引擎工作区与项目目录须 ASCII
  路径，启动前清 LC_CTYPE/LANG/LC_ALL。
- 纯 R 库函数（MGF 解析/打分/工具）见独立 R 包
  [ReadMGF](https://github.com/Songhao-Jiang/ReadMGF)。

---

# English

```bat
cd /d "G:\isoPSV_multiengine"
set CONFIG=G:\project\config_step1.txt
"%RSCRIPT%" code\user_run.R
```

- Interact only with the master entry `code/user_run.R`; the `CONFIG`
  environment variable points at the absolute config_step1.txt path
  (non-ASCII paths always via the environment, never as CLI arguments).
- **Three rounds**: round 1 = step-1 unlabeled extraction + step-2 search
  strategies (stop point asks for config_step2 search_file in manual mode);
  round 2 = merge candidates + generate config_step3 (stop point asks for
  predict_file with manual pDeep2); round 3 = stage-A/B arbitration + all
  reports.
- Checkpoint resume: rerun the same command after an interruption; completed
  steps are skipped via done markers/fingerprints.
- See [parameters.md](parameters.md) for every key; engine workspaces and
  project directories must be ASCII paths; clear LC_CTYPE/LANG/LC_ALL first.
- Pure R library functions (MGF parsing / scoring / utilities) live in the
  standalone R package
  [ReadMGF](https://github.com/Songhao-Jiang/ReadMGF).
