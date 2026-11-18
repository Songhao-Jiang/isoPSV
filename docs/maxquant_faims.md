# isoPSV MaxQuant 自动化搜库方法（含 FAIMS 兼容性说明）

> 适用版本：isoPSV v2.0.178+（MaxQuant 集成路径）／ MaxQuant 2.8.1.0+
> **v2.0.185 修订**：更正 config 键位置（由工作站实弹复测发现原文档不准确）；新增 FAIMS 兼容性实测结论。

## ⚠ FAIMS 兼容性结论（工作站实弹 2026-11-18）

**Thermo Exploris 480 FAIMS 数据在 MaxQuant 2.8.1.0 上不可用**：
- 92% MS2 precursor mass = NaN（MaxQuant 无法解析该仪器的 FAIMS CV 元数据）
- 拆分 raw 后仍 28% NaN
- raw→mzML 后 first search 挂起 >1h
- **推荐替代**：MSFragger（同机已验证 1,378,952 PSM / Trypsin）

isoPSV 的 MaxQuant 集成在**非 FAIMS 数据**上正常工作。FAIMS 问题出在 MaxQuant 对特定仪器 FAIMS 元数据的解析，不是 isoPSV 的缺陷。

## 一、原理

isoPSV 通过模板替换自动化 MaxQuant 搜库——用户提供 `.raw` 文件路径，软件改写 `mqpar.xml` 模板后调 `MaxQuantCmd.exe`，搜库完成后自动转换 `msms.txt` → `.spectra` 进入竞争仲裁。

## 二、前提条件

| 项 | 要求 |
|---|---|
| MaxQuant | ≥ 2.0（本模板基于 2.8.1.0） |
| 原始数据 | Thermo `.raw` 格式；**FAIMS 数据须先确认 MaxQuant 兼容性**（见上方结论） |
| .NET | ≥ 6.0 运行时 |

## 三、isoPSV 配置（v2.0.185 更正）

**实际键读取位置**：`maxquant_params` 和 `maxquant_mode` 由 `user_step2.R` 的通用引擎循环读取（L3118-3119），来源为 config_step2 或 config_step1（step1 优先）。`maxquant_exe` 在 config_step1；`search_file_raw` 在 config_step1（网页版写入）或 config_step2。

### config_step1.txt

```
maxquant_exe    D:/Bioinformatics/MaxQuant_2.8.1.0/bin/MaxQuantCmd.exe
search_file_raw    D:/data/sample/
```

### config_step2.txt（或 config_step1，step1 优先）

```
maxquant_params    G:/path/to/maxquant_con.mqpar.xml
maxquant_mode    1
```

| 键 | 位置 | 说明 |
|---|---|---|
| `maxquant_exe` | step1 | MaxQuantCmd.exe 完整路径 |
| `search_file_raw` | step1 | `.raw` 文件夹或分号分隔的文件列表 |
| `maxquant_params` | step1 或 step2 | mqpar.xml 模板路径（默认用自带 `code/engine_templates/maxquant_con.mqpar.xml`） |
| `maxquant_mode` | step1 或 step2 | 1＝自动（0＝跳过） |

**注意**：不填 `maxquant_params` 时使用内置默认模板。`maxquant_mode=1` 时必须有 `search_file_raw` 指向 `.raw` 文件。

### 模板占位符（软件自动替换，用户不用管）

| 占位符 | 替换为 |
|---|---|
| `{FASTA}` | 用户 FASTA 路径（反斜杠 Windows 格式） |
| `{OUT}` | 结果目录 |
| `{RAWS}` | raw 文件路径 XML 块 |
| `{EXP}` | 实验名 XML 块 |
| `{FRAC}` | 分数标签 XML 块 |

### 酶切类型（自动改写）

`enzyme` 键（step1）自动改写 mqpar `<enzymes>` 块：

| isoPSV 键 | MaxQuant 名 |
|---|---|
| trypsin | Trypsin |
| lysc | LysC |
| lysn | LysN |
| argc | ArgC |
| aspn | AspN |
| gluc | GluC |
| chymotrypsin | Chymotrypsin |
| lysarginase | LysargiNase |

不支持的酶 → 该槽位跳过＋提示。

## 四、FAIMS 注意事项（仅适用于已确认兼容的仪器）

如果 MaxQuant 在您的仪器上能正确解析 FAIMS raw（msmsScans.txt 无 NaN precursor mass）：

- MaxQuant 自动按 CV 分段搜库
- `matchBetweenRuns` 保持 `False`（FAIMS 的 RT 跨 CV 不可比）
- `secondPeptide` 保持 `True`（FAIMS 碎片谱更纯）

**验证方法**：先用 MaxQuant GUI 跑一个 raw → 检查 `msmsScans.txt` 的 `Precursor mass` 列 → 如果有大量 NaN，说明该仪器的 FAIMS 与 MaxQuant 不兼容，改用 MSFragger。

## 五、产物与日志

| 产物 | 位置 |
|---|---|
| mqpar.xml | `<result>/maxquant_auto/mqpar.xml` |
| MaxQuant 日志 | `<result>/maxquant_auto/maxquant.log` |
| msms.txt | `<第一个 raw 目录>/combined/txt/msms.txt` |
| 转换后 .spectra | `<result>/maxquant_auto/maxquant_conv.spectra` |

## 六、常见问题

| 症状 | 原因 | 解法 |
|---|---|---|
| `raw_folder/search_file_raw 未配置` | 没填或路径下无 `.raw` | 填 `search_file_raw` |
| `无 msms.txt` | MaxQuant 执行失败 | 查看 maxquant.log 尾部 |
| **msmsScans.txt 大量 NaN precursor** | **该仪器 FAIMS 与 MaxQuant 不兼容** | **改用 MSFragger** |
| first search 挂起 | Andromeda 遇到 NaN precursor 死循环 | 同上 |
| `[WARN] MaxQuant 槽位将缺席合并` | v2.0.184 E7 的醒目告警（msms.txt 缺失） | 检查上游原因 |

## 七、自定义模板

1. 在 GUI MaxQuant 中调好参数 → 另存 mqpar.xml
2. 把 fasta/raw/输出路径替换为 `{FASTA}`/`{RAWS}`/`{OUT}` 占位符
3. `maxquant_params` 键指向自定义模板
