# isoPSV MaxQuant FAIMS 自动化搜库方法

> 适用版本：isoPSV v2.0.178+（MaxQuant 集成路径）／ MaxQuant 2.8.1.0+（FAIMS 原生支持）

## 一、原理

isoPSV 通过模板替换自动化 MaxQuant 搜库——用户提供 `.raw` 文件路径，软件改写 `mqpar.xml` 模板（FASTA/输出目录/文件列表/实验名/酶切类型）后调 `MaxQuantCmd.exe`，搜库完成后自动转换 `msms.txt` → `.spectra` 格式进入竞争仲裁。FAIMS 数据不需要额外配置——MaxQuant 2.x 从 Thermo `.raw` 文件头**自动读取补偿电压（CV）**并按 CV 分开搜库。

## 二、前提条件

| 项 | 要求 |
|---|---|
| MaxQuant | ≥ 2.0（本模板基于 **2.8.1.0**） |
| 原始数据 | Thermo `.raw` 格式（FAIMS 采集模式），每 CV 段完整 |
| .NET | ≥ 6.0 运行时 |
| 内存 | FAIMS 多 CV 数据内存需求 ≈ 常规 × CV 数（16GB+ 推荐） |

## 三、isoPSV 配置步骤

### config_step1.txt

```
maxquant_exe    D:/Bioinformatics/MaxQuant_2.8.1.0/bin/MaxQuantCmd.exe
maxquant_params    G:/path/to/maxquant_con.mqpar.xml
maxquant_mode    1
search_file_raw    D:/data/FAIMS_sample/
```

| 键 | 说明 |
|---|---|
| `maxquant_exe` | MaxQuantCmd.exe 完整路径（不是 MaxQuant.exe GUI 版） |
| `maxquant_params` | mqpar.xml 模板路径（可用 isoPSV 自带的 `code/engine_templates/maxquant_con.mqpar.xml`） |
| `maxquant_mode` | 1＝自动运行（0＝手动；2＝手动待办） |
| `search_file_raw` | `.raw` 文件夹或分号分隔的文件列表 |

### 模板占位符（软件自动替换，用户不用管）

| 占位符 | 替换为 |
|---|---|
| `{FASTA}` | 用户 FASTA 路径（反斜杠 Windows 格式） |
| `{OUT}` | 结果目录 |
| `{RAWS}` | raw 文件路径 XML 块 |
| `{EXP}` | 实验名 XML 块（取文件名去 .raw） |
| `{FRAC}` | 分数标签 XML 块（默认 32767＝不分馏） |

### 酶切类型（自动改写）

`config_step1` 的 `enzyme` 键自动改写 mqpar `<enzymes>` 块。映射表：

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

不支持 MaxQuant 的酶（如 none）→ 该槽位跳过＋提示。

## 四、FAIMS 专属注意事项

### 4.1 MaxQuant 自动处理（用户无需配置）

- **CV 检测**：MaxQuant 自动从 `.raw` 元数据读取 FAIMS 补偿电压列表
- **分段搜库**：每个 CV 段独立搜库（等价于把每个 CV 当独立 fraction）
- **肽段匹配**：跨 CV 段不做 Match Between Runs（默认 `matchBetweenRuns=False`，与 isoPSV 模板一致）

### 4.2 需要注意的 mqpar 参数

| 参数 | 默认值（模板） | FAIMS 建议 | 说明 |
|---|---|---|---|
| `matchBetweenRuns` | False | 保持 False | FAIMS 的 RT 跨 CV 不可比；MBR 会引入假阳性 |
| `lfqStabilizeLargeRatios` | True | 保持 | LFQ 在 FAIMS 下按 CV 独立计算 |
| `secondPeptide` | True | 保持 | FAIMS 碎片谱更好（ Selected ion 更纯），二肽鉴定收益高 |
| `minScoreModifiedPeptides` | 40 | 保持 | AndroMeda score 阈值 |

### 4.3 isoPSV 仲裁侧的 FAIMS 考虑

isoPSV 的竞争仲裁在**谱图级别**操作——每个 `msms.txt` 行对应一张谱图，FAIMS 数据与常规数据在仲裁层无区别。但需要注意：

- **成对标记验证**（Stage B）需要同一 peptide 的 light/heavy 在同一 CV 段。如果标记对在不同 CV 段洗脱（极端情况），Stage B 的 RT 窗口可能找不到配对——此时该谱标记为 `single`（无配对），不影响其他判据
- **定量**：isoPSV 不做定量（竞争仲裁软件），FAIMS 的定量问题由 MaxQuant 自己的 LFQ 处理

## 五、完整操作流程

```
1. 准备 FASTA（含标记氨基酸的同位素质量，如果做成对标记验证）
2. 把 .raw 文件放到一个文件夹（如 D:/data/FAIMS_sample/）
3. config_step1.txt 填入：
   maxquant_exe     → MaxQuantCmd.exe 路径
   maxquant_params  → mqpar 模板路径（可用自带模板）
   maxquant_mode    1
   search_file_raw  → .raw 文件夹路径
   enzyme           trypsin（按实际酶切类型）
   fasta            → FASTA 路径
   ...
4. 网页版点"运行"或 CLI 直接跑 user_run.R
5. isoPSV 自动：
   → 生成 mqpar.xml（替换占位符＋酶切类型）
   → 调 MaxQuantCmd.exe（后台异步，日志落盘）
   → 等待完成（进度心跳在主日志可见）
   → 从 combined/txt/msms.txt 转换为 .spectra 格式
   → 自动登记为竞争槽位（策略字母按引擎字母序分配）
   → 进入五策略合并 → 竞争仲裁
```

## 六、产物与日志

| 产物 | 位置 |
|---|---|
| mqpar.xml（改写后） | `<result>/maxquant_auto/mqpar.xml` |
| MaxQuant 日志 | `<result>/maxquant_auto/maxquant.log` |
| msms.txt（原始） | `<第一个 raw 所在目录>/combined/txt/msms.txt` |
| 转换后 .spectra | `<result>/maxquant_auto/maxquant_conv.spectra` |
| 参数侧车 | `<result>/maxquant_auto/maxquant_auto.run_params` |
| msms.txt→spectra 转换日志 | 同 maxquant.log（追加） |

## 七、常见问题

| 症状 | 原因 | 解法 |
|---|---|---|
| `[maxquant 自动] raw_folder/search_file_raw 未配置或无 .raw 文件，跳过` | 没填 `search_file_raw` 或路径下无 `.raw` | 填入 raw 文件夹路径 |
| `[maxquant 自动] 无 msms.txt（exit N）` | MaxQuant 执行失败 | 查看 `maxquant.log` 尾部；常见＝.NET 版本、raw 文件损坏 |
| MaxQuant 找不到 FAIMS 参数 | MaxQuant 版本 < 2.0 | 升级到 2.x+ |
| msms.txt 为空 | FAIMS raw 文件无鉴定 | 检查 raw 文件是否正确采集；MaxQuant FDR 过严 |
| `[maxquant 自动][WARN] 适配器转换失败` | msms.txt 列名/格式异常 | 查看 maxquant.log 尾部；确保 MaxQuant 版本与模板匹配 |
| 内存不足 | FAIMS 多 CV × 大 FASTA | 减小 FASTA 或升级内存 |

## 八、自定义模板

自带模板（`maxquant_con.mqpar.xml`）覆盖常规场景。如需自定义（如加修饰、改质量容差）：

1. 复制 `code/engine_templates/maxquant_con.mqpar.xml` 为自定义模板
2. 在 GUI 版 MaxQuant 中调好参数 → 另存 mqpar.xml → 用文本编辑器把 fasta/raw/输出路径替换为 `{FASTA}`/`{RAWS}`/`{OUT}` 占位符
3. `config_step1` 的 `maxquant_params` 指向自定义模板

**注意**：`{RAWS}`/`{EXP}`/`{FRAC}` 是块级占位符（整块 `<string>...</string>` 行替换），不要只替换路径文本。
