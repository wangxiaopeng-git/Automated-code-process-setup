# 参数、参考来源与结果解释

本文与仓库首页构成完整使用手册。安装、目录规划和逐步操作见 [首页](../README.md)。所有命令在仓库根目录执行。

## 1. 选择参考方式

| 命令参数 | 参考来源 | 默认检查 |
|---|---|---|
| 不提供参考 | 自备模板 | none |
| `--reference genome.fa` | 本地完整基因组，可选 `--gtf` 或 `--gff` | blast |
| `--reference references/species` | 已导入的参考包 | blast |
| `--assembly GCF_000005845.2` | NCBI 元数据、所需区域、在线数据库 | remote |

`--reference` 和 `--assembly` 互斥。`--specificity none` 是显式跳过检查；不允许把本地 `blast/exact` 与在线组装混用，也不允许用 `remote` 假装检查本地参考。`--gtf/--gff` 只用于原始本地 FASTA，互斥且不能另外附加在已导入的参考包上。

## 2. 输入与坐标

- FASTA、CSV/TSV、XLSX、GenBank：模板序列输入，不能仅凭文件名称推断其物种或坐标。
- BED：0-based 半开区间，列为 contig、start、end、可选 ID、score、strand；对目标加 `--flank`（默认 300 bp）后取模板。不能把不同版本的坐标混用。
- TXT：本地模式每行一个注释 gene/transcript ID；在线模式每行一个数字 NCBI Gene ID。Gene 摘要的坐标必须与所选组装的版本化序列匹配，否则要求改用 BED。
- 表格 `target_start/target_length`：模板上的 1-based 起点和长度，要求产物覆盖目标，不是引物固定起点。
- 基因/转录本 ID 取其完整基因组跨度，包含内含子，不自动剪接 cDNA。长跨度需要合理的产物长度，或改为特定外显子/区间。

局部引物结合区输出为 1-based 闭区间；基因组字段以 `_genomic_start_0based`、`_genomic_end_exclusive` 标明半开区间。负链模板映射回原参考后保留正确的链方向。

## 3. Primer3 与检索参数

```bash
# 读取 JSON 参数，再用命令行覆盖候选数和产物长度
primer-workflow design data/project1/targets.fasta --config data/examples/design-config.json --num-return 3 --product-min-size 150 --product-max-size 400 -o results/project1_parameters
```

JSON 字段对应 `DesignConfig`。长度单位 nt/bp，Tm 为 °C，单价盐、Mg²⁺ 和 dNTP 为 mM，引物浓度为 nM。最小/最优/最大值必须有序，浓度必须有限且在允许范围内。JSON 不支持 `#` 注释。

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `--num-return` | 5 | 每条模板最多返回的候选对数，最大 100 |
| `--flank` | 300 | 参考区间两侧取序列的长度 |
| `--product-min-size / --product-max-size` | 100 / 300 | 设计产物范围 |
| `--specificity-max-product` | 2000 | 潜在脱靶产物的最大长度，必须覆盖设计上限 |
| `--max-mismatches` | 2 | BLAST 全长无 gap 结合位点的最多错配数，0–5 |
| `--max-region-bp` | 200000 | 在线目标加 flanks 的大小上限，最大 1 Mb |
| `--remote-max-queries` | 40 | 去重后在线引物查询上限，最大 100 |
| `--remote-hit-limit` | 100 | 在线每条查询的返回命中上限，最大 500 |
| `--remote-timeout` | 900 | 在线任务的轮询等待预算，60–3600 秒；请求也有独立网络超时 |
| `--cache-dir` | 项目内隐藏目录 | 本地默认 `.primer-cache`，在线默认 `.ncbi-cache` |

较大的脱靶产物上限会增加组合与资源开销；小的上限会漏掉范围以外的产物。增加候选数不保证得到更好的引物，在线方式还受查询限额影响。

## 4. 本地参考包

日常可以直接传 FASTA 和注释，由程序缓存。希望记录明确的物种/组装/注释版本时，可手动导入：

```bash
# 建立经过验证的本地参考包；请把路径、TaxID 和版本换成自己的数据
primer-workflow reference import --genome data/project1/genome.fna.gz --annotation data/project1/genes.gff3.gz --taxid 9606 --assembly GCF_000001405.40 --annotation-release release_name -o data/project1/reference

# 校验参考包文件是否改变
primer-workflow reference verify data/project1/reference

# 使用参考包，不需要重复传注释
primer-workflow design data/project1/targets.bed --reference data/project1/reference -o results/project1_bundle
```

手动导入的 TaxID/版本由使用者负责准确填写；自动准备原始文件时未知 TaxID 记录为 null，使用内容哈希标记本地参考，不会猜测物种身份。检查覆盖 contig、重复序列 ID、碱基和注释越界；注释索引用于定位，不是完整的 GFF3 图关系审查。

## 5. NCBI 检索、区域获取与完整下载

在线功能为实验性接口，尚未完成公开服务端到端验证；请先核对已知目标的版本和坐标。

### 只查询元数据

```bash
# 确保项目目录已经建立，按物种查找候选组装
primer-workflow reference search --taxon 'Escherichia coli' --limit 20 > data/project3/assemblies.json

# 获取指定版本的组装报告和全部分页序列清单，不下载完整基因组
primer-workflow reference inspect --assembly GCF_000005845.2 --cache-dir data/project3/.ncbi-cache > data/project3/assembly.json
```

检索的 `limit` 为 1–1000，保留 NCBI 原始分页提示，不把第一批结果称为全部组装。`inspect` 会检查元数据的精确版本及序列清单分页完整性；返回结果含原始组装报告、序列数量和清单摘要，缓存保留原始响应。

### 在线轻量获取

```bash
# 输入是目标版本上的 BED：只获取目标附近序列和可用区域 features
primer-workflow design data/examples/ncbi-targets.bed --assembly GCF_000005845.2 --specificity none -o results/ncbi_regions

# 输入是自备模板：不取整个基因组，提交候选做在线筛查
primer-workflow design data/project3/targets.fasta --assembly GCF_000005845.2 --num-return 2 -o results/ncbi_screen
```

EFetch 使用版本化核酸 accession、`seq_start/seq_stop` 和 `gbwithparts` 获取局部 GenBank。区域大小、返回 accession/version 和长度会核对；GCA 组装优先使用 GenBank accession，GCF 优先 RefSeq。区域注释可能为空或局部截断，也不一定等同于最新版 Datasets GFF 发布；不自动扩展成完整基因图谱。

GET 请求使用受限次数的重试与退避；每次实际尝试都经过速率限制。超过一分钟的服务端等待要求会报“稍后重试”，而不是无限阻塞。缓存校验和不匹配时停止使用该响应；可在程序停止后移除相应缓存文件再重试。默认响应缓存有效期 7 天，BLAST 已完成结果/等待任务缓存有效期 1 天；时间不等于数据发布版本。

### 确实需要完整基因组时

```bash
# 下载指定组装的整个基因组和可用注释，导入本地参考包；可能占用大量空间
primer-workflow reference download --accession GCF_000005845.2 --taxid 511145 -o data/project3/reference

# 后续离线设计，使用本地 BLAST+
primer-workflow design data/project3/targets.fasta --reference data/project3/reference -o results/project3_local
```

完整下载与轻量路径是不同选择：完整下载默认压缩包上限 10 GB（`--max-download-gb`），解压上限是其 4 倍。检查路径、组装和 TaxID，记录 SHA-256；官方提供 MD5 清单时核验并记录覆盖情况。不会覆盖已有参考包；不实现压缩包 HTTP 断点续传。

## 6. 特异性证据的边界

**本地 exact**：最多 5 Mb 参考，搜索完整引物及反向互补的完全匹配，不查错配。

**本地 blast**：BLAST+ `blastn-short`、word size 4、关闭低复杂度屏蔽。无 gap 的局部种子比对延伸至完整引物后，重新计算错配。BLAST 仍是启发式，不保证发现所有可能结合。

**在线 remote**：NCBI URLAPI `blastn`、nt 数据库、物种 TaxID 过滤、word size 7。批量去重提交引物，按 RID 查询，记录数据库和时间。返回命中只在版本化序列 ID 与所选组装匹配时使用；同一个合并命中列出多个组装 contig 时逐个核对。对局部 HSP 获取完整短结合区，计算从 5′ 端开始的 1-based 错配位置。即使有预期产物，仍不能证明整个组装已被在线数据库覆盖；这种结果固定使用 `remote_limited`。

所有方式按同一 contig 上相向的 L/R 引物和反向布局枚举潜在产物，受最大产物长度限制。不模拟同一条引物自身配对扩增、带 gap 结合、群体变异、未提供的细胞器/替代组装或污染背景。物种、多倍体和同源拷贝的实验目标需要使用者判断。

| 状态 | 解读 |
|---|---|
| not_run | 没检查 |
| single_product_exact / blast | 本地检索条件内单一预期产物 |
| single_product_unassigned_* | 没有预期坐标，只找到一个潜在产物；先核对位置 |
| single_product_remote_limited | 在线证据中单一预期产物，全组装覆盖未验证 |
| multiple_products / multiple_products_remote | 多潜在产物，存在脱靶风险 |
| expected_not_found* / no_product_found* | 预期位置或产物未被确认，不代表“无脱靶” |
| incomplete_remote | 命中达到上限、查询缺失或无法核对序列版本 |
| error | 检查失败，见 specificity_error |

候选保持 Primer3 的 `pair_rank`，不会用远程分数冒充 Primer3 罚分。`screening_priority` 给出复核类别，不是实验合格结论。在线证据有限的候选应优先核对坐标；必要时用本地完整参考复核，再验证条带、产物序列等。

## 7. 输出文件

| 文件 | 内容 |
|---|---|
| results.csv / tsv / xlsx | 候选、引物方向、Tm/GC、长度、罚分、状态和诊断 |
| results.json | 完整参数、热力学、潜在产物、参考来源及输入哈希 |
| primers.fasta | 订购方向 5′→3′ 的引物 |
| amplicons.fasta | 模板方向的预测扩增子 |
| primers.bed | 有映射时的参考结合区，0-based 半开 |
| report.html | 离线预览前 1000 行 |
| online_reference.json | 在线组装信息、区域 features、请求校验记录和范围限制 |
| online_templates.fasta | BED/Gene ID 在线取到的局部模板，可供复用 |
| region_features.tsv | 区域注释，compound feature 拆成多行并保留 part_index 和 qualifiers；不是完整 GTF/GFF |
| manifest.json | 输出文件 SHA-256 |

`online_templates.fasta` 与 `region_features.tsv` 只在确实在线获取了目标区域时生成；没有区域 feature 时注释表只有表头，JSON 有 annotation_status。FASTA 自备模板不被说成在线下载的序列。

发卡、自二聚体和异二聚体输出包含 Tm、ΔG（cal/mol）以及计算条件，没有统一合格阈值。表格对公式前缀转义，原始字段请用 JSON。退出码 0 表示设计成功且无检查运行错误，1 表示有记录或检查失败/不完整，2 表示致命配置/输入问题；“在线有限证据”和“发现多产物”本身并不等于程序运行失败。

## 8. 服务与依据

- NCBI Datasets REST API：<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/>
- 官方接口定义：<https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml>
- NCBI E-utilities 参数：<https://www.ncbi.nlm.nih.gov/books/NBK25499/>
- BLAST Common URL API：<https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html>

本程序只实现上述 NCBI 在线来源。Ensembl、UCSC 等提供的本地参考可导入，但不能混用版本名或坐标。公共接口可能变化；服务报错时保留错误与 RID，不绕过 TLS 或网络策略，不将错误解释为检查成功。
