# 进阶用法与计算限制

建议先完成仓库首页的两个场景，再阅读本文。这里保留手动参考包管理方法；日常使用可直接传 `--reference genome.fa --gtf genes.gtf` 或 `--gff genes.gff3`，无需手动导入。所有命令仍在仓库根目录执行。

# 批量引物设计工作流

基于学生原始 Python/Primer3 项目的扩展版。支持多格式模板、按参考区间或基因 ID 取序列、本地参考包、NCBI 检索下载、两种本地特异性检索，以及带版本和参数记录的多格式报告。

**适用于常规基因组 PCR 候选设计。计算候选和单一预测产物不等于实验验证。** 当前不实现剪接 RT-qPCR、简并引物、探针、多重 PCR 或跨物种保守引物自动设计。

## 安装

Python 3.10+，建议 Python 3.12。Linux、macOS、Windows 可使用同一 Python CLI；BLAST+ 需额外安装并把 `blastn` 放入 PATH。

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -e .
primer-workflow --help
python -m unittest discover -s tests -v
```

Linux 安装可选 BLAST+：`sudo apt-get update && sudo apt-get install ncbi-blast+`。macOS 可使用 `brew install blast`；Windows 使用 NCBI 官方 BLAST+ 安装包。

## 1. 已有模板：直接批量设计

```bash
primer-workflow design data/examples/input.fasta -o results/templates
primer-workflow design data/examples/input.csv -o results/table --num-return 3
```

高级命令支持：FASTA/FASTA.gz、CSV/TSV（也可 gzip）、XLSX 的第一工作表、GenBank/GenBank.gz。表格必需 `id,sequence`，可选 `target_start,target_length`，目标坐标为 **1-based 闭区间的起点加长度**。不执行 Excel 公式。模板接受 A/C/G/T/N；参考中的其他 IUPAC 模糊碱基取序列时转成 N，不自动产生简并引物。

输出目录必须不存在，避免覆盖结果；重复 ID 标记为错误，其他记录继续执行。

所有 Primer3 配置可以写成 JSON，用 `--config data/examples/design-config.json` 加载；完整字段见 `DesignConfig`（`src/batch_primer_design/core.py`）。CLI 可单独覆盖候选数及产物长度。盐、Mg²⁺、dNTP 单位为 mM，引物浓度为 nM。

## 2. 本地部署参考基因组与注释

```bash
primer-workflow reference import \
  --genome /path/to/genome.fna.gz \
  --annotation /path/to/annotation.gff3.gz \
  --taxid 9606 --assembly GCF_000001405.40 \
  --annotation-release '填写实际注释发布版本' \
  -o references/human/GRCh38
primer-workflow reference verify references/human/GRCh38
```

也支持 GTF；无注释时省略 `--annotation`，仍可使用 BED 和特异性检查。检查 contig 对应、重复 contig、碱基和注释区间越界，生成标准 FASTA、注释 SQLite 和 SHA-256 清单。原文件保持不变。文件名别名不自动猜测，例如 `chr1` 和 `NC_000001.11` 不会静默混用。

每个物种/组装/注释版本使用独立目录，可离线运行；本地导入的 TaxID 和版本由使用者负责提供，NCBI 下载会与归档中的物种信息核对。验证只确认文件完整性，不替代生物学注释审查。注释索引针对 gene/transcript 定位，不是完整的 GFF3 图关系验证器。

```bash
# BED 为 0-based、右端不包含；第 4 列 ID，第 6 列可指定 +/-
primer-workflow design targets.bed --reference references/human/GRCh38 \
  --flank 300 --specificity blast -o results/bed

# txt 每行一个精确 gene/transcript ID，或无歧义的名称
primer-workflow design genes.txt --reference references/human/GRCh38 \
  --product-max-size 3000 --specificity-max-product 5000 -o results/genes
```

**ID 输入扩增的是完整基因/转录本的基因组跨度，包含内含子，不拼接 cDNA。** 长基因在默认 100–300 bp 范围下通常无结果；针对某个外显子或短区段，请提供 BED。GTF 仅有 exon/CDS 时会按 gene_id/transcript_id 聚合基因组跨度。名称有歧义必须改用精确唯一 ID；不默认挑选“最长转录本”。

仓库内提供了 10 kb **人工合成**示例（不代表真实物种）：

```bash
primer-workflow reference import --genome data/examples/synthetic-genome.fa \
  --annotation data/examples/synthetic.gff3 --taxid 32630 --assembly synthetic-v1 \
  -o references/synthetic
primer-workflow design data/examples/targets.bed --reference references/synthetic \
  --specificity exact --num-return 2 -o results/synthetic
```

## 3. NCBI 检索和下载

```bash
# 学名或 TaxID；返回前 20 条报告，保留 next_page_token 等 NCBI 原始信息
primer-workflow reference search --taxon 'Escherichia coli' --limit 20
# 明确选择版本，不能只写物种名或省略 accession 的版本后缀
primer-workflow reference download --accession GCF_000005845.2 --taxid 511145 \
  -o references/ecoli/GCF_000005845.2
```

通过官方 NCBI Datasets v2 HTTPS API，使用证书校验、超时与退避重试。下载前可用 `--max-download-gb` 限制压缩包大小（默认 10 GB；解压上限为其 4 倍），需自行预留足够磁盘。临时目录完整验证后才生成最终参考包。检查 ZIP 路径、组装与 TaxID；保留可用注释/转录本、下载记录及 SHA-256，NCBI 附带 MD5 清单时核验清单。物种没有注释时只提供基因组能力。

已有参考包可直接复用离线运行；下载不会覆盖现有包。下载中断后重新运行会重新下载，当前不实现 HTTP 断点续传。检索显示的 TaxID 可能是菌株级，与上级物种不同，应使用所选组装的 TaxID。

需要网络允许 `api.ncbi.nlm.nih.gov`，若服务发生重定向还需允许其实际官方目标；代理/TLS 错误不会被绕过。当前开发云环境阻止了该域名的 CONNECT，因此本次完成了受控 HTTP/归档集成测试，**尚未完成 NCBI 实站端到端下载验证**。API 官方参考：<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/>。

## 4. 特异性检查的含义

- `--specificity none`：显式跳过检查；不提供参考时默认采用此模式，报告 `not_run`。提供 `--reference` 时默认采用 `blast`。
- `exact`：最多 5 Mb 的参考，对两条引物及反向互补作精确全长搜索，适合小基因组和演示；**不检查错配结合**。
- `blast`：BLAST+ `blastn-short`，word size 4、关闭低复杂度屏蔽，将无 gap 的种子比对延伸到完整引物长度，重新核对错配，默认最多 2 个错配（`--max-mismatches`）。报告每个错配从 5′ 端开始的 1-based 位置。

两种模式都按同一 contig 上相向引物枚举左右与反向布局的异引物扩增子。`--specificity-max-product` 默认 2000 bp，必须覆盖设计的产物长度上限；这也限定了脱靶检索的产物范围。命中或产物达到上限、超时、工具缺失时报告 `error`，不声称完成。

`single_product_exact/blast` 表示在所记录的条件下只预测到预期产物；自备 FASTA 缺少基因组映射时使用 `single_product_unassigned_*`。`multiple_products` 表示多个潜在产物；`expected_not_found` 表示未找到映射预期产物。

BLAST 是启发式搜索，会遗漏无足够种子或带 gap 的结合位点。此处不模拟同一条引物自身配对扩增，不自动纳入未提供的细胞器、替代组装、群体变异、其他物种或污染背景；同源/多倍体拷贝必须由实验目标和参考覆盖范围决定。结果只作为筛选证据，不能作为“保证无脱靶”的判定。

## 5. 输出

每次任务输出：

| 文件 | 内容 |
|---|---|
| results.csv / results.tsv / results.xlsx | 引物、Tm/GC/长度/罚分、扩增长度、独立记录号、状态、诊断 |
| results.json | 完整候选、热力学、脱靶产物、参数、参考清单和输入校验和 |
| primers.fasta | 引物订购方向 5′→3′，以记录号和候选序号命名 |
| amplicons.fasta | 模板方向的预测扩增子 |
| primers.bed | 有参考映射的结合区；0-based 半开区间，含正负链 |
| report.html | 离线预览前 1000 行，完整结果见 JSON/表格 |
| manifest.json | 所有输出文件的 SHA-256 |

局部引物结合坐标为 1-based 闭区间；基因组字段名明确标记 0-based/右端不包含。反向模板可映射回参考链。热力学输出发卡、自/异二聚体的 Tm 和 ΔG（cal/mol），附盐条件和计算温度，属于辅助指标，并未新增统一合格阈值。

表格对公式前缀作转义，机器处理请用 JSON 保留原始 ID。`status=ok` **仅表示候选生成成功**，特异性另看 `specificity_status`。退出码：0=所有记录设计成功且无检查错误；1=有记录失败或检查错误；2=任务配置/输入等致命错误。多潜在产物会保存在报告中供筛选，退出码不将它解释为任务运行失败。

## 原有简易界面

```bash
batch-primer-gui
batch-primer-design data/examples/input.csv -o primer_results.csv
```

网页界面保留原来的粘贴序列、FASTA/CSV/TSV 和 CSV 下载功能，默认仅监听 127.0.0.1。**新增参考管理和完整报告使用 `primer-workflow` CLI**。不要将开发服务器直接暴露到公网。原始说明见 [original-readme.md](original-readme.md)，审核见 [audit.md](audit.md)。

## 验证与后续

```bash
python -m unittest discover -s tests -v
```

测试包含真实 Primer3、负链/GTF、输入错误、重复 ID、输出文件、受控 NCBI HTTP/ZIP，以及有 BLAST+ 时的真实比对。GitHub Actions 对 Python 3.10/3.12 执行测试。实际实验还应验证单一条带、产物序列、扩增效率等。

后续能力：跨外显子 RT-qPCR、VCF 避变异、跨物种保守/简并引物、探针、多重 PCR 池优化，以及大型批处理队列。它们需要独立的设计逻辑和实验验收。
