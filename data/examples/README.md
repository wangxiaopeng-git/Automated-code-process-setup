# 教学示例：先运行，再照着准备自己的输入

本目录是原 `examples/` 目录的新位置。所有序列用于演示软件，不代表真实物种实验结论。建议不要直接改这里的文件，复制到 `data/project1` 后再编辑。下面命令都从**仓库根目录**执行，且每次使用不同的结果目录。

## 各个文件是什么

| 文件 | 内容 | 可以如何替换成自己的数据 |
|---|---|---|
| `single.fasta` | 一条完整 DNA 模板 | 用文本编辑器替换标题和 DNA |
| `input.fasta` | gene_A、gene_B 两条模板 | 每条序列前写 `>唯一名称`，可追加更多记录 |
| `input.csv` | 两个样本及可选目标坐标 | 保留表头，每行一条模板 |
| `input.tsv` | 与 CSV 相同，列间用 Tab | 用表格工具另存为制表符分隔文本 |
| `input.xlsx` | Excel 第一张工作表 | 保留列名，替换每行名称和完整 DNA；不用公式 |
| `input.gb` | GenBank 格式的两条序列 | 使用真实 GenBank 文件，不靠改后缀转换格式 |
| `design-config.json` | 引物数、长度范围、Tm 和浓度等参数 | 复制后修改数字，JSON 不能加入 `#` 注释 |
| `synthetic-genome.fa` | 10 kb 人工参考，contig 名为 synthetic_contig | 真实项目换成目标物种的整个参考基因组 |
| `synthetic.gff3` | 对应参考中的两个短目标 | 必须与参考序列配套；真实项目使用同一组装注释 |
| `genes.txt` | 注释中存在的 target_plus、target_minus ID | 每行一个在你注释中确实存在的 ID |
| `targets.bed` | 上述小参考中的两个目标区间 | 按自己的参考坐标填写，注意坐标规则 |

## 不同模板格式的运行命令

这些命令相互独立，任选一个格式练习。不提供参考，所以结果里的特异性状态为 `not_run`。

```bash
# 一条 DNA，先体验最简单的设计
primer-workflow design data/examples/single.fasta -o results/example_single

# 多条 DNA，批量处理
primer-workflow design data/examples/input.fasta -o results/example_fasta

# 逗号分隔表格
primer-workflow design data/examples/input.csv -o results/example_csv

# Tab 分隔表格
primer-workflow design data/examples/input.tsv -o results/example_tsv

# Excel 第一张工作表
primer-workflow design data/examples/input.xlsx -o results/example_excel

# GenBank 格式
primer-workflow design data/examples/input.gb -o results/example_genbank
```

表格必需 `id` 和 `sequence` 列。`target_start` 和 `target_length` 可一起为空；示例第一条记录的 `80,20` 要求扩增产物覆盖模板的第 80–99 位，位置从 1 开始。它不是基因组坐标，也不是引物起点。

## 按基因 ID 输入：genes.txt

文件内容为：

```text
target_plus
target_minus
```

```bash
# 指定原始参考与配套注释，自动定位两个目标并设计
# exact 仅为了让小型演示无需额外安装 BLAST；它不允许错配检索
primer-workflow design data/examples/genes.txt --reference data/examples/synthetic-genome.fa --gff data/examples/synthetic.gff3 --specificity exact -o results/example_ids
```

真实基因组通常需要 BLAST：安装 BLAST 后删除 `--specificity exact` 即使用默认 BLAST 检查，并换一个新的结果目录。注释用于定位，不会自动把多个外显子剪接成 cDNA。

## 按区间输入：targets.bed

BED 用 Tab 分隔。示例第一行为：

```text
synthetic_contig	2000	2010	plus	0	+
```

这里显示的列间空白代表一个 Tab。含义：参考名称为 `synthetic_contig`，目标为 **0-based、右端不包含**的 `[2000,2010)`，即通常说的第 2001–2010 位；样本名为 `plus`，正链 `+`。第 5 列 `0` 是占位分值，不影响设计。文件第二条记录演示负链。

```bash
# 在目标两侧各多取默认 300 bp 模板；参考末端会自动裁剪
primer-workflow design data/examples/targets.bed --reference data/examples/synthetic-genome.fa --gff data/examples/synthetic.gff3 --specificity exact -o results/example_bed
```

BED 提供的已经是参考区间，因此注释可不提供。染色体/contig 名必须与 FASTA 标题中的 ID 完全一致，不要直接把 `chr1` 当作 `NC_000001.11`。

## JSON 参数的单位

`design-config.json` 中：`num_return` 是每条模板最多返回的候选对数；`product_min_size/product_max_size` 是产物 bp 长度；`primer_min_tm/primer_opt_tm/primer_max_tm` 是引物最低、最优、最高 Tm（°C）；`salt_monovalent_mm/salt_divalent_mm/dntp_mm` 是 mM；`primer_dna_nm` 是 nM。不要混用 mM 和 nM。

```bash
# 使用随附参数，但保存到独立结果目录
primer-workflow design data/examples/input.fasta --config data/examples/design-config.json -o results/example_config
```

每个输出目录中可先查看 `report.html` 和 `results.xlsx`；完整参数与特异性信息在 `results.json`。
