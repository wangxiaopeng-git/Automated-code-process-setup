# 批量 PCR 引物设计工具

这是一个基于 Primer3 的 Python 命令行工具。它从 FASTA、CSV 或 TSV 批量读取模板序列，为每条模板输出若干对 PCR 引物，并把无结果样本及 Primer3 的筛选原因一并写入 CSV。

## 1. 安装

需要 Python 3.9 或更高版本。在本目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -e .
```

## 2. 准备输入

### CSV/TSV

必填列为 `id` 和 `sequence`。可选的 `target_start`、`target_length` 用来指定扩增产物必须覆盖的目标区域。`target_start` 是 **1-based**，即模板第一位为 1；两个目标字段必须同时填写或同时留空。

```csv
id,sequence,target_start,target_length
gene_A,ACGT...,120,25
gene_B,TGCA...,,
```

### FASTA

FASTA 支持多条序列，但不携带目标区域。程序用标题行第一个空白字符前的内容作为样本 ID。

```text
>gene_A optional description
ACGT...
>gene_B
TGCA...
```

模板只允许 `A/C/G/T/N`（不区分大小写）。默认不允许引物中出现 `N`；可用 `--max-ns` 调整。

## 3. 运行

### 图形界面（推荐）

macOS 可直接双击 `启动引物设计.command`。分享给新电脑时，首次启动会自动创建独立环境并安装锁定的 `primer3-py 2.4.0`，需要联网且只执行一次。也可以在终端执行：

```bash
source .venv/bin/activate
batch-primer-gui
```

启动后浏览器会打开本机界面。序列可以直接复制粘贴（支持单条裸序列或多条 FASTA），也可以上传 FASTA、CSV、TSV 文件。设置参数并点击“开始设计”，完成后可预览并下载 CSV。序列只发送到本机 `127.0.0.1`，不会上传到互联网。使用期间请保留启动时出现的终端窗口。

### 命令行

使用默认设置（引物 18–25 nt、Tm 57–63°C、产物 100–300 bp，每个样本最多 5 对）：

```bash
batch-primer-design data/examples/input.csv -o primer_results.csv
```

常见的参数调整：

```bash
batch-primer-design data/examples/input.fasta \
  -o primer_results.csv \
  --num-return 3 \
  --primer-min-tm 58 \
  --primer-opt-tm 60 \
  --primer-max-tm 62 \
  --product-min-size 150 \
  --product-max-size 400
```

查看全部参数：

```bash
batch-primer-design --help
```

## 4. 输出说明

每一行代表一对候选引物，`pair_rank=1` 为 Primer3 综合罚分最低的首选。左右引物序列都按订购时使用的 **5'→3'** 方向输出。结合区起止坐标均为模板上的闭区间 **1-based** 坐标。

- `status=ok`：找到引物对。
- `status=no_pair`：约束下没有合格引物，`message` 包含 Primer3 解释。
- `status=error`：该样本输入无效、ID 重复或 Primer3 计算异常。单条错误不会中断其他样本。

`pair_penalty` 越低越优，但程序不做全基因组特异性检索。正式订购前仍应使用 Primer-BLAST 或本地比对数据库检查物种特异性，并根据实际 PCR 体系调整盐、Mg²⁺、dNTP 和引物浓度参数。

## 5. 测试

```bash
python -m unittest discover -s tests -v
```

## 默认参数

| 参数 | 默认值 |
|---|---:|
| 候选引物对数 | 5 |
| 引物长度 | 18 / 20 / 25 nt（最小/最优/最大） |
| 引物 Tm | 57 / 60 / 63°C |
| GC 含量 | 35–65% |
| 产物长度 | 100–300 bp |
| 单价盐 | 50 mM |
| Mg²⁺ | 1.5 mM |
| dNTP | 0.6 mM |
| 引物浓度 | 50 nM |

## 参考

- Primer3 手册：https://primer3.org/manual.html
- primer3-py 文档：https://libnano.github.io/primer3-py/

