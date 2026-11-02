# 批量引物设计：从准备文件到查看结果

这是一个使用 Primer3 设计 **常规 PCR 引物候选**的工具。支持单序列与批量输入、本地基因组检查，以及按 NCBI 组装版本获取目标区域和在线 BLAST 筛查。

从安装到结果，按以下顺序使用：**安装 → 建立项目目录 → 准备输入 → 选择设计方式 → 查看候选和特异性证据**。

| 方式 | 适合的情况 | 需要下载整个基因组吗 | 检查范围 |
|---|---|---|---|
| 场景一：只有模板 | 快速生成候选，不做脱靶检索 | 不需要 | 未检查特异性 |
| 场景二：本地参考 | 重复分析、大批量、需要可控参考范围 | 需要本地参考 | 提供的完整参考序列 |
| 场景三：在线轻量 | 少量目标、磁盘有限、按组装版本取序列 | 不需要 | NCBI 在线数据库返回命中，按组装序列版本筛选；覆盖不保证完整 |

在线方式不意味着必然比本地更准确。它减少本地存储，但受公共服务速度、数据库收录、命中数上限和启发式搜索影响；可靠性需要看检查证据，而不是只看找到几对引物。

推荐第一次使用时，先完整运行下面的“场景一”，再做“场景二”。**示例序列只用于练习软件，不是已经实验验证的引物或真实物种研究结论。**

## 一、先看懂文件放在哪里

以后每个实验项目建一个文件夹。`project1` 只是示范名字，可以换成 `rice_geneA`、`mouse_pcr` 等，建议使用英文字母、数字和下划线。

```text
Batch-primer-design-workflow/       ← 下载下来的整个代码文件夹（仓库根目录）
├── data/
│   ├── examples/                  ← 软件附带的示例，不需要自己准备
│   │   ├── single.fasta           ← 只有一条序列，适合第一次练习
│   │   ├── input.fasta            ← 两条序列，演示批量设计
│   │   ├── input.csv / input.tsv  ← 相同用途的表格示例
│   │   ├── input.xlsx / input.gb  ← Excel / GenBank 示例
│   │   ├── synthetic-genome.fa    ← 人工合成的小型参考基因组
│   │   ├── synthetic.gff3         ← 与这个小基因组配套的注释
│   │   ├── genes.txt              ← 上述注释中的基因 ID 示例
│   │   ├── ncbi-targets.bed       ← NCBI 在线组装的区间示例
│   │   ├── ncbi-genes.txt         ← 在线 NCBI Gene ID 示例
│   │   └── targets.bed            ← 上述基因组上的目标区间示例
│   └── project1/                  ← 你自己新建，放输入序列及配套参考文件
│       ├── targets.fasta
│       ├── genome.fa              ← 场景二才需要
│       └── genes.gtf              ← 场景二可选，也可以是 genes.gff3
├── results/                       ← 程序运行后自动生成
│   └── project1/                  ← 这个项目的结果
│       ├── results.xlsx
│       ├── report.html
│       └── ...
└── README.md                      ← 这份使用说明
```

`data/project1` 和 `results/project1` **不是同一个文件夹**：前者放你准备的输入，后者放程序计算的输出。不要把输入文件放进结果目录。默认会忽略用户项目数据、缓存和结果，避免用 Git 时意外提交大基因组或私人数据；`data/examples` 中的示例仍随代码分发。

## 二、安装软件（第一次使用时做一次）

需要 Python **3.10 或更高版本**，推荐 3.12，并已安装 Git。`#` 开头的是注释，解释下一行命令的用途；不用输入终端前面的 `$` 或 `>>>`。以下 Bash 命令适用于 macOS/Linux；Windows 的安装命令单独列在后面。

### macOS / Linux

```bash
# 检查 Python 版本；若提示不存在，请先安装 Python 3.12
python3 --version

# 下载本项目的可运行代码分支
git clone --branch feature/primer-workflow https://github.com/wangxiaopeng-git/Batch-primer-design-workflow.git

# 进入仓库根目录；后续所有相对路径都以这里为起点
cd Batch-primer-design-workflow

# 创建独立的 Python 环境，避免与其他项目的软件包混在一起
python3 -m venv .venv

# 激活这个环境；以后重新打开终端时，也要先 cd 到仓库并执行这一行
source .venv/bin/activate

# 安装本项目和所需依赖；点号 . 表示当前代码文件夹，首次安装需要联网
python -m pip install -e .

# 确认命令安装成功；会显示可用命令和参数说明
primer-workflow --help
```

### Windows PowerShell

先安装 Python 3.12 和 Git，然后在 PowerShell 中执行：

```powershell
# 下载代码
git clone --branch feature/primer-workflow https://github.com/wangxiaopeng-git/Batch-primer-design-workflow.git

# 进入仓库根目录
cd Batch-primer-design-workflow

# 使用 Python 3.12 创建独立环境
py -3.12 -m venv .venv

# 激活环境
.\.venv\Scripts\Activate.ps1

# 安装程序和依赖
python -m pip install -e .

# 显示帮助，确认安装完成
primer-workflow --help
```

如果 PowerShell 不允许激活脚本，可以直接使用 `.\.venv\Scripts\python.exe -m pip install -e .` 安装；后续把 `primer-workflow` 替换成 `.\.venv\Scripts\primer-workflow.exe` 即可。下面设计命令都写成一行，Windows 也可以直接使用。

## 三、场景一：只提供模板序列，直接设计引物

**你需要准备什么？** 一条或多条 DNA 序列。不需要基因组、不需要注释，也不需要安装 BLAST。程序会计算引物，但不会检查全基因组脱靶。

### 第一步：先用随附的单序列示例跑通

```bash
# 读取已经提供的单序列示例，把结果保存到 results/demo_single
primer-workflow design data/examples/single.fasta -o results/demo_single
```

运行结束后，用文件管理器进入 `results/demo_single`，双击 `report.html` 查看概览，或用 Excel 打开 `results.xlsx`。如果成功，会显示 `successful_records: 1`；默认最多返回 5 对候选。同一个目录不允许覆盖，重复练习时把输出改成 `results/demo_single_run2`。

### 第二步：复制示例，建立自己的项目

macOS/Linux：

```bash
# 只建立输入目录；不要手动建立 results/project1，程序会自动建立
mkdir -p data/project1

# 复制包含两条模板的示例；之后可编辑这个副本，不要修改原始示例
cp data/examples/input.fasta data/project1/targets.fasta

# 用你的项目文件设计引物，并明确指定项目结果目录
primer-workflow design data/project1/targets.fasta -o results/project1
```

Windows PowerShell 的前两行改为：

```powershell
# 新建项目输入目录；已存在时可继续使用
New-Item -ItemType Directory -Force data/project1

# 复制示例，准备自己的输入文件副本
Copy-Item data/examples/input.fasta data/project1/targets.fasta

# 开始设计，与 macOS/Linux 使用相同的参数
primer-workflow design data/project1/targets.fasta -o results/project1
```

不写 `-o` 时，程序也会按输入文件所在文件夹的名字选择输出：例如 `data/project1/targets.fasta` 对应 `results/project1`。这与上面的显式写法等价；**如果这个结果目录已经存在，先换一个新的输出目录再运行**，例如 `-o results/project1_run2`，不需要删除旧结果。

### 第三步：把示例换成自己的序列

用纯文本编辑器打开 `data/project1/targets.fasta`，不要用 Word 保存。文件的结构是“标题行 + DNA 序列”；下面这条是可直接识别的完整示例：

```fasta
>my_gene_1
ACGTTGCAAGTCGATCGTACGATCGGCTAGCTAGGCTAGCTACGATCGATCGTACGATCGATGCTAGCTAGGATCGATCGTACGTAGCTAGCTGACTGATCGATCGTAGCTAGCTACGTAGCTGATCGTAGCTAGCTACGATCGTAGCTGACTGATCGATCGTAGCTAGCTGACTGATCGATCGTAGCTAGCTAG
```

- `>` 必须在每条标题行开头，后面是唯一的样本名称，例如 `my_gene_1`。名称尽量不带空格，否则只使用第一个空格前的内容。
- 下一行替换成你实际的 DNA 模板序列；较长序列可以分成多行。支持 A/C/G/T/N，不区分大小写。
- 要做第二个基因，就继续添加 `>my_gene_2` 和它的序列。多条序列会逐条设计，名称不能重复。
- **基因名称不是 DNA 序列**，例如只输入 `BRCA1` 无法在这个场景中设计。若只有基因 ID，请看场景二的基因 ID 用法。
- 不要把说明中的省略号 `...` 写进 DNA 序列。模板需要有足够的长度放置两条引物和扩增区；默认产物长度为 100–300 bp。

换成自己的序列后，用相同命令运行，并选择一个新结果目录：

```bash
# 相同输入文件可以反复编辑、重新分析，但每次使用不同的结果目录
primer-workflow design data/project1/targets.fasta -o results/project1_run2
```

### 如果你更习惯 Excel 或表格

可以复制 `data/examples/input.xlsx` 作为模板，用 Excel 编辑：第一行列名必须保留 `id` 和 `sequence`，每一行放一个样本。完整 DNA 放在 `sequence` 单元格，不要插入省略号或公式。`target_start` 和 `target_length` 是可选列，可一起留空。

```bash
# 以 Excel 格式为例；请先把编辑后的文件放到指定位置
primer-workflow design data/project1/targets.xlsx -o results/project1_excel

# 不想准备文件时，也可以先直接跑提供的 Excel 示例
primer-workflow design data/examples/input.xlsx -o results/demo_excel
```

其他格式只要更换输入路径即可：FASTA（`.fa/.fna/.fasta`，可 gzip）、CSV、TSV（可 gzip）、Excel `.xlsx`、GenBank `.gb/.gbk`。每种格式的示例和字段解释见 [示例说明](data/examples/README.md)。

`target_start=80,target_length=20` 表示产物必须覆盖模板第 80–99 位，坐标从 1 开始；**它不是指定引物必须从第 80 位开始**。两个字段必须一起填写或一起留空。

## 四、场景二：提供本地基因组，提高特异性筛选能力

### 第一步：准备目标物种的文件

把文件下载、解压后放在自己的项目目录，例如：

```text
data/project1/
├── targets.fasta     ← 希望设计引物的模板序列
├── genome.fa         ← 目标物种的参考基因组（不是某一个基因的序列）
└── genes.gtf         ← 与这个基因组同一组装版本的注释；也可以是 genes.gff3
```

这些文件可以有其他名字，只要命令中的路径与真实文件名相同即可。支持 `genome.fna.gz`、`genes.gtf.gz`、`genes.gff3.gz`，不要求手动解压 gzip 文件。

**去哪里下载？** 在 [NCBI Datasets Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) 搜索目标物种学名，选择一个明确的组装版本，下载 Genome FASTA 和可用的 GFF3 注释。GTF 也可以使用其他权威来源的配套文件，但必须与基因组属于**同一个组装版本**，染色体名称也要一致。不要把 GRCh37 的基因组与 GRCh38 的注释混用，也不要自行把 GFF3 改名成 GTF 来转换格式。

**注释是可选的。** 检查引物在全基因组中是否有其他潜在扩增位置，只需要参考基因组。注释的作用是按基因/转录本 ID 找到目标区域；如果输入已经是 FASTA 模板，就不强制提供注释。没有注释的物种也能使用这个场景。

### 第二步：安装用于特异性检查的 BLAST+

这是独立工具，不包含在 Python 的普通安装中。安装后请重新打开终端、激活 `.venv`，并确认 `blastn -version` 能显示版本号。

```bash
# Ubuntu/Debian：刷新系统软件列表
sudo apt-get update

# Ubuntu/Debian：安装 NCBI BLAST+（需要管理员权限）
sudo apt-get install ncbi-blast+

# 检查是否安装成功
blastn -version
```

macOS 已安装 Homebrew 时，可执行 `brew install blast`。Windows 请从 [NCBI 官方 BLAST+ 安装目录](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 下载对应安装包，并将含 `blastn.exe` 的文件夹加入 PATH；重新打开终端后执行 `blastn -version`。

### 第三步：直接传入本地文件，无需手动建立参考包

以下是几种**替代写法**，根据自己有的文件选择一条运行，不需要依次全部执行：

```bash
# 只有基因组也可以：设计候选后，自动用 BLAST 检查潜在脱靶产物
primer-workflow design data/project1/targets.fasta --reference data/project1/genome.fa -o results/project1_reference

# 如果有 GTF：自动检查基因组与注释是否匹配，并建立本地缓存
primer-workflow design data/project1/targets.fasta --reference data/project1/genome.fa --gtf data/project1/genes.gtf -o results/project1_gtf

# 如果注释是 GFF/GFF3，改用 --gff；不要同时使用 --gtf 和 --gff
primer-workflow design data/project1/targets.fasta --reference data/project1/genome.fa --gff data/project1/genes.gff3 -o results/project1_gff
```

**默认行为很简单：**

| 你提供的参数 | 程序做什么 |
|---|---|
| 没有 `--reference` | 只按输入模板设计，特异性状态为 `not_run` |
| 有 `--reference genome.fa` | 自动准备参考，并默认启用 BLAST 特异性检查 |
| 再加 `--gtf` 或 `--gff` | 同时导入配套注释，可支持基因 ID 输入 |
| 显式写 `--specificity none` | 按你的选择跳过检查，即使提供了参考也不会声称已验证 |

程序会在 `data/project1/.primer-cache/` 自动生成可复用参考；这是内部缓存，不需要手工编辑。基因组和注释的**内容**变化后会准备新的缓存，不能靠改文件名冒充新版本。需要节省空间时可以在程序停止后删除自己项目的 `.primer-cache`，下次会重建；原始输入文件不会被删除或改写。大型基因组首次处理会更慢，并需要额外磁盘存放标准化参考。

如果 BLAST 未安装，会给出错误提示，不会自动退回“不检查”。如果你明确只想试运行一个不超过 5 Mb 的小参考，可以写 `--specificity exact`，但该模式只查完全匹配，**不查错配结合，不能代替正式 BLAST 筛选**。

### 不用下载真实基因组，也能先练习场景二

仓库已提供一个 **10 kb 人工小基因组**、配套 GFF3 和基因 ID，便于练习。macOS/Linux：

```bash
# 为第二个练习单独建立输入目录，不影响场景一的 project1
mkdir -p data/project2

# 复制配套的示例基因组、注释和目标 ID；这三个文件必须配套使用
cp data/examples/synthetic-genome.fa data/project2/genome.fa
cp data/examples/synthetic.gff3 data/project2/genes.gff3
cp data/examples/genes.txt data/project2/genes.txt

# 自动读取基因 ID 并定位目标；为方便第一次运行，显式使用小参考的精确匹配模式
primer-workflow design data/project2/genes.txt --reference data/project2/genome.fa --gff data/project2/genes.gff3 --specificity exact -o results/project2
```

Windows PowerShell 把 `mkdir -p` 换成 `New-Item -ItemType Directory -Force data/project2`，把 `cp` 换成 `Copy-Item`；设计命令保持相同。安装好 BLAST 后，可用上面同一组文件删除 `--specificity exact`，改用新的输出目录 `-o results/project2_blast`，测试默认 BLAST 流程。

`genes.txt` 中每行一个注释中的 ID，例如示例里的 `target_plus` 和 `target_minus`。换成真实项目时，必须换成你自己注释中存在的 ID。

**目前按基因/转录本 ID 设计时，要求覆盖其完整基因组跨度，包含内含子，不会自动拼接 cDNA。** 很长的基因在默认 100–300 bp 产物范围下可能没有候选；若只想设计某一个外显子或短区段，请用 [BED 区间示例](data/examples/README.md#按区间输入targetsbed)。这不是 RT-qPCR 的跨外显子设计功能。

## 五、场景三：按 NCBI 组装版本在线获取和筛查（实验性）

在线接口为实验性功能，尚未完成 NCBI 实站端到端验证。使用时应先对已知目标核对版本、序列和坐标。

适合少量目标。**在线模式不需要安装本地 BLAST，也不下载整条染色体或整个基因组。** 首先获取组装和染色体/contig 元数据，再按输入区间获取局部序列及可用注释；需要检查时，把去重后的候选引物批量提交给 NCBI BLAST。

**在线 BLAST 会将引物序列发送到 NCBI。** 如果序列不能提交外部服务，请使用本地参考方式。公共数据库不一定收录所有物种和所有组装版本，程序也不保证适用于每一个物种。

### 第一步：查找并确认组装版本

```bash
# 建立输入项目目录；Windows PowerShell 用 New-Item -ItemType Directory -Force data/project3
mkdir -p data/project3

# 按物种学名检索组装，将返回报告保存在项目目录；这里只获取元数据
primer-workflow reference search --taxon 'Escherichia coli' --limit 20 > data/project3/assemblies.json

# 确认选择的精确版本，以及对应的染色体 accession 和长度；仍不下载整个基因组
primer-workflow reference inspect --assembly GCF_000005845.2 --cache-dir data/project3/.ncbi-cache > data/project3/assembly.json
```

`GCF_000005845.2` 是组装号；`NC_000913.3` 是其中一条核酸序列的版本号，二者不是同一种 ID。选择组装时应核对物种/菌株、组装名称、注释来源和版本。不要只写 `GCF_000005845` 或一个学名来代替精确版本。

### 第二步：选择自己的输入方式

**A. 已有 DNA 模板：直接在线筛查。** 把自己的 FASTA 放在 `data/project3/targets.fasta`，运行：

```bash
# 使用指定组装的序列清单核对在线 BLAST 返回的命中；最多生成 2 对/模板
primer-workflow design data/project3/targets.fasta --assembly GCF_000005845.2 --num-return 2 -o results/project3_online
```

这不会自动证明你的 FASTA 来自该组装。若没有目标基因组坐标，即使找到一个产物，也会显示 `single_product_unassigned_remote_limited`，要求核对位置。FASTA 本身不会提供完整基因注释；需要区域注释时，使用下面的 BED 输入。

**B. 已知基因组位置：仅获取目标附近序列和区域注释。** 先用提供的真实参考区间示例练习：

```bash
# 复制与 GCF_000005845.2 对应的区间示例；真实项目请修改成自己的组装坐标
cp data/examples/ncbi-targets.bed data/project3/targets.bed

# 仅获取目标两侧各 300 bp 及可用区域注释，设计引物，默认执行在线 BLAST
primer-workflow design data/project3/targets.bed --assembly GCF_000005845.2 --num-return 2 -o results/project3_regions

# 如果只想取序列并设计，不提交 BLAST，则明确关闭检查，并另存一次结果
primer-workflow design data/project3/targets.bed --assembly GCF_000005845.2 --specificity none -o results/project3_fetch_only
```

Windows PowerShell 把 `cp` 换成 `Copy-Item`。示例 BED 第一行是 `NC_000913.3` 的 `[363300,363310)` 区间，即常用的 1-based 第 363301–363310 位。真实项目必须同时替换**组装号、序列 ID 和坐标**，不能只替换物种名。文件为 Tab 分隔。

**C. 只有 NCBI Gene ID：按版本匹配坐标。** `.txt` 每行一个数字 Gene ID；基因符号和 Ensembl ID 不会被默认猜测成 NCBI ID。

```bash
# 复制数字 Gene ID 示例；真实项目换成自己基因的 NCBI Gene ID
cp data/examples/ncbi-genes.txt data/project3/genes.txt

# 要覆盖完整基因跨度，应按真实基因长度设置产物上限；此处仅演示提取与设计
primer-workflow design data/project3/genes.txt --assembly GCF_000005845.2 --specificity none --num-return 1 --product-max-size 6000 --specificity-max-product 8000 -o results/project3_gene
```

Gene ID 示例 `945006` 需要在线元数据能够将该基因的坐标明确映射到所选组装中的序列版本。旧组装可能不在当前 Gene 摘要坐标中，程序会报“没有唯一匹配坐标”，不会偷偷改用新组装。这种情况请提供已确认的版本对应 BED。长基因、重复区等也可能因设计约束返回 `no_pair`，示例不承诺每次公共数据更新后都能产生候选。

### 第三步：了解在线结果和资源限制

- 每个目标加两侧序列默认最多 **200 kb**，用 `--max-region-bp` 可调整，硬上限 1 Mb；超过时改用短 BED 或本地方式。注释来自该区域 EFetch 返回的 feature，不是完整 GTF/GFF 下载，也不保证与最新 Datasets 注释发布相同。没有 feature 时仍可获取序列。
- 缓存在项目的 `.ncbi-cache/` 中，仅存元数据、目标区域、短引物结合区域和 BLAST 任务/结果；响应有校验和和有效期。它不替代本地全基因组数据库。不要把缓存提交到 Git。
- 引物去重后默认最多 **40 条查询/任务**（不是 40 对），可以用 `--remote-max-queries` 调整，硬上限 100。在线方式面向小批量；目标很多时用本地参考。
- 每条查询默认最多返回 100 个命中，用 `--remote-hit-limit` 调整，硬上限 500；达到上限或缺少查询结果时显示 `incomplete_remote`，不能当作没有脱靶。
- 公共服务请求会自动限速。BLAST 可能等待数分钟；`--remote-timeout` 默认 900 秒，等待超时保留任务 RID。对于仍在等待的任务，使用**相同输入、参数和缓存**重新执行，可继续查询仍有效的 RID。若已生成错误报告，先换一个新的输出目录；明确失败/未知 RID 任务需检查记录并使用新缓存目录重新提交。超过一天的任务/结果缓存过期后会重新提交。
- 可用 `--ncbi-email you@example.org` 提供你自己的联系邮箱；示例地址必须换成自己的真实地址，不要求密钥。
- 访问需要允许 `api.ncbi.nlm.nih.gov`、`eutils.ncbi.nlm.nih.gov`、`blast.ncbi.nlm.nih.gov`。网络失败会记录错误，不会改成“检查通过”。

在线筛查使用 **NCBI nt + 物种过滤**，再把返回的版本化核酸序列与指定组装的序列清单取交集。**这不是把整套在线数据库严格限制在某个组装，也不能证明该组装所有序列都已被搜索。** 新组装未收录、nt 合并条目、短引物种子不足、命中上限等都影响覆盖。单一产物状态因此固定带 `remote_limited`，不会当作本地全参考验证。

目前在线服务实现 NCBI Datasets、E-utilities 和 BLAST URLAPI。其他来源的基因组/注释可用于本地方式；本程序不会将一个 NCBI 组装号直接当作 Ensembl/UCSC 的版本号。

## 六、结果怎么看、常见问题怎么处理

第一次建议先打开这两个文件：

| 文件 | 用途 |
|---|---|
| `results/project1/report.html` | 双击用浏览器离线查看样本、候选引物和状态 |
| `results/project1/results.xlsx` | 用 Excel 查看、筛选和保存引物序列；引物按订购的 5′→3′ 方向输出 |

其他输出包括 CSV/TSV、完整 JSON、引物与扩增子 FASTA、有参考映射时的 BED，以及带校验和的清单。在线任务还会输出 `online_reference.json`，包含组装信息、获取的目标区域注释、请求来源和校验记录。确实在线获取了目标区域时，还会生成 `online_templates.fasta`（局部模板）和 `region_features.tsv`（区域注释表），方便复用。字段、坐标、热力学单位和高级参数见 [进阶说明](docs/advanced-usage.md)。

- `status=ok`：找到候选引物，不等于已经证明特异性或实验成功。
- `status=no_pair`：当前约束下没有候选；先查看 `message`，再检查模板长度、目标跨度、产物长度和 Tm 范围。
- `status=error`：这条记录有输入或运行问题，查看 `message`。
- `specificity_status=not_run`：没有执行特异性检查；场景一这是正常情况。
- `single_product_exact` / `single_product_blast`：在已记录的检索条件下找到单一预期产物。
- `single_product_unassigned_*`：只找到一个潜在产物，但原始 FASTA 没有预期基因组坐标，尚不能确认它就是你想要的位置。
- `multiple_products` / `multiple_products_remote`：存在多个潜在扩增产物，需要进一步筛选。
- `single_product_remote_limited`：在线返回的命中中只预测到单一预期产物，但全组装覆盖未验证。
- `single_product_unassigned_remote_limited`：单一在线潜在产物尚未确认就是你的目标位置。
- `incomplete_remote`：命中数达到上限、查询缺失或命中无法核对版本，检查不完整。
- `expected_not_found` 或特异性 `error`：检查参考版本、坐标或工具信息，不要直接把这对引物当成通过检查。

`screening_priority` 辅助区分本地单产物、在线有限证据、脱靶风险、未确认目标和检查不完整；它不是实验合格证。程序**保留所有候选及其检查结果，不会把有脱靶风险的候选悄悄删掉**；下单前应结合 JSON 中的检索范围、错配和产物信息人工复核，并做实验验证。BLAST 是启发式搜索，结果受参考完整性、错配模型和最大脱靶产物长度等限制。

常见终端问题：

| 提示/现象 | 怎么处理 |
|---|---|
| `primer-workflow: command not found` | 回到仓库目录激活 `.venv`，确认执行过安装命令 |
| 找不到 `data/project1/...` | 确认终端当前在仓库根目录，文件名和后缀与命令一致 |
| `output exists` | 换新输出目录，例如 `-o results/project1_run2`，保留旧结果 |
| `BLAST+ is required` | 按场景二安装并检查 `blastn -version`；不要以为已经完成检查 |
| `contig mismatch` / `coordinates outside genome` | 检查基因组和注释是否来自同一组装、染色体 ID 是否一致 |
| `gene/transcript not found` | 使用配套注释中的 ID；不要把其他数据库版本的 ID 直接混用 |
| 长基因没有引物 | 检查完整目标跨度；选择合理短区段，或有依据地调整产物长度 |

如果希望每条模板返回 3 对、产物长度为 150–400 bp，可以这样运行：

```bash
# 修改候选数和扩增长度；这里仍是场景一，不执行全基因组特异性检查
primer-workflow design data/project1/targets.fasta --num-return 3 --product-min-size 150 --product-max-size 400 -o results/project1_custom

# 读取已提供的参数 JSON；复制到项目后可以编辑，参数含义见文件旁的示例说明
primer-workflow design data/project1/targets.fasta --config data/examples/design-config.json -o results/project1_config
```

## 七、其他入口与验证

```bash
# 打开简易本机网页界面，适合粘贴序列或上传 FASTA/CSV/TSV
batch-primer-gui

# 运行自动测试；测试通过不代表你的真实样本已经实验验证
python -m unittest discover -s tests -v
```

**本地参考、在线功能和完整报告通过 `primer-workflow` 命令使用**，简易网页不提供本地参考或在线组装入口。需要将全基因组长期保存到本地时，也可以使用 NCBI 完整下载命令；手动参考包、完整下载、错配限制和其他进阶说明见 [进阶说明](docs/advanced-usage.md)。

支持不同物种的本地参考输入，不提供“保证无脱靶”的承诺。当前实现常规基因组 PCR；剪接 RT-qPCR、简并引物、探针和多重 PCR 尚未实现。
