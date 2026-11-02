# 安装版本和命令路径排查

## 判断问题属于哪一层

如果提示 `unrecognized arguments: --gff`、`unrecognized arguments: --assembly`，或 `reference inspect` 不在命令列表中，程序还没有进入序列读取或 NCBI 请求阶段。首先核对实际加载的 CLI，不能靠重新下载基因组解决。

- 0.2.0 没有这些参数，不能执行相应命令。
- 0.3.0 支持 `--gff`、`--assembly`、`reference inspect`。
- 0.3.1 起提供 `primer-workflow --version`。较早版本不识别 `--version` 时，用下面的 Python 命令核对。

终端显示 `(.venv)` 只说明激活过一个同名环境，不证明它属于当前项目。`(.venv) (base)` 同时显示本身不是错误。多个项目副本（例如一个名为 `Batch-primer-design-workflow`，另一个名为 `Batch-primer-design-workflow-main`）各自的环境和可编辑安装可能指向不同目录。

## 1. 显示实际使用的程序和源码

在报错的同一终端运行：

```bash
# 显示当前所在目录，确认这里有 pyproject.toml 和 src 文件夹
pwd

# 显示 shell 当前会调用哪个命令，是否存在多个同名命令
type -a primer-workflow

# 显示实际解释器、包版本和源码路径；不会访问网络
python -c "import sys, batch_primer_design, batch_primer_design.advanced_cli as cli; print('python:', sys.executable); print('version:', batch_primer_design.__version__); print('cli:', cli.__file__)"
```

源码应来自你准备运行的项目的 `src/batch_primer_design/advanced_cli.py`；如果解释器或源码指向另一份项目，即使当前目录里的示例和说明是新的，命令仍会按照另一份源码解析参数。

## 2. 让当前项目与它自己的环境绑定（macOS/Linux）

以下操作保留 `data`、`results`，不删除环境，不重置 Git。请先进入你实际准备使用、含 `pyproject.toml` 的那个项目根目录。

```bash
# 确认本地分支和文件改动；如果使用源码 ZIP 而非 Git 下载，跳过 Git 命令
git status --short --branch

# Git 用户：确认当前是 main 后再执行；需要切换时先保存改动，再用 git switch main
# 只允许快进，不覆盖本地提交；若提示有冲突/不能快进，先保留现场，不要执行 reset --hard
git pull --ff-only origin main

# 如果本目录还没有 .venv，才运行这一行；需要 Python 3.10+，建议 3.12
# python3 -m venv .venv

# 明确使用本目录的解释器重新安装，不依赖 PATH 指向的另一个 python/pip
./.venv/bin/python -m pip install --upgrade --no-cache-dir -e .

# 用同一个解释器确认版本和源码来自本项目
./.venv/bin/python -c "import sys, batch_primer_design, batch_primer_design.advanced_cli as cli; print(sys.executable); print(batch_primer_design.__version__); print(cli.__file__)"

# 验证本地参考参数存在：应看到 --gff 和 --assembly
./.venv/bin/python -m batch_primer_design.advanced_cli design --help

# 验证在线元数据命令存在：应列出 inspect
./.venv/bin/python -m batch_primer_design.advanced_cli reference --help
```

如果使用 ZIP，应从含所需参数的代码版本取得源码，并在该目录安装；`pip install -e .` 只会安装当前文件夹，不会自动替你从 GitHub 获取新源码。不要直接复制另一份项目的 `.venv`。

## 3. 用明确解释器先运行，再恢复简短命令

```bash
# 本地示例：绕过可能仍指向旧目录的 primer-workflow 命令
./.venv/bin/python -m batch_primer_design.advanced_cli design data/project2/genes.txt --reference data/project2/genome.fa --gff data/project2/genes.gff3 --specificity exact -o results/project2
```

如果结果目录已存在，换成新目录如 `results/project2_check`，不要为了重跑而删除旧结果。

确认成功后，若当前有已激活的旧环境，先执行 `deactivate`，再在当前项目目录执行：

```bash
# 激活当前项目的环境
source .venv/bin/activate

# macOS 默认 zsh 刷新命令路径；Bash 使用 hash -r
rehash

# 再核对命令路径，应该指向当前项目的 .venv
type -a primer-workflow
```

Windows PowerShell 可以用 `.\.venv\Scripts\python.exe -m batch_primer_design.advanced_cli ...` 做同样的显式解释器检查。

## 4. 参数恢复后仍有错误时

`input file does not exist`、网络错误、`no_pair` 等属于后续阶段，需要根据新错误处理。不要把参数识别恢复等同于 NCBI 实站或实验验证成功。

复制命令时避免文件名末尾出现不换行空格等隐藏字符。例如 `genes.txt` 和 `genes.txt ` 可能是两个不同文件。可以在 Python 中显示文件名的精确形式：

```bash
# 输出项目目录内文件名的 repr，便于识别末尾隐藏空格
./.venv/bin/python -c "from pathlib import Path; print([repr(p.name) for p in Path('data/project3').iterdir()])"
```
