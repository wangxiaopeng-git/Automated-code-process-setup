# 输入数据目录

- `examples/` 是程序提供的教学示例，不是你需要自己生成的文件。先运行仓库首页的示例，再复制一份到自己的项目文件夹中修改。
- 建议新建 `data/project1/`，放自己的 `targets.fasta`、参考 `genome.fa` 和可选的 `genes.gtf` / `genes.gff3`。
- 在仓库根目录执行 `primer-workflow design data/project1/targets.fasta -o results/project1`。
- 有本地参考时追加 `--reference data/project1/genome.fa`，可再追加 `--gtf data/project1/genes.gtf` 或 `--gff data/project1/genes.gff3`。不需要手动建立参考包；默认自动用 BLAST 检查特异性。
- 不需要事先建立结果目录，程序会生成 `results/project1/`；如果它已存在，请换新名字，例如 `results/project1_run2`。
- 用户项目数据和隐藏的 `.primer-cache/` 不会默认提交到 Git。仅 `data/examples/` 与本说明随代码分发。

- 在线模式可以用 `--assembly GCF_000005845.2` 代替本地 `--reference`，按版本获取目标区域或在线筛查。输入 BED/数字 NCBI Gene ID 时只获取小范围序列和可用注释；缓存位于项目 `.ncbi-cache/`。在线覆盖范围有限，不能当作完整本地参考验证。
