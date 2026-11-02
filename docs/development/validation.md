# 本次交付验证

- 环境：Linux，Python 3.12；primer3-py 2.4.0；NCBI BLAST+ 2.16.0。
- unittest：39 项全部通过，无跳过（BLAST 已提供给测试进程）。
- 独立代码复核发现的末端错配漏检、短脱靶产物下界、空/部分校验清单误判，均补充回归测试并修正。
- 真实 Primer3 和 BLAST 集成：正负链参考坐标与引物序列回查、预期扩增子匹配、末端及两位错配位点恢复。
- CLI：10 kb 人工合成参考导入成功；BED+BLAST、基因 ID+精确检索均处理 2 个目标，生成 4 对候选，无检查错误；四对 BLAST 候选分别只有一个映射预期产物。
- 本机网页 HTTP 验证：页面读取成功，提交真实示例后返回 2 个成功样本、10 对候选及 CSV。
- 包：editable 安装、sdist 与 wheel 构建完成。
- NCBI：HTTP/归档导入、TaxID 不匹配、路径穿越、校验清单测试通过。云代理拒绝到 api.ncbi.nlm.nih.gov 的 CONNECT（403）；未进行实站下载验证。
- GitHub Actions 已配置 Python 3.10/3.12。这里仅验证本地 Python 3.12，不将配置 CI 等同于远程 CI 已通过。
- 跨系统安装、真实物种基因组的大规模性能、真实样本实验效果未验证。功能与检索限制见 README。

## 新手流程改进验证（2026-11-02）

- 自动测试增加到 50 项，Linux/Python 3.12 下全部通过，无跳过，真实 BLAST+ 2.16.0 可用。
- 新增覆盖：原始 FASTA + GFF/GTF 自动导入、负链、内容更新后的缓存失效、缓存复用、无参考默认跳过检查、有参考默认 BLAST、缺少 BLAST 的明确错误、显式跳过检查、注释参数冲突、项目默认输出。
- 逐条运行了首页/示例说明中的 11 个设计流程：单序列、项目复制、FASTA、CSV、TSV、Excel、GenBank、JSON 配置、基因 ID + exact、基因 ID + 默认 BLAST、BED。全部生成成功结果。
- 独立复核未发现重要问题；相对/绝对输入路径的项目命名边界已统一，并有回归测试。
- sdist/wheel 构建通过；检查源码包包含随附 Excel/GenBank/TSV/FASTA 和说明，不包含 `.primer-cache`、用户项目数据或结果。
- Windows PowerShell 和 macOS 命令提供了操作说明，但本次未在这两个系统上实测。NCBI 实站验证限制保持不变。
