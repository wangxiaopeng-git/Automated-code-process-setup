from __future__ import annotations

import argparse
import csv
import io
import json
import re
import tempfile
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from .core import DesignConfig, InputError, RESULT_FIELDS, design_batch, read_records


MAX_REQUEST_BYTES = 25 * 1024 * 1024


INDEX_HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>批量 PCR 引物设计</title>
  <style>
    :root { color-scheme: light; --ink:#172033; --muted:#657086; --line:#dce2eb;
      --blue:#1769e0; --blue2:#0d54bc; --bg:#f3f6fb; --card:#fff; --good:#087a50;
      --bad:#b13b32; }
    * { box-sizing:border-box; }
    body { margin:0; color:var(--ink); background:var(--bg); font:15px/1.45 -apple-system,
      BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif; }
    .top { background:linear-gradient(120deg,#102e63,#1769e0); color:#fff; padding:34px 24px 44px; }
    .wrap { max-width:1050px; margin:auto; }
    h1 { margin:0 0 5px; font-size:30px; letter-spacing:.2px; }
    .subtitle { opacity:.82; }
    main { margin-top:-22px; padding:0 24px 50px; }
    .card { background:var(--card); border:1px solid rgba(23,50,90,.09); border-radius:14px;
      box-shadow:0 8px 30px rgba(29,56,96,.09); padding:22px; margin-bottom:16px; }
    h2 { font-size:17px; margin:0 0 16px; }
    .drop { display:block; border:1.5px dashed #9eb3cf; border-radius:11px; padding:20px;
      background:#f8faff; text-align:center; cursor:pointer; transition:.15s; }
    .drop:hover { border-color:var(--blue); background:#f1f6ff; }
    input[type=file] { width:100%; max-width:580px; }
    .source-tabs { display:flex; gap:8px; margin-bottom:14px; }
    .source-tab { background:#edf2f9; color:#526079; font-weight:600; }
    .source-tab.active { background:var(--blue); color:#fff; }
    .source-panel { border:1px solid var(--line); border-radius:11px; padding:16px; background:#fafcff; }
    .paste-head { display:flex; justify-content:space-between; gap:12px; align-items:end; margin-bottom:9px; }
    .paste-head label { flex:1; max-width:300px; }
    input[type=text],textarea { width:100%; border:1px solid var(--line); border-radius:8px;
      padding:9px 10px; font:inherit; outline:none; background:#fff; }
    textarea { min-height:180px; resize:vertical; font:13px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
    input[type=text]:focus,textarea:focus { border-color:var(--blue); box-shadow:0 0 0 3px #1769e01a; }
    .char-count { color:var(--muted); font-size:12px; white-space:nowrap; padding-bottom:9px; }
    .hint { color:var(--muted); font-size:13px; margin:8px 0 0; }
    .grid { display:grid; grid-template-columns:repeat(4,minmax(130px,1fr)); gap:14px 16px; }
    label span { display:block; color:#46536a; font-size:13px; margin-bottom:5px; }
    input[type=number] { width:100%; border:1px solid var(--line); border-radius:8px;
      padding:9px 10px; font:inherit; outline:none; }
    input[type=number]:focus { border-color:var(--blue); box-shadow:0 0 0 3px #1769e01a; }
    .actions { display:flex; gap:10px; align-items:center; flex-wrap:wrap; }
    button,.download { border:0; border-radius:9px; padding:10px 18px; font:600 15px/1 inherit;
      cursor:pointer; text-decoration:none; display:inline-flex; align-items:center; min-height:40px; }
    button.primary { background:var(--blue); color:#fff; }
    button.primary:hover { background:var(--blue2); }
    button:disabled { opacity:.55; cursor:wait; }
    .download { color:var(--blue); background:#eaf2ff; }
    .hidden { display:none !important; }
    #status { color:var(--muted); }
    .spinner { width:17px;height:17px;border:2px solid #ffffff66;border-top-color:#fff;
      border-radius:50%;animation:spin .7s linear infinite;margin-right:8px; }
    @keyframes spin { to { transform:rotate(360deg); } }
    .summary { display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin:4px 0 18px; }
    .metric { padding:13px; background:#f6f8fc; border-radius:9px; }
    .metric b { display:block; font-size:23px; }
    .metric span { color:var(--muted); font-size:12px; }
    .table-wrap { overflow:auto; border:1px solid var(--line); border-radius:10px; }
    table { border-collapse:collapse; width:100%; min-width:920px; font-size:13px; }
    th,td { padding:9px 10px; border-bottom:1px solid #e8ecf2; text-align:left; white-space:nowrap; }
    th { position:sticky;top:0;background:#f4f7fb;color:#43506a;font-weight:600; }
    tr:last-child td { border-bottom:0; }
    .ok { color:var(--good);font-weight:600; }.error,.no_pair { color:var(--bad);font-weight:600; }
    .error-box { color:#8c2923;background:#fff0ee;border:1px solid #f1c5c1;border-radius:9px;
      padding:12px 14px;margin-top:14px;white-space:pre-wrap; }
    .privacy { color:var(--muted); font-size:12px; text-align:center; margin-top:20px; }
    @media(max-width:760px){.grid{grid-template-columns:repeat(2,1fr)}.summary{grid-template-columns:repeat(2,1fr)}}
  </style>
</head>
<body>
  <header class="top"><div class="wrap"><h1>批量 PCR 引物设计</h1>
    <div class="subtitle">Primer3 Core 2.6.1 · primer3-py 2.4.0 · 本地运行</div></div></header>
  <main class="wrap">
    <form id="designForm">
      <section class="card"><h2>1. 输入模板序列</h2>
        <div class="source-tabs">
          <button id="pasteTab" class="source-tab active" type="button">复制粘贴序列</button>
          <button id="fileTab" class="source-tab" type="button">上传文件</button>
        </div>
        <div id="pastePanel" class="source-panel">
          <div class="paste-head"><label><span>样本名称（粘贴单条裸序列时使用）</span>
            <input id="sampleId" type="text" value="pasted_sequence" maxlength="80"></label>
            <div id="charCount" class="char-count">0 个字符</div></div>
          <textarea id="sequenceText" spellcheck="false" placeholder=">sample_1&#10;ACGTTGCA...&#10;&#10;也可以直接粘贴一条不带 FASTA 标题的 DNA 序列"></textarea>
          <div class="hint">支持单条裸序列，或粘贴含多个样本的 FASTA 内容；空格和换行会自动处理。</div>
        </div>
        <div id="filePanel" class="source-panel hidden">
          <label class="drop"><input id="template" type="file" accept=".fasta,.fa,.fna,.csv,.tsv">
            <div class="hint">支持 FASTA、CSV 和 TSV；CSV/TSV 必须包含 id、sequence 列</div></label>
        </div>
      </section>
      <section class="card"><h2>2. 设置设计参数</h2><div class="grid">
        <label><span>候选引物对数</span><input name="num_return" type="number" min="1" value="5" required></label>
        <label><span>产物最短 (bp)</span><input name="product_min_size" type="number" min="1" value="100" required></label>
        <label><span>产物最长 (bp)</span><input name="product_max_size" type="number" min="1" value="300" required></label>
        <label><span>最长连续相同碱基</span><input name="max_poly_x" type="number" min="1" value="4" required></label>
        <label><span>引物最短 (nt)</span><input name="primer_min_size" type="number" min="1" value="18" required></label>
        <label><span>引物最优 (nt)</span><input name="primer_opt_size" type="number" min="1" value="20" required></label>
        <label><span>引物最长 (nt)</span><input name="primer_max_size" type="number" min="1" value="25" required></label>
        <label><span>引物允许 N 数</span><input name="max_ns" type="number" min="0" value="0" required></label>
        <label><span>最低 Tm (°C)</span><input name="primer_min_tm" type="number" step="0.1" value="57" required></label>
        <label><span>最优 Tm (°C)</span><input name="primer_opt_tm" type="number" step="0.1" value="60" required></label>
        <label><span>最高 Tm (°C)</span><input name="primer_max_tm" type="number" step="0.1" value="63" required></label>
        <label><span>单价盐 (mM)</span><input name="salt_monovalent_mm" type="number" step="0.1" value="50" required></label>
        <label><span>最低 GC (%)</span><input name="primer_min_gc" type="number" step="0.1" value="35" required></label>
        <label><span>最高 GC (%)</span><input name="primer_max_gc" type="number" step="0.1" value="65" required></label>
        <label><span>Mg²⁺ (mM)</span><input name="salt_divalent_mm" type="number" step="0.1" value="1.5" required></label>
        <label><span>dNTP (mM)</span><input name="dntp_mm" type="number" step="0.1" value="0.6" required></label>
      </div><p class="hint">CSV/TSV 可添加 target_start 和 target_length，指定扩增产物必须覆盖的 1-based 目标区。</p></section>
      <section class="card"><div class="actions">
        <button id="run" class="primary" type="submit"><span id="spin" class="spinner hidden"></span><span id="runText">开始设计</span></button>
        <a id="download" class="download hidden">下载完整 CSV</a><span id="status">等待开始</span>
      </div><div id="error" class="error-box hidden"></div></section>
    </form>
    <section id="results" class="card hidden"><h2>3. 设计结果</h2>
      <div id="summary" class="summary"></div><div class="table-wrap"><table><thead><tr>
        <th>样本</th><th>排名</th><th>正向引物 5'→3'</th><th>反向引物 5'→3'</th>
        <th>正向 Tm</th><th>反向 Tm</th><th>产物 bp</th><th>状态</th><th>说明</th>
      </tr></thead><tbody id="tbody"></tbody></table></div>
      <p class="hint">页面最多预览前 100 行；完整结果请下载 CSV。</p></section>
    <div class="privacy">所有序列仅在本机 127.0.0.1 处理，不会上传到互联网。</div>
  </main>
<script>
const form=document.getElementById('designForm'),run=document.getElementById('run'),spin=document.getElementById('spin');
const runText=document.getElementById('runText'),statusEl=document.getElementById('status'),err=document.getElementById('error');
const results=document.getElementById('results'),tbody=document.getElementById('tbody'),download=document.getElementById('download');
const pasteTab=document.getElementById('pasteTab'),fileTab=document.getElementById('fileTab');
const pastePanel=document.getElementById('pastePanel'),filePanel=document.getElementById('filePanel');
const sequenceText=document.getElementById('sequenceText'),template=document.getElementById('template');
let downloadUrl=null,inputMode='paste';
function textCell(row,value,cls=''){const td=document.createElement('td');td.textContent=value??'';if(cls)td.className=cls;row.appendChild(td)}
function setMode(mode){inputMode=mode;const paste=mode==='paste';pasteTab.classList.toggle('active',paste);fileTab.classList.toggle('active',!paste);pastePanel.classList.toggle('hidden',!paste);filePanel.classList.toggle('hidden',paste);statusEl.textContent=paste?'可直接粘贴 DNA 序列或 FASTA':'请选择 FASTA、CSV 或 TSV 文件'}
pasteTab.addEventListener('click',()=>setMode('paste'));fileTab.addEventListener('click',()=>setMode('file'));
template.addEventListener('change',()=>{if(template.files.length)setMode('file')});
sequenceText.addEventListener('input',()=>{document.getElementById('charCount').textContent=`${sequenceText.value.length.toLocaleString()} 个字符`});
form.addEventListener('submit',async(e)=>{e.preventDefault();
  run.disabled=true;spin.classList.remove('hidden');runText.textContent='正在设计…';statusEl.textContent='正在读取并计算';err.classList.add('hidden');results.classList.add('hidden');download.classList.add('hidden');
  try{const fd=new FormData(form),params={};for(const [k,v] of fd.entries())params[k]=v;
    let payload={input_mode:inputMode,params};
    if(inputMode==='paste'){const content=sequenceText.value.trim();if(!content)throw new Error('请粘贴模板序列');payload.content=content;payload.sample_id=document.getElementById('sampleId').value}
    else{const file=template.files[0];if(!file)throw new Error('请选择模板文件');payload.filename=file.name;payload.content=await file.text()}
    const response=await fetch('/api/design',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
    const data=await response.json();if(!response.ok||!data.ok)throw new Error(data.error||'设计失败');
    document.getElementById('summary').innerHTML=`<div class="metric"><b>${data.total}</b><span>总样本</span></div><div class="metric"><b>${data.success}</b><span>成功样本</span></div><div class="metric"><b>${data.failed}</b><span>未成功样本</span></div><div class="metric"><b>${data.row_count}</b><span>结果行数</span></div>`;
    tbody.replaceChildren();for(const item of data.rows){const tr=document.createElement('tr');textCell(tr,item.sample_id);textCell(tr,item.pair_rank);textCell(tr,item.left_sequence_5to3);textCell(tr,item.right_sequence_5to3);textCell(tr,item.left_tm_c);textCell(tr,item.right_tm_c);textCell(tr,item.product_size_bp);textCell(tr,item.status,item.status);textCell(tr,item.message);tbody.appendChild(tr)}
    if(downloadUrl)URL.revokeObjectURL(downloadUrl);downloadUrl=URL.createObjectURL(new Blob(['\ufeff',data.csv],{type:'text/csv;charset=utf-8'}));download.href=downloadUrl;download.download=data.output_name;download.classList.remove('hidden');results.classList.remove('hidden');statusEl.textContent=`完成：${data.success} 个成功，${data.failed} 个未成功`;results.scrollIntoView({behavior:'smooth'});
  }catch(ex){err.textContent=ex.message;err.classList.remove('hidden');statusEl.textContent='设计失败'}finally{run.disabled=false;spin.classList.add('hidden');runText.textContent='开始设计'}
});
</script></body></html>"""


def _config(params: Dict[str, Any]) -> DesignConfig:
    def integer(name: str, default: int) -> int:
        try:
            return int(params.get(name, default))
        except (TypeError, ValueError) as exc:
            raise InputError(f"参数 {name} 必须是整数") from exc

    def decimal(name: str, default: float) -> float:
        try:
            return float(params.get(name, default))
        except (TypeError, ValueError) as exc:
            raise InputError(f"参数 {name} 必须是数字") from exc

    config = DesignConfig(
        num_return=integer("num_return", 5),
        primer_min_size=integer("primer_min_size", 18),
        primer_opt_size=integer("primer_opt_size", 20),
        primer_max_size=integer("primer_max_size", 25),
        primer_min_tm=decimal("primer_min_tm", 57),
        primer_opt_tm=decimal("primer_opt_tm", 60),
        primer_max_tm=decimal("primer_max_tm", 63),
        primer_min_gc=decimal("primer_min_gc", 35),
        primer_max_gc=decimal("primer_max_gc", 65),
        product_min_size=integer("product_min_size", 100),
        product_max_size=integer("product_max_size", 300),
        max_poly_x=integer("max_poly_x", 4),
        max_ns=integer("max_ns", 0),
        salt_monovalent_mm=decimal("salt_monovalent_mm", 50),
        salt_divalent_mm=decimal("salt_divalent_mm", 1.5),
        dntp_mm=decimal("dntp_mm", 0.6),
    )
    config.validate()
    return config


def _input_from_payload(payload: Dict[str, Any]) -> tuple[str, str]:
    mode = str(payload.get("input_mode", "file"))
    content = payload.get("content")
    if not isinstance(content, str) or not content.strip():
        raise InputError("请输入模板序列或选择模板文件")
    content = content.strip()

    if mode == "paste":
        if content.startswith(">"):
            return "pasted_sequences.fasta", content + "\n"
        sample_id = str(payload.get("sample_id", "pasted_sequence")).strip()
        sample_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id).strip("_.-")
        sample_id = sample_id or "pasted_sequence"
        sequence = re.sub(r"\s+", "", content).upper()
        return f"{sample_id}.fasta", f">{sample_id}\n{sequence}\n"

    if mode != "file":
        raise InputError("未知的序列输入方式")
    filename = Path(str(payload.get("filename", ""))).name
    if not filename:
        raise InputError("请选择有效的模板文件")
    if Path(filename).suffix.lower() not in {".fasta", ".fa", ".fna", ".csv", ".tsv"}:
        raise InputError("文件类型必须是 FASTA、CSV 或 TSV")
    return filename, content


def _design(payload: Dict[str, Any]) -> Dict[str, Any]:
    filename, content = _input_from_payload(payload)
    if len(content.encode("utf-8")) > MAX_REQUEST_BYTES:
        raise InputError("模板内容不能超过 25 MB")

    config = _config(payload.get("params") or {})
    from primer3.bindings import design_primers

    with tempfile.TemporaryDirectory(prefix="primer-design-") as temporary:
        input_path = Path(temporary) / filename
        input_path.write_text(content, encoding="utf-8")
        records = list(read_records(input_path))
        if not records:
            raise InputError("输入文件中没有样本")
        rows = list(design_batch(records, config, design_primers))

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
    writer.writeheader()
    from .exporters import spreadsheet_safe

    writer.writerows({k: spreadsheet_safe(v) for k, v in row.items()} for row in rows)
    success_ids = {row["record_id"] for row in rows if row["status"] == "ok"}
    failed_ids = {row["record_id"] for row in rows if row["status"] != "ok"}
    return {
        "ok": True,
        "total": len(records),
        "success": len(success_ids),
        "failed": len(failed_ids),
        "row_count": len(rows),
        "rows": rows[:100],
        "csv": stream.getvalue(),
        "output_name": f"{Path(filename).stem}_primers.csv",
    }


class Handler(BaseHTTPRequestHandler):
    def _json(self, status: int, payload: Dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        if self.path != "/":
            self.send_error(404)
            return
        body = INDEX_HTML.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self) -> None:
        if self.path != "/api/design":
            self.send_error(404)
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > MAX_REQUEST_BYTES + 1_000_000:
                raise InputError("请求为空或文件过大")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            self._json(200, _design(payload))
        except (InputError, OSError, UnicodeError, json.JSONDecodeError) as exc:
            self._json(400, {"ok": False, "error": str(exc)})
        except Exception as exc:
            self._json(500, {"ok": False, "error": f"Primer3 运行失败：{exc}"})

    def log_message(self, format: str, *args: object) -> None:
        return


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Start the local primer design web GUI"
    )
    parser.add_argument(
        "--port", type=int, default=0, help="local port; 0 chooses a free port"
    )
    parser.add_argument(
        "--no-browser", action="store_true", help="do not open a browser"
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    url = f"http://127.0.0.1:{server.server_port}/"
    print(f"批量引物设计界面：{url}", flush=True)
    print("使用期间请保留此窗口；按 Control+C 可停止。", flush=True)
    if not args.no_browser:
        threading.Timer(0.3, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
