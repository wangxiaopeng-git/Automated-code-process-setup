# Primer workflow implementation plan

Goal: Extend the uploaded Python/Primer3 prototype into a reproducible local-first, multi-species workflow.
Architecture: Keep the existing GUI and legacy CLI. Add an advanced CLI using shared parsers, reference bundles, annotation resolution, Primer3, explicit specificity checking and report exporters.
Tech stack: Python >=3.10, primer3-py 2.4.0, Biopython, openpyxl, requests; optional NCBI BLAST+.
Spec: docs/audit.md; user requests direct implementation and GitHub delivery.

## Constraints and decisions
Use the existing isolated checkout; do not create worktrees. Do not add large genomes or credentials to Git. Preserve provenance. Genomic intervals are internally 0-based half-open; exported primer coordinates explicitly identify 1-based conventions. No specificity run means not_run. BLAST-based checks enumerate inward-facing pairs and report their search limits; they are not experimental guarantees. First release supports genomic PCR; spliced RT-qPCR, multiplex, degenerate primers and VCF are documented future modes, not claimed implemented.

## Review focus
Malformed FASTA and annotations; negative-strand mappings; duplicate input identities; off-targets on either strand; untrusted remote ZIPs and partial downloads.

## Task 1: Import and harden existing prototype
- [x] Write regression tests for empty FASTA headers, zero records and nonfinite/negative chemistry parameters; run RED.
- [x] Fix core.py, cli.py and GUI diagnostic/count handling; run full unittest suite GREEN.

## Task 2: Inputs and reference bundles
- [x] Write tests for gzip FASTA, XLSX, GenBank, BED, GFF3/GTF, negative strand genes and validation failures; run RED.
- [x] Implement inputs.py and reference.py. Reference import validates identifiers and coordinates, copies source data into immutable bundles and records SHA256 manifest. Resolve targets into TemplateRecord plus coordinate context; run GREEN.

## Task 3: Specificity and report workflow
- [x] Write tests for paired hits, reverse strand amplicons, no check state, real Primer3 coordinates and output artifacts; run RED.
- [x] Implement specificity.py using BLAST+ short-query hits with exhaustive pairing within bounded configured product length; cap hits and mark incomplete results. Implement workflow.py, exporters.py and advanced CLI. Export CSV/TSV/XLSX/JSON/FASTA/HTML/BED and provenance; run GREEN.

## Task 4: NCBI and delivery
- [x] Write tests around HTTP boundary for NCBI search/download, unsafe archives and failure cleanup; run RED.
- [x] Implement ncbi.py with verified HTTPS, retry/backoff, cache checksum and assembly identity. Add CLI entry points, README, examples, CI and installation instructions.
- [x] Run full suite, package install/build, CLI functional checks and optional BLAST integration. Record live online limitations; obtain independent review, fix significant findings, push authorized feature branch and attempt PR.

## Execution record
39/39 tests passed with real BLAST+ and Primer3. Independent review findings fixed with failing-then-passing regressions. NCBI live access is blocked by the cloud proxy; transport/archive tests use controlled HTTP responses. Implementation performed inline per user request for direct generation. GitHub delivery uses feature/primer-workflow; remote push/PR outcome is reported in chat.
