# Online workflow verification record

Version: 0.3.0. User instructions are README.md and docs/advanced-usage.md; this is a maintainer verification record.

- Linux / Python 3.12 / primer3-py 2.4.0 / local BLAST+ 2.16.0: 77 tests passed, no skips.
- Controlled HTTP integration covers assembly/sequence metadata, pagination completeness, exact sequence versions, bounded GenBank regions and regional features, numeric Gene IDs, positive/negative strand mapping, deduplicated query batches, remote XML and L/R products, cache reuse/corruption, rate limits and resumed RID polling deadlines.
- Independent review found merged nt Hit exclusion, incorrect persisted poll dispatch time, rate-limit sleeps beyond deadlines, and automatic GET retries bypassing throttling. Each received a failing-then-passing regression and a full green suite.
- Online metadata inspection now returns actual contig accessions, aliases and lengths, not only their hash. Region outputs include reusable template FASTA and feature TSV; no full-genome or whole-assembly coverage claim is made.
- Existing installed template CLI and local-reference CLI smoke checks both processed two records and four candidate pairs successfully.
- Editable installation 0.3.0, sdist and wheel built successfully. Distribution inspection confirms online modules and examples are present; runtime caches, private inputs/results and maintainer records are excluded.
- Official NCBI Datasets OpenAPI from ncbi/datasets/master was inspected for endpoint/field schemas. api.ncbi.nlm.nih.gov, eutils.ncbi.nlm.nih.gov and blast.ncbi.nlm.nih.gov were rejected by the environment proxy. Additive domain requirements were saved to the environment draft but remain unapplied (runtime spec revision 1). No claim of live NCBI/BLAST end-to-end compatibility is made; online mode remains explicitly experimental.
- Windows/macOS and remote GitHub Actions were not independently validated in this execution.
