## Source/evidence mapping fix in 0.5.0

OctopuSV 0.5.0 fixes cases in earlier releases where affected merged SVCF
records could contain incorrect `SOURCES` / `SOURCE_IDS` ordering, or could
lose additional evidence from the same source in caller-mode output.

Existing merged SVCF files are not rewritten automatically. To check an
existing file, run:

```bash
octopusv validate-svcf -i file.svcf
```

If validation reports `E_SRC_004` or `E_FMT_002` on an older merged SVCF,
re-run the original merge with OctopuSV 0.5.0. Do not manually repair the
positional source/evidence fields.
