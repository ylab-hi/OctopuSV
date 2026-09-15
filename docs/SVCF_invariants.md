# SVCF source/evidence invariants for OctopuSV 0.5.0

This short file records implementation invariants used by the 0.5.0 code and
validator. It is not the SVCF 1.1 specification.

1. In caller-merge SVCF, the number of entries in `SOURCES` must equal the
   number of evidence columns on that record.
2. When `SOURCE_IDS` is present in caller-merge SVCF, `SOURCES[i]`,
   `SOURCE_IDS[i]`, and evidence column `i` refer to the same source record.
   A `.` entry in `SOURCE_IDS` is a positional missing-value placeholder; it
   occupies that position and must not be removed when parsing the list.
3. `SOURCES` may contain duplicate labels when one source contributes multiple
   evidence records to the same merged event.
4. Evidence count is not unique-source count. Caller/source-level voting and
   source-level statistics use unique sources; caller-mode evidence records are
   preserved individually.
5. Keys beginning with `_octopusv_` are runtime-only metadata and must never be
   serialized into SVCF output.
6. Sample mode is a per-sample summary layout. If one input sample contributes
   multiple evidence blocks to one merged event, OctopuSV keeps the first
   evidence column of that input file and reports the collapse. Therefore, if
   per-sample caller merges used different caller input orders, the retained
   caller may differ across samples.
