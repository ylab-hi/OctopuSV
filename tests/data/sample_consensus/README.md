# OctopuSV 1.0 sample-consensus integration fixture

This fixture models the intended cohort workflow:

1. Each biological sample is first caller-merged (one SVCF file per sample).
2. Each per-sample SVCF contains multiple caller evidence blocks.
3. Those per-sample SVCFs are then merged with `octopusv merge --mode sample`.

The two inputs intentionally use different caller orders so sample-level genotype
synthesis can be tested for caller-order independence.

`baseline_old/sample_union_old.svcf` is a forensic snapshot of the historical
first-evidence sample-mode behavior. It is NOT the future OctopuSV 1.0 golden.
It exists only so intentional changes can be explained record by record.

## Events

| POS | Case | Future 1.0 consensus per sample |
| ---: | --- | --- |
| 1000 | all callers HET | GT=0/1, UC=3, UV=3 |
| 2000 | all callers HOM-alt | GT=1/1, UC=3, UV=3 |
| 3000 | HET/HOM disagreement; caller order reversed between samples | GT=1/., UC=3, UV=3 |
| 4000 | one caller missing | GT=0/1, UC=2, UV=2 |
| 5000 | same caller contributes two conflicting non-ref evidence blocks | GT=1/., UC=2, UV=2 |
| 6000 | haploid chrY ALT | GT=1, UC=2, UV=2 |
| 7000 | single-caller control | GT=0/1, UC=1, UV=1 |
| 8000 | explicit carrier/absent tie (adversarial) | GT=./., UC=1, UV=2 |

Notes:
- A caller missing the event is NO_VOTE, not ABSENT.
- POS 8000 intentionally includes an explicit 0/0 caller evidence block to test
  the rare presence-conflict path. It is adversarial, not the expected common case.
- sample-level AD will become unavailable in SVCF 1.1 synthesis; that change is
  not represented in the historical baseline.
