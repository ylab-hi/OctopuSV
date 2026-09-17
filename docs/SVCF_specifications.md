# SVCF 1.1 Specification: A VCF-Based Intermediate Format for Structural Variant Processing and Integration

**Specification version:** 1.1
**Status:** Locked stable contract for OctopuSV 1.0
**Reference implementation:** OctopuSV 1.0
**Recommended extension:** `.svcf`


## 1. Introduction

SVCF is a VCF-based intermediate format for structural variant normalization, integration, source tracking, sample synthesis, and downstream analysis. It is defined by this specification and implemented by OctopuSV.

SVCF keeps the line-oriented text structure and core columns of the [Variant Call Format (VCF) Version 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf), but it adds conventions that standard VCF does not provide for representing the relationships among:

- source caller records (observations);
- merged structural-variant events;
- biological sample-level calls.

SVCF is therefore an OctopuSV intermediate representation rather than a replacement for the official VCF specification. **SVCF 1.1 is not guaranteed to be a conforming VCF file**, especially in caller mode: caller-mode records may contain a record-specific number of evidence blocks after `FORMAT` even though the `#CHROM` header contains one trailing label. General VCF software may therefore read only parts of an SVCF file or may reject it. Software that requires conventional VCF sample columns should consume the VCF produced by `octopusv svcf2vcf`, not SVCF directly.

The recommended file extension is:

```text
.svcf
```

Before using SVCF with software that expects conventional VCF sample columns, convert it with:

```bash
octopusv svcf2vcf -i input.svcf -o output.vcf
```

This document defines **SVCF 1.1**, the versioned SVCF contract implemented by OctopuSV 1.0.

The words **must**, **should**, and **may** are used in their ordinary specification sense:

- **must** means that the rule is required;
- **should** means that the rule is recommended unless there is a clear reason not to follow it;
- **may** means that the item is optional.

### 1.1 Terminology and semantic layers

SVCF 1.1 distinguishes three concepts that must not be conflated:

- **Observation / evidence record:** one structural-variant record reported by one source caller.
- **Merged event:** the OctopuSV event obtained after record matching/grouping. A merged event may contain multiple observations, including multiple observations from the same caller.
- **Sample call:** the biological-sample-level state synthesized for one merged event in multi mode.

Caller mode stores observations associated with a merged event. Multi mode stores synthesized sample calls. A representative event record is not a substitute for the complete caller evidence, and a caller evidence block is not a substitute for a synthesized sample call.

A **source** is the explicit caller/input identity bound to an evidence block. A **caller vote** is a caller-level state used during sample synthesis; multiple evidence records from the same source still contribute at most one caller vote.

---

## 2. File identity and versioning

### 2.1 Versioned SVCF 1.1

A versioned SVCF 1.1 file must declare both its SVCF version and its data model:

```text
##SVCFVersion=1.1
##OctopuSV_mode=caller
```

or:

```text
##SVCFVersion=1.1
##OctopuSV_mode=multi
```

The first meta-information line should remain:

```text
##fileformat=VCFv4.2
```

OctopuSV normally also writes:

```text
##source=OctopuSV
##fileDate=...
##OctopuSV_WARNING=This is SVCF format. Use 'octopusv svcf2vcf' to change back to standard VCF format before bcftools/vcftools
```

Once an SVCF version is declared, readers must not infer a different data model from the number of columns, the presence of `SOURCES`, filenames, ID prefixes, or other heuristics.

For SVCF 1.1:

- `##OctopuSV_mode=caller` requires the caller FORMAT defined in Section 6.1;
- `##OctopuSV_mode=multi` requires the synthesized sample FORMAT defined in Section 6.2;
- the declared version and mode must agree with every record in the file;
- conflicting version or mode declarations are invalid;
- a file declaring an unsupported future version (for example `1.2`) must not be silently interpreted as SVCF 1.1.

### 2.2 Legacy SVCF

Files without `##SVCFVersion` are **legacy SVCF**. This includes the unversioned SVCF layout used by earlier OctopuSV releases and described in the original OctopuSV publication. Legacy SVCF predates the versioned SVCF 1.1 contract; this specification does not retroactively assign those files a formal `SVCFVersion=1.0` identity. OctopuSV may continue to read legacy files through compatibility paths, but an unversioned file is not guaranteed to satisfy the SVCF 1.1 contract.

In particular, historical multi-sample output from `octopusv correct` may contain:

```text
##OctopuSV_mode=multi
```

while still using caller-style evidence blocks rather than the synthesized SVCF 1.1 sample schema. Such files remain legacy/unversioned and must not be interpreted as versioned SVCF 1.1 multi files.

### 2.3 Compatibility and version-bump policy

SVCF version numbers describe the file contract, not the OctopuSV software version. Within SVCF 1.1:

- the two mode-specific FORMAT schemas and their field order are stable;
- the meanings of the core fields defined by this specification are stable;
- additional optional meta-information or INFO annotations may be added only when they do not change the interpretation of existing SVCF 1.1 records;
- a change that alters a required FORMAT field, field order, positional binding rule, or the meaning of an existing core field requires a new SVCF specification version.

Readers may ignore unknown optional meta-information or INFO annotations when doing so is safe, but they must not ignore an unsupported `SVCFVersion` or a mode/FORMAT conflict.

**SVCF 1.1 is a locked contract.** If OctopuSV implementation behavior disagrees with a normative SVCF 1.1 rule, the implementation must be investigated first; the SVCF 1.1 specification must not be changed merely to match implementation drift. An intentional change to a required schema rule, positional binding rule, coordinate rule, missing-value meaning, or other core semantic defined here requires a new SVCF specification version.

---

## 3. SVCF 1.1 data models

SVCF 1.1 has two versioned data models.

| Mode | Marker | Data model | Typical producer |
|---|---|---|---|
| Caller | `##OctopuSV_mode=caller` | Evidence-preserving caller observations associated with one merged event | single-sample `octopusv correct`; `octopusv merge --mode caller` |
| Multi | `##OctopuSV_mode=multi` | Fixed biological-sample matrix of synthesized sample-level calls | `octopusv merge --mode sample` |

The mode describes the **data model**, not merely the command that produced the file.

### 3.1 Caller mode

Caller mode is the evidence-preserving layer.

The `#CHROM` header must contain exactly one trailing column after `FORMAT`. That header label is not an enumeration of all caller evidence blocks. Individual caller-mode records may contain one or more evidence blocks, depending on how many source records support that merged event.

Two common caller-mode forms are valid:

1. **single-evidence caller SVCF**, normally produced by single-sample `octopusv correct`;
2. **caller-merge SVCF**, where one event may retain evidence from multiple source records/callers.

A versioned caller file must use the caller FORMAT defined in Section 6.1 on every record.

A caller SVCF should not mix direct single-evidence records that omit `SOURCES` with caller-merge records that use `SOURCES` within the same file. OctopuSV validates these layouts separately.

### 3.2 Multi mode

Multi mode is the synthesized sample layer.

The names after `FORMAT` in the `#CHROM` header are biological sample names and define one fixed sample-column order for the entire file. Sample names should be unique within the header so that each column has an unambiguous biological-sample identity.

Every record must contain exactly one sample block for each declared sample column and must use the multi FORMAT defined in Section 6.2.

A multi-mode sample block is not a raw caller record. It is a sample-level summary derived after caller evidence has already been grouped into the merged event.

`##OctopuSV_mode=multi` in **versioned SVCF 1.1** therefore means the synthesized output of `octopusv merge --mode sample`. Historical unversioned files carrying the same marker remain legacy and are not covered by this rule.

---

## 4. File structure

An SVCF file contains:

1. meta-information lines beginning with `##`;
2. one header line beginning with `#CHROM`;
3. tab-delimited data lines.

The first nine columns are fixed and must appear in this order:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
```

At least one column must follow `FORMAT`.

A versioned writer must emit VCF meta-information definitions for SVCF-specific INFO/FORMAT fields that it writes, with Number/Type compatible with the semantics in this specification. Unknown optional meta-information lines are permitted and should be preserved by structure-preserving tools.

### 4.1 Caller-mode header shape

SVCF 1.1 caller mode requires exactly one trailing header label:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
```

The name may be an original input sample name rather than the literal word `SAMPLE`.

Caller-mode data records may nevertheless contain more than one evidence block. This variable-width record layout is intentional and is not a conventional VCF sample matrix.

### 4.2 Multi-mode header shape

SVCF 1.1 multi mode requires at least one biological sample column:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample_1  sample_2  ...
```

Every data record must contain exactly that number of sample blocks.

---

## 5. Fixed fields and INFO fields

The first nine columns describe the normalized or merged event. Evidence/sample blocks retain source- or sample-specific information.

| Field | SVCF rule |
|---|---|
| `CHROM` | Contig containing the event start or first breakpoint. |
| `POS` | Positive, 1-based coordinate of the event start or first breakpoint. |
| `ID` | Identifier of the representative SVCF event record. |
| `REF` | Reference allele of the representative record. |
| `ALT` | Alternate allele of the representative record. `BND` requires valid VCF breakend notation. `TRA` may use valid VCF breakend notation when orientation is known, or symbolic `<TRA>` when `CHR2` and numeric `END` define the remote breakpoint and orientation is unknown. |
| `QUAL` | Representative quality value, or `.` when unavailable. |
| `FILTER` | Filter status of the representative event. |
| `INFO` | Event-level annotations and source/sample relationship fields. |
| `FORMAT` | The exact mode-specific SVCF schema. |

Every SVCF record must contain these INFO keys, even when the value is `.`:

```text
SVTYPE
END
SVLEN
CHR2
SUPPORT
SVMETHOD
RTID
AF
STRAND
RNAMES
```

Current OctopuSV writers also use:

```text
SOURCES
SOURCE_IDS
```

### 5.1 INFO field meanings

| Field | Meaning |
|---|---|
| `SVTYPE` | Structural variant type: `DEL`, `DUP`, `INV`, `INS`, `TRA`, or `BND`. |
| `CHR2` | Contig containing the second coordinate or mate breakpoint. |
| `END` | Second coordinate used by the SVCF record. Its meaning depends on `SVTYPE`. |
| `SVLEN` | Positive event length when applicable; `TRA`/`BND` use `.`. |
| `SUPPORT` | Read-support value associated with the representative event. It is not the number of callers or samples. |
| `SVMETHOD` | Method that produced the current SVCF event, normally `OctopuSV`. |
| `RTID` | Related/reciprocal record ID when available. |
| `AF` | Allele frequency when available. |
| `STRAND` | Event strand/orientation value when available. |
| `RNAMES` | Supporting read names when available. |
| `SOURCES` | Ordered source labels associated with a merged record. |
| `SOURCE_IDS` | Ordered original record IDs used for source-level traceability. |

`SUPPORT` may be `.` or a non-negative integer.

---

## 6. Mode-specific FORMAT schemas

`CO` must remain the final FORMAT field in both SVCF 1.1 schemas. This permits robust parsing of source IDs and ALT strings that themselves contain colons.

### 6.1 Caller FORMAT

SVCF 1.1 caller mode uses exactly:

```text
GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO
```

The order must not change.

| Field | Meaning in caller mode |
|---|---|
| `GT` | Genotype from the source evidence record. |
| `AD` | Source reference/alternate allele depths when available. |
| `LN` | Absolute source-event length when available. |
| `ST` | Source strand/orientation value. |
| `QV` | Source quality value. |
| `TY` | Source structural-variant type. |
| `ID` | Original source record ID. |
| `SC` | Source caller/method label. |
| `REF` | Original source REF. |
| `ALT` | Original source ALT, including original BND syntax when applicable. |
| `CO` | Source coordinates. |

### 6.2 Multi FORMAT

SVCF 1.1 multi mode uses exactly:

```text
GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO
```

The order must not change.

The final five fields remain:

```text
ID:SC:REF:ALT:CO
```

so the shared SVCF block parser can continue to parse colon-containing IDs, BND ALT values, and symbolic ALT strings safely.

| Field | Meaning in multi mode |
|---|---|
| `GT` | Synthesized sample-level genotype. |
| `AD` | `.,.` for synthesized sample calls; caller allele depths are not composable across callers. |
| `UC` | Number of unique callers supporting carrier presence. |
| `UV` | Number of unique callers contributing a valid presence vote. |
| `LN` | Length associated with the sample/event representation. |
| `ST` | Representative strand/orientation value. |
| `QV` | Representative quality value. |
| `TY` | Structural-variant type. |
| `ID` | Sample input-event ID retained by OctopuSV. |
| `SC` | `OctopuSV` for synthesized sample calls. |
| `REF` | Representative REF. |
| `ALT` | Representative ALT. |
| `CO` | Representative/source coordinate string when available. |

---

## 7. Source/evidence binding

### 7.1 Caller-merge positional contract

In caller-merge records, source identity is explicit and positional.

If a record contains:

```text
SOURCES=sourceA,sourceA,sourceB
SOURCE_IDS=idA1,idA2,idB1
```

then the record means:

```text
SOURCES[1]    <-> SOURCE_IDS[1]    <-> evidence block 1
SOURCES[2]    <-> SOURCE_IDS[2]    <-> evidence block 2
SOURCES[3]    <-> SOURCE_IDS[3]    <-> evidence block 3
```

Duplicate `SOURCES` values are legal. One caller may contribute multiple evidence records to one merged event.

The number of evidence blocks is therefore **not** the number of unique callers.

When `SOURCES` is present on a caller-merge record:

- its item count must equal the number of evidence blocks;
- duplicate source labels must be preserved;
- source order is evidence-block order.

OctopuSV normally emits caller evidence in deterministic input-source order. Consumers must nevertheless rely on the positional binding itself rather than assigning biological meaning to the order.

When `SOURCE_IDS` is present:

- its item count must equal the number of `SOURCES` items;
- `.` is a positional missing-ID placeholder and must not be dropped;
- each `SOURCE_IDS[i]` must match the `ID` stored in evidence block `i`.

Software must not reconstruct source identity from filenames, ID prefixes, sample ordering, or other heuristics when explicit source/evidence binding is available.

`SOURCES` and `SOURCE_IDS` are delimiter-based positional lists in SVCF 1.1. The format does not define an escaping mechanism for list/INFO delimiters. Therefore:

- a `SOURCES` item **must** be non-empty and **must not** be `.`;
- `.` is permitted in `SOURCE_IDS` only as the positional missing-ID placeholder;
- non-missing `SOURCES` and `SOURCE_IDS` items **must not** contain `,`, `;`, `=`, or whitespace.

These characters are reserved because comma separates positional list items, semicolon separates INFO fields, equals separates INFO keys from values, and whitespace is not a portable VCF token character. A conforming writer must reject such values rather than silently sanitize, truncate, split, or escape them. Readers interpret delimiters structurally; because SVCF 1.1 defines no escaping layer, an intended identifier that already contains one of these reserved delimiters cannot be losslessly reconstructed after serialization. A future encoding that permits such values would require an explicit specification change.

### 7.2 Single-evidence caller records

A direct single-sample `octopusv correct` output contains one caller evidence block and may omit `SOURCES`/`SOURCE_IDS`. There is no cross-evidence source ambiguity in such a record.

### 7.3 Multi-mode `SOURCES`

In SVCF 1.1 multi mode, the `#CHROM` header defines the fixed biological-sample column order.

`SOURCES` is event-level metadata listing samples/input files that contributed evidence to the merged event. It is not a replacement for the fixed sample-column order and its item count does not define the number of sample columns.

When `SOURCE_IDS` is present in multi mode, its positions correspond to the positions in `SOURCES`. Multi mode is a synthesis layer: these IDs identify the retained per-sample input-event representation used by OctopuSV and are **not an exhaustive lossless list of all caller evidence IDs** that contributed upstream. Applications that require complete caller evidence must retain/use the caller-mode SVCF.

---

## 8. Sample synthesis semantics

SVCF 1.1 multi mode is a defined synthesis layer rather than a raw evidence layer.

### 8.1 One vote per unique caller

Multiple evidence records from the same caller are first reduced to one caller state. They do not count as multiple independent caller votes.

Caller-level states include the concepts of:

```text
NO_VOTE
ABSENT
CARRIER_HET
CARRIER_HOM
CARRIER_UNKNOWN
```

Examples:

- all usable evidence from one caller is `0/0` -> `ABSENT`;
- all is `0/1` -> `CARRIER_HET`;
- all is `1/1` -> `CARRIER_HOM`;
- `0/1` and `1/1` from the same caller -> carrier present but zygosity unresolved;
- contradictory absence and carrier evidence from the same caller -> that caller contributes no valid presence vote.

Haploid `0` and `1` are supported as absence and carrier states respectively.

The current consensus model is biallelic. Genotypes using other ALT allele indices are not promoted into a biallelic consensus claim.

### 8.2 Presence consensus

Let:

```text
C = number of valid carrier caller states
A = number of valid absent caller states
```

Then:

```text
C > A  -> carrier
A > C  -> absent
C == A and C+A > 0 -> ./.
no valid votes       -> ./.
```

`UC` records the number of unique carrier callers (`C`).

`UV` records the total number of unique callers with a valid presence vote (`C + A`).

### 8.3 Zygosity consensus

When carrier presence is established, zygosity is resolved separately.

If carrier callers agree on heterozygosity, OctopuSV writes:

```text
0/1
```

If they agree on homozygous ALT, OctopuSV writes:

```text
1/1
```

If carrier presence is established but HET/HOM evidence conflicts, OctopuSV writes:

```text
1/.
```

This means that at least one ALT allele is supported while the second allele is unresolved. OctopuSV does not replace this with `0/1` or `./.` because either would assert more or less than the evidence supports.

If all supporting carrier states are haploid, the synthesized genotype may be:

```text
1
```

### 8.4 Sample-level allele depth

Synthesized multi-mode calls use:

```text
AD=.,.
```

Caller allele depths are not summed, averaged, or otherwise combined because different callers may use overlapping reads and different support definitions.

---

## 9. Unobserved sample placeholders

SVCF 1.1 multi mode uses fixed-width sample columns. When a sample has no input event contributing to a merged event, OctopuSV writes an internal layout placeholder with the characteristic values:

```text
GT=0/0
AD=.,.
UC=0
UV=0
ID=.
SC=.
```

The remaining unavailable fields are written as `.`.

This placeholder is **not** an evidence-backed homozygous-reference call. It represents an unobserved event in that sample within the fixed SVCF matrix.

This is distinct from an evidence-backed absence call such as:

```text
GT=0/0
UC=0
UV=2
SC=OctopuSV
```

and from an unresolved call with evidence, such as:

```text
GT=./.
UC=0
UV=0
SC=OctopuSV
```

Downstream conversion policy is described in Section 12.

### 9.1 Missing and partial genotype semantics

SVCF 1.1 uses VCF-style missing-value notation deliberately:

| Value | Meaning in SVCF 1.1 |
|---|---|
| `.` | Scalar value unavailable / not defined. |
| `.,.` | Two-component allele depth unavailable; it does **not** mean zero depth. |
| `./.` | Diploid genotype unresolved / unavailable. |
| `1/.` | At least one ALT allele is supported; the second allele is unresolved. |
| `0/0` with `UV>0` | Evidence-backed absence call in synthesized multi mode. |
| internal `0/0` with `UC=0`, `UV=0`, `ID=.`, `SC=.` | Fixed-width unobserved-event placeholder; not an evidence-backed genotype. |

Software must not silently convert these missing/partial states into numeric zero or a more specific genotype unless an explicit conversion policy says to do so.

---

## 10. Coordinate field (`CO`) and colon-containing values

`CO` has the form:

```text
startChrom_startPos-endChrom_endPos
```

Examples:

```text
1_10889-1_10936
1_3845267-hs37d5_32469995
```

Contig names may contain underscores or hyphens. Parsers must not split coordinate strings using the first underscore or first hyphen without validating the resulting coordinates.

SVCF 1.1 does **not** permit `:` inside contig names used by record-level `CHROM` or `INFO/CHR2`. The fixed evidence block uses `:` as its field delimiter and `CO` embeds those contig names without an escaping layer, so a colon-bearing contig cannot currently be serialized and parsed losslessly. A conforming writer must reject such records rather than emit an ambiguous evidence block. Full support for colon-bearing contig names requires a future explicit encoding/specification change.

Source record IDs may contain colons. BND ALT strings and symbolic ALT values may also contain colons, for example:

```text
MantaDEL:469174:0:1:0:0:0
N]chr2:12345]
<INS:ME:ALU>
```

Software must therefore not assume that every colon inside an evidence/sample block is a FORMAT separator. OctopuSV keeps `ID:SC:REF:ALT:CO` as the final five FORMAT fields and uses a shared structure-aware parser.

---

## 11. Structural variant representation

Legal SVCF structural-variant types are:

```text
DEL
DUP
INV
INS
TRA
BND
```

`END` is the second coordinate used by the SVCF record. Its meaning depends on `SVTYPE`.

For `BND`, the mate coordinate is encoded in breakend `ALT` and must agree with `CHR2` and `END`. For `TRA`, the remote breakpoint may be encoded either in breakend `ALT` or, for symbolic `<TRA>`, by `CHR2` and `END`.

### 11.1 DEL, DUP, and INV

For `DEL`, `DUP`, and `INV`:

- `END` must be a numeric value;
- `END` must be greater than or equal to `POS`;
- `CHR2` is normally the same as `CHROM`;
- `SVLEN` is a positive event length when available.

Current OctopuSV output uses the absolute length rather than a negative deletion length.

### 11.2 INS

For `INS`, the keys `END`, `SVLEN`, and `CHR2` must be present.

When the insertion length is known, current OctopuSV SVCF output normally represents its internal span as:

```text
END = POS + SVLEN
```

For example:

```text
POS=10889
END=10936
SVLEN=47
```

When the value is unavailable, `END` or `SVLEN` may be `.`. The current validator requires the keys but does not require a numeric relationship for `INS`.

The internal SVCF endpoint is used by OctopuSV processing and merging. During `svcf2vcf` conversion, an insertion is written with the conventional VCF endpoint:

```text
END = POS
```

### 11.3 TRA and BND

`TRA` and `BND` both represent events with two breakpoints, but SVCF 1.1 does not require the same ALT representation for both types.

#### 11.3.1 BND

`BND` requires valid VCF breakend notation in `ALT`. The accepted forms are:

```text
t[chr:pos[
t]chr:pos]
[chr:pos[t
]chr:pos]t
```

where `t` is sequence placed before or after the breakend expression.

For `BND`:

- `CHR2` must contain the mate contig;
- `END` must contain the numeric mate position;
- the mate contig and position encoded in `ALT` must agree with `CHR2` and `END`;
- `SVLEN` must be `.`.

A retained `BND` record is valid SVCF and is not considered a conversion failure.

#### 11.3.2 TRA

`TRA` requires two known breakpoint coordinates but does not require known breakend orientation. Orientation is additional evidence rather than a prerequisite for representing the event.

A `TRA` record may therefore use either of two representations.

**Orientation known: breakend ALT**

```text
CHROM=1
POS=3845267
ALT=C[hs37d5:32469995[
CHR2=hs37d5
END=32469995
SVTYPE=TRA
SVLEN=.
```

When `TRA` uses breakend notation:

- `CHR2` must contain the mate contig;
- `END` must contain the numeric mate position;
- the mate contig and position encoded in `ALT` must agree with `CHR2` and `END`;
- `SVLEN` must be `.`.

**Orientation unknown: symbolic `<TRA>` ALT**

```text
CHROM=1
POS=3845267
ALT=<TRA>
CHR2=hs37d5
END=32469995
SVTYPE=TRA
SVLEN=.
STRAND=.
```

When `TRA` uses symbolic `<TRA>`:

- `CHR2` must contain the mate contig;
- `END` must contain the numeric mate position;
- `SVLEN` must be `.`;
- orientation may remain unknown and must not be invented solely to construct breakend notation.

This representation preserves caller outputs in which both breakpoint coordinates are known but breakend orientation is not reported. OctopuSV may compare or merge such TRA records using the known breakpoint coordinates; when orientation is available in both compared records, software may additionally use it as supporting evidence.

OctopuSV uses `TRA` when an event can be represented as a translocation from the available breakpoint evidence. Missing orientation alone does not invalidate an otherwise well-defined `TRA`.

---

## 12. Conversion to VCF

SVCF preserves relationships that conventional VCF cannot always represent directly. `svcf2vcf` therefore performs a defined, potentially lossy conversion.

For `TRA`, VCF export must preserve the remote breakpoint. Breakend-form `TRA` records retain the mate coordinate in `ALT`. Symbolic `<TRA>` records must retain `CHR2` and `END`, because those fields carry the remote breakpoint when orientation is unknown.

```bash
octopusv svcf2vcf -i input.svcf -o output.vcf
```

### 12.1 Caller-mode conversion

For a caller-mode record with exactly one evidence block, OctopuSV treats the block as a direct call rather than a synthesis step. The source genotype is preserved and source AD/LN may be retained in the VCF sample column.

For a caller-mode record with multiple evidence blocks, OctopuSV synthesizes one sample-level call using the same unique-caller consensus model described in Section 8.

For multi-evidence synthesis:

```text
AD=.,.
DP=.
```

because caller allele depths are not composable.

The VCF FORMAT for a synthesized multi-evidence caller record includes:

```text
GT:AD:DP:UC:UV:LN
```

where `UC`/`UV` make the synthesis interpretable.

A single-evidence caller record may use:

```text
GT:AD:DP:LN
```

VCF permits record-specific FORMAT layouts; SVCF-to-VCF export does not add empty UC/UV fields to records that were not synthesized.

For synthesized multi-evidence calls, `LN` follows the merged event length (`abs(INFO/SVLEN)` when available) rather than borrowing the length from an arbitrarily selected caller record.

### 12.2 Multi-mode conversion

SVCF 1.1 multi-mode sample columns export as:

```text
GT:AD:DP:UC:UV:LN
```

Synthesized calls normally have:

```text
AD=.,.
DP=.
```

### 12.3 Unobserved-sample export policy

For SVCF 1.1 multi files, `svcf2vcf` supports an explicit policy for true unobserved-event placeholders:

```text
--unobserved-sample-gt missing
--unobserved-sample-gt ref
```

The default is:

```text
missing
```

which exports a true unobserved placeholder as:

```text
./.
```

With:

```text
--unobserved-sample-gt ref
```

only true unobserved placeholders are exported as:

```text
0/0
```

Evidence-backed unresolved calls are never converted to `0/0` by this option.

When the policy is applicable, the output VCF records the selected behavior in a meta-information line:

```text
##OctopuSV_unobserved_sample_gt=missing
```

or:

```text
##OctopuSV_unobserved_sample_gt=ref
```

Using `ref` is an explicit operational choice for presence/absence cohort analysis. It is not equivalent to joint genotyping or a reference-confidence likelihood.

If `ref` is requested on an input that does not contain the SVCF 1.1 synthesized sample schema (`UC`/`UV`), OctopuSV warns that the option has no effect.

### 12.4 DP and SUPPORT are different quantities

`INFO/SUPPORT` is merged-event/representative support information. It must not be interpreted as a direct substitute for sample-level `FMT/DP`.

The original SVCF should be retained whenever caller-level evidence, source-specific read support, or the exact evidence relationships are needed.

---

## 13. Producer classification

SVCF identity follows the output data model.

| Workflow | Output classification |
|---|---|
| single-sample raw VCF -> `octopusv correct` | SVCF 1.1 caller |
| caller SVCFs -> `octopusv merge --mode caller` | SVCF 1.1 caller |
| per-sample caller SVCFs -> `octopusv merge --mode sample` | SVCF 1.1 multi |
| multi-sample raw VCF -> `octopusv correct` | legacy/unversioned SVCF |

The historical multi-sample `correct` layout is intentionally not assigned a third SVCF 1.1 mode. Its columns are biological samples, but its blocks still use caller evidence FORMAT rather than synthesized `UC`/`UV` sample calls.

**Reference OctopuSV workflow note (non-normative):** current OctopuSV merge preflight does not accept the legacy multi-sample `correct` representation as merge input. For cohort workflows, split the original multi-sample VCF into biological samples, run `octopusv correct` separately for each sample to obtain SVCF 1.1 caller files, and then combine those per-sample caller SVCFs with `octopusv merge --mode sample`.

---

## 14. Validation

The reference validation command is:

```bash
octopusv validate-svcf -i input.svcf
```

For SVCF 1.1, validation includes:

- supported `SVCFVersion`;
- explicit legal `OctopuSV_mode`;
- agreement between declared mode and `#CHROM` column shape;
- exact mode-specific FORMAT on every record;
- required leading header columns;
- required INFO keys;
- duplicate INFO-key detection;
- legal `SVTYPE`;
- valid `SUPPORT` syntax;
- caller/source/evidence column counts;
- `SOURCE_IDS` positional consistency when present;
- representable SVCF 1.1 `SOURCES` / `SOURCE_IDS` item syntax, including reserved-character and missing-value rules;
- sample-column count consistency in multi mode;
- structural-variant coordinate checks;
- BND breakend ALT/CHR2/END agreement;
- TRA breakpoint validation, including breakend ALT/CHR2/END agreement when breakend notation is used and `CHR2` plus numeric `END` when symbolic `<TRA>` is used;
- parseable `CO` values according to validator policy.

A versioned file that declares `caller` but uses the multi FORMAT is invalid. A versioned file that declares `multi` but uses the caller FORMAT is invalid. A versioned caller file with more than one trailing `#CHROM` column is invalid.

Unversioned legacy files remain readable through compatibility behavior where supported, but they are not automatically upgraded to the SVCF 1.1 contract.

By default, failure to find a parseable `CO` may be reported as a warning. It can be treated as an error with:

```bash
octopusv validate-svcf -i input.svcf --strict-co
```

---

## 15. Intermediate tools

A tool that rewrites a versioned SVCF file must preserve its explicit identity and must not emit a file that still claims SVCF 1.1 while violating the corresponding schema.

OctopuSV regression tests require structure-preserving intermediate operations such as filtering, subsetting, and contig normalization to retain the SVCF version/mode declarations and produce output that still passes `octopusv validate-svcf`.

For source-aware filtering, implementations must distinguish two explicit naming contexts. In merged records, `INFO/SOURCES` contains the user-facing source labels assigned to merge inputs (for example, labels derived from filenames or supplied with `--caller-names`). A direct single-evidence caller record may omit `SOURCES`; in that case its explicit FORMAT `SC` value identifies the caller/method software recorded by the source VCF. `SOURCES` labels and `SC` values are therefore not required to share a naming namespace. Readers must not substitute record-ID prefixes, filenames, or `#CHROM` labels when neither explicit representation is available.

---

## 16. Conformance and implementation guidance

A **conforming SVCF 1.1 writer** must emit files satisfying the required identity, header shape, FORMAT, and positional-binding rules in this specification. A **conforming SVCF 1.1 reader** must enforce the versioned identity and must not silently reinterpret an invalid 1.1 file through legacy heuristics.

The OctopuSV reference checker is:

```bash
octopusv validate-svcf -i input.svcf
```

Third-party implementations do not need to reproduce OctopuSV internals, but they should use the validator during development/interoperability testing.

Software that writes SVCF 1.1 should:

- declare `##SVCFVersion=1.1` and exactly one supported mode;
- use the exact FORMAT associated with that mode;
- keep caller evidence and sample synthesis as separate data models;
- preserve explicit source/evidence bindings;
- preserve duplicate source labels when multiple records from one caller support the same event;
- preserve positional `.` placeholders in `SOURCE_IDS`;
- reject `SOURCES` / non-missing `SOURCE_IDS` items containing reserved delimiters or whitespace rather than silently rewriting them;
- avoid reconstructing known source identity from filenames, record-ID prefixes, or column order;
- fail rather than silently downgrade when a versioned header/schema cannot be written correctly.

Software that reads SVCF 1.1 should:

- treat the declared version/mode as authoritative;
- reject unsupported versions rather than silently interpreting them as 1.1;
- reject a FORMAT that disagrees with the declared mode;
- parse colon-containing IDs/ALT values with a structure-aware SVCF block parser;
- distinguish raw caller evidence from synthesized sample calls;
- distinguish unobserved sample placeholders from evidence-backed absence or unresolved calls.

Files without `##SVCFVersion` may be handled through explicit legacy compatibility paths, but legacy inference must not override an explicit versioned identity.

The specification is the normative description of SVCF 1.1; the shared schema module and regression tests are the reference-implementation safeguards. When the contract changes, the specification, shared schema definition, writers, readers, validator, converters, and regression tests must be reviewed and updated together. A code change that intentionally changes a normative SVCF 1.1 rule must either remain backward-compatible with this specification or introduce a new SVCF version.
---

## Appendix A. Changes from legacy SVCF (non-normative)

This appendix summarizes migration-relevant differences between unversioned legacy SVCF and versioned SVCF 1.1. It is guidance for users and implementers; the normative requirements are defined in the main sections above.

1. **Explicit file identity.** SVCF 1.1 declares both `##SVCFVersion=1.1` and `##OctopuSV_mode=caller|multi`. Legacy files may omit the version declaration.
2. **Two explicit data models.** Caller mode is the evidence-preserving layer; multi mode is the synthesized biological-sample layer. Historical multi-sample `correct` output remains legacy rather than being assigned a third 1.1 mode.
3. **Mode-specific FORMAT schemas.** Caller mode uses `GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO`; multi mode uses `GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO`.
4. **Explicit source/evidence binding.** In caller-merge records, `SOURCES[i]`, `SOURCE_IDS[i]`, and evidence block `i` are positionally bound. Duplicate `SOURCES` values are legal because one caller may contribute multiple evidence records.
5. **Positional missing IDs.** `SOURCE_IDS=.` is a real positional placeholder and must not be dropped during parsing or rewriting.
6. **Representable positional atoms.** SVCF 1.1 defines no escaping for `SOURCES` / `SOURCE_IDS`; reserved delimiters and whitespace are rejected rather than silently rewritten.
7. **Unique-caller synthesis.** Multi-mode sample calls count each unique caller at most once. Carrier presence and zygosity are resolved separately.
8. **Partial carrier genotype.** `1/.` means at least one ALT allele is established while the second allele is unresolved; it is distinct from `./.`.
9. **Synthesized allele depth is unavailable.** Caller allele depths are not composable across callers, so synthesized sample calls use `AD=.,.` and VCF export uses `DP=.`.
10. **Unobserved sample policy is explicit.** Internal `UV=0` placeholders remain distinguishable in SVCF and VCF export records the selected `missing|ref` interpretation in the output header.
11. **Colon-containing values require structure-aware parsing.** Record IDs and ALT representations may contain colons; readers must not parse SVCF sample blocks with naive positional `split(":")` logic.
12. **Legacy compatibility is explicit, not authoritative.** Unversioned files may still be read through compatibility paths, but legacy inference must never override an explicit versioned SVCF identity.
13. **TRA orientation is optional.** `TRA` requires two known breakpoint coordinates, not necessarily known breakend orientation. Breakend-form `TRA` records encode the remote breakpoint and orientation in `ALT`; symbolic `<TRA>` records use `CHR2` and numeric `END` for the remote breakpoint and leave orientation unresolved. `BND` remains breakend-ALT only.
