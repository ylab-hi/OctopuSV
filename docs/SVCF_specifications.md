# SVCF: A VCF-Based Format for Structural Variant Processing and Integration

## 1. Introduction

SVCF is a VCF-based format for structural variant normalization, integration, provenance tracking, and downstream analysis. It is defined by this specification and implemented by OctopuSV.

SVCF keeps the line-oriented text structure and core columns of the [Variant Call Format (VCF) Version 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf). It adds a fixed evidence schema and a small set of conventions needed to preserve structural variant records from different callers and samples.

SVCF does not replace or revise the official VCF specification. A general VCF program may be able to read parts of an SVCF file, but it is not expected to understand every SVCF-specific field or layout. In particular, caller-merge SVCF records may contain a variable number of evidence columns.

Software other than OctopuSV may read or write SVCF when it follows the rules in this document.

The recommended file extension is:

```text
.svcf
```

Before using an SVCF file with general VCF tools such as bcftools or vcftools, convert it to VCF:

```bash
octopusv svcf2vcf -i input.svcf -o output.vcf
```

This revision describes the SVCF contract implemented by OctopuSV v0.4.1.

The words **must**, **should**, and **may** are used in their ordinary specification sense:

- **must** means that the rule is required;
- **should** means that the rule is recommended unless there is a clear reason not to follow it;
- **may** means that the item is optional.

## 2. File format

An SVCF file contains:

1. meta-information lines beginning with `##`;
2. one header line beginning with `#CHROM`;
3. tab-delimited data lines.

Missing values are written as a single dot (`.`), following VCF convention.

### 2.1 Meta-information lines

The first meta-information line should be:

```text
##fileformat=VCFv4.2
```

OctopuSV normally writes:

```text
##source=OctopuSV
```

It may also write:

```text
##fileDate=...
##OctopuSV_WARNING=This is SVCF format. Use 'octopusv svcf2vcf' to change back to standard VCF format before bcftools/vcftools
```

A sample/multi SVCF is identified by:

```text
##OctopuSV_mode=multi
```

The header should define the contigs, symbolic ALT alleles, INFO fields, FILTER values, and FORMAT fields used in the file.

OctopuSV may preserve metadata definitions from an input VCF or SVCF. When an input header already defines the same INFO, FORMAT, FILTER, or ALT ID, the current OctopuSV header writer keeps the input definition. Preserved metadata does not change the SVCF record rules in this specification.

### 2.2 Header line

The first nine columns are fixed and must appear in this order:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
```

At least one evidence column must follow `FORMAT`.

The complete header therefore has the form:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  evidence-column-1  ...
```

The first nine columns describe the normalized or merged SVCF record. The columns after `FORMAT` preserve evidence from source caller records or samples.

The exact meaning and number of evidence columns depend on the SVCF mode described in Section 3.

### 2.3 Data lines

Each data line must:

- be tab-delimited;
- contain at least ten columns;
- use the nine fixed columns in the required order;
- use the SVCF FORMAT string defined in Section 6;
- contain the number of evidence columns required by its mode.

A compact caller-merge example is shown below:

```text
1	10889	svim.INS.3	T	TCCAGGGGAGGAGGCGTGGCACAGGCGCAGAGACACATGCTAGCGCGC	34	PASS	SVTYPE=INS;END=10936;SVLEN=47;CHR2=1;SUPPORT=28;SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.;SOURCES=svim,pbsv;SOURCE_IDS=svim.INS.3,pbsv.INS.1	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	1/1:1,28:47:.:34:INS:svim.INS.3:SVIM-v2.0.0:T:TCCAGGGGAGGAGGCGTGGCACAGGCGCAGAGACACATGCTAGCGCGC:1_10889-1_10936	1/1:1,25:47:.:.:INS:pbsv.INS.1:pbsv:G:GAGGAGGCGTGGCACAGGCGCAGAGACACATGCTAGCGCGCCCAGGGG:1_10896-1_10943
```

## 3. SVCF modes

SVCF has three modes.

| Mode | Mode marker | `SOURCES` | Evidence columns |
|---|---|---|---|
| Single-caller | No `##OctopuSV_mode=multi` | Must be absent | Exactly one |
| Caller-merge | No `##OctopuSV_mode=multi` | Must be present and non-empty on every record | One column for each entry in `SOURCES`; the number may vary by record |
| Sample/multi | `##OctopuSV_mode=multi` | Must be present and non-empty on every record | A fixed number matching the names after `FORMAT` in the `#CHROM` header |

A file without the multi-sample marker must not mix single-caller records and caller-merge records. In a no-marker file, either all records contain `SOURCES` or none of them do.

### 3.1 Single-caller mode

Single-caller SVCF is normally produced by `octopusv correct`.

It has one evidence column for each record. The column name normally comes from the input sample name or an OctopuSV fallback name.

A single-caller record must not contain a non-empty `SOURCES` field.

### 3.2 Caller-merge mode

Caller-merge SVCF is produced when records from multiple callers are merged.

Every record must contain `SOURCES`. The number of evidence columns on that record must equal the number of comma-separated entries in `SOURCES`.

The number of evidence columns may differ between records. A record supported by two callers has two evidence columns; a record supported by three callers has three.

This variable-width layout is part of caller-merge SVCF. It is not a conventional VCF sample layout.

### 3.3 Sample/multi mode

Sample/multi SVCF contains:

```text
##OctopuSV_mode=multi
```

The names after `FORMAT` in the `#CHROM` header define the evidence-column order for the whole file.

Every data line must contain exactly that number of evidence columns.

## 4. Fixed fields

The fixed fields follow the general VCF 4.2 field model, with the SVCF rules below.

| Field | SVCF rule |
|---|---|
| `CHROM` | Contig containing the event start or first breakpoint. It must not be empty. |
| `POS` | Positive, 1-based coordinate of the event start or first breakpoint. |
| `ID` | Identifier of the representative SVCF record. |
| `REF` | Reference allele of the representative record. |
| `ALT` | Alternate allele of the representative record. For `TRA` and `BND`, it must use a valid VCF breakend bracket form. |
| `QUAL` | Quality value of the representative record, or `.` when unavailable. |
| `FILTER` | Filter status of the representative record. `PASS` means that the record passed the filters applied by its source workflow. |
| `INFO` | Event-level annotations and provenance fields. |
| `FORMAT` | The fixed SVCF evidence schema. |

The representative record is the normalized or selected record written in the first nine columns. Its values may differ from the original values kept in the evidence columns.

## 5. INFO fields

SVCF uses the following INFO definitions:

```text
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="The software used to identify the SV">
##INFO=<ID=RTID,Number=1,Type=String,Description="Associated ID for reciprocal translocations if available">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Source caller/sample labels supporting this merged SV record">
##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original IDs of merged SVs from different callers or samples">
```

Every SVCF record must contain these keys:

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

The key must be present even when its value is `.`.

`SOURCES` is mode-dependent:

- it must be absent or empty in single-caller mode;
- it must be present and non-empty in caller-merge and sample/multi modes.

`SOURCE_IDS` is written by the current OctopuSV merge implementation to preserve original record IDs. It is not used by the current validator to determine the SVCF mode.

### 5.1 INFO field meanings

| Field | Meaning |
|---|---|
| `SVTYPE` | Structural variant type. Allowed values are `DEL`, `DUP`, `INV`, `INS`, `TRA`, and `BND`. |
| `CHR2` | Contig containing the second coordinate or mate breakpoint. |
| `END` | Second coordinate used by the SVCF record. Its meaning depends on `SVTYPE`. |
| `SVLEN` | Positive event length when applicable. `TRA` and `BND` use `.`. |
| `SUPPORT` | Read-support value of the representative record. It is not the number of callers or samples. |
| `SVMETHOD` | Method that produced the current SVCF record, normally `OctopuSV`. |
| `RTID` | Related or reciprocal record ID when available; otherwise `.`. |
| `AF` | Allele frequency when available; otherwise `.`. |
| `STRAND` | Strand or orientation value retained by OctopuSV; otherwise `.`. |
| `RNAMES` | Supporting read names when available; otherwise `.`. |
| `SOURCES` | Comma-separated labels of the sources supporting a merged record. |
| `SOURCE_IDS` | Comma-separated original record IDs written for source-level traceability. |

`SUPPORT` may be `.` or a non-negative integer.

For example:

```text
SOURCES=svim,pbsv
SOURCE_IDS=svim.INS.3,pbsv.INS.1
```

In caller-merge mode, the source order is also the evidence-column order:

```text
SOURCES item 1  <-> evidence column 1
SOURCES item 2  <-> evidence column 2
```

## 6. FORMAT and evidence columns

The SVCF FORMAT string is fixed:

```text
GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO
```

The order must not change.

The corresponding FORMAT definitions are:

```text
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">
##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV (e.g., +, -, -+, ++)">
##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Quality value">
##FORMAT=<ID=TY,Number=1,Type=String,Description="Type of SV (e.g., TRA, DEL, INS)">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Unique identifier for the SV">
##FORMAT=<ID=SC,Number=1,Type=String,Description="Source from which SV was identified">
##FORMAT=<ID=REF,Number=1,Type=String,Description="Reference allele sequence">
##FORMAT=<ID=ALT,Number=1,Type=String,Description="Alternate allele sequence">
##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate information of the SV">
```

Each evidence column describes one source record or one sample entry.

| Field | Meaning |
|---|---|
| `GT` | Genotype retained from the source record. |
| `AD` | Reference and alternate allele depths when available. |
| `LN` | Absolute source-event length when available. |
| `ST` | Source strand or orientation value when available. |
| `QV` | Source quality value when available. |
| `TY` | Source structural variant type. |
| `ID` | Original source record ID. |
| `SC` | Source caller, caller version, file label, or another source name written by OctopuSV. |
| `REF` | Original source REF value. |
| `ALT` | Original source ALT value, including the original BND expression when applicable. |
| `CO` | Source coordinates in the form described below. |

`CO` must remain the last FORMAT field.

The `CO` value has the form:

```text
startChrom_startPos-endChrom_endPos
```

Examples:

```text
1_10889-1_10936
1_3845267-hs37d5_32469995
```

Contig names may contain underscores or hyphens. A parser should separate each coordinate from the rightmost underscore before the numeric position.

Source IDs and BND ALT strings may contain colons. Software must not assume that every colon inside an evidence column is a FORMAT separator. The current OctopuSV parser keeps `CO` as the final field and uses bounded or right-side parsing where needed.

For non-`TRA` and non-`BND` records, the current SVCF parser first looks for a usable `CO` value in the evidence columns. If no usable `CO` is found, it falls back to `CHROM`, `POS`, `CHR2`, `END`, and, when needed, `SVLEN`.

For `TRA` and `BND`, the parser first reads the mate coordinate from `ALT` and falls back to `CHR2` and `END` when needed.

## 7. Structural variant representation

The legal SVCF structural variant types are:

```text
DEL
DUP
INV
INS
TRA
BND
```

### 7.1 DEL, DUP, and INV

For `DEL`, `DUP`, and `INV`:

- `END` must be a numeric value;
- `END` must be greater than or equal to `POS`;
- `CHR2` is normally the same as `CHROM`;
- `SVLEN` is a positive event length when available.

Current OctopuSV output uses the absolute length rather than a negative deletion length.

### 7.2 INS

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

### 7.3 TRA and BND

`TRA` and `BND` use VCF breakend notation in `ALT`.

The accepted forms are:

```text
t[chr:pos[
t]chr:pos]
[chr:pos[t
]chr:pos]t
```

where `t` is sequence placed before or after the breakend expression.

For `TRA` and `BND`:

- `CHR2` must contain the mate contig;
- `END` must contain the numeric mate position;
- the mate contig and position in `ALT` must agree with `CHR2` and `END`;
- `SVLEN` must be `.`.

Example:

```text
CHROM=1
POS=3845267
ALT=C[hs37d5:32469995[
CHR2=hs37d5
END=32469995
SVTYPE=TRA
SVLEN=.
```

OctopuSV uses `TRA` when the breakend can be represented as a translocation with sufficient confidence.

`BND` is retained when OctopuSV cannot safely convert a breakend to another supported SV type. A retained `BND` record is valid SVCF and is not considered a conversion failure.

## 8. Provenance

SVCF separates the representative event from the source evidence.

The fixed fields and INFO describe the representative normalized or merged event. Evidence columns retain source-specific values such as genotype, quality, coordinates, caller ID, REF, and ALT.

The representative values do not have to match every source record exactly. For example, two insertion records may be merged while retaining slightly different source positions, lengths, or inserted sequences in their evidence columns.

In caller-merge mode:

- `SOURCES` identifies the source columns;
- `SOURCE_IDS` preserves original source record IDs;
- the evidence columns follow the `SOURCES` order;
- each source event remains traceable through its evidence column.

In sample/multi mode, the `#CHROM` header defines a fixed global order for the evidence columns.

## 9. Missing values

A missing value must be written as `.` rather than `NA`.

Common examples include:

```text
QUAL=.
SUPPORT=.
SVLEN=.
GT=./.
AD=.,.
QV=.
CO=.
```

`./.` represents a missing genotype reported in an evidence column.

A dot in another FORMAT field means that the corresponding source value is unavailable.

A missing value does not remove the field from the fixed SVCF FORMAT layout.

## 10. Conversion to VCF

SVCF is intended to preserve normalized structural variant data and source evidence. VCF is the exchange format for software that does not implement SVCF.

Convert an SVCF file with:

```bash
octopusv svcf2vcf -i input.svcf -o output.vcf
```

Conversion may be lossy because a conventional VCF does not preserve the full caller-merge evidence layout.

In the current caller-merge conversion used by the OctopuSV regression tests:

- the representative event remains one VCF record;
- `SOURCES` and `SOURCE_IDS` remain in INFO;
- the variable SVCF evidence columns are reduced to a conventional sample column;
- the output FORMAT is `GT:AD:DP:LN`;
- `DP` is calculated from the available allele depths;
- insertion `END` is changed from the internal SVCF endpoint to `POS`.

The original SVCF should be retained whenever source-level provenance is needed.

## 11. Validation

The reference validation command is:

```bash
octopusv validate-svcf -i input.svcf
```

The current validator checks:

- the required leading header columns;
- legal SVCF mode;
- the fixed FORMAT string and order;
- required INFO keys;
- legal `SVTYPE`;
- valid `SUPPORT`;
- mode-specific evidence-column counts;
- numeric `END` values for `DEL`, `DUP`, and `INV`;
- BND syntax for `TRA` and `BND`;
- agreement among BND `ALT`, `CHR2`, and `END`;
- the presence of a parseable `CO` value.

By default, failure to find a parseable `CO` is reported as a warning. It can be treated as an error with:

```bash
octopusv validate-svcf -i input.svcf --strict-co
```

The exit codes are:

```text
0  PASS or PASS_WITH_WARNINGS
1  validation failed
2  file unreadable or empty
```

## 12. Implementing SVCF

Software that writes SVCF should follow the required columns, INFO fields,
FORMAT layout, SVCF modes, and structural-variant coordinate rules described
in this document.

Software that reads SVCF should recognize the three SVCF modes, preserve
source evidence and provenance, and correctly handle TRA and BND breakend
records.

`octopusv validate-svcf` can be used to check whether a file follows the
current SVCF rules implemented by OctopuSV.

When the SVCF format changes, the specification, writers, parsers, validator,
converter, and regression tests should be updated together.

