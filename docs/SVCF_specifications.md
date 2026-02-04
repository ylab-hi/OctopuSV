# SVCF: A Customized VCF Format for Structural Variants

## 1. Introduction

SVCF is simply a customized version of the [Variant Call Format (VCF) Version 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) adapted for working with structural variants (SVs). This is not a new format standard, but merely our internal implementation with minor adjustments to better handle structural variant data.

This document describes how we use the standard VCF format with some custom fields and conventions to support structural variant analysis. The customization primarily adds specific fields for representing structural variants while maintaining full compatibility with standard VCF tools. 

These adaptations were created to improve integration with the SVoctopus software and facilitate operations such as merging, statistical analysis, and benchmarking of structural variants.

## 2. File Structure

Our customized VCF follows the standard VCF structure:

- **Meta-information lines**: Start with "##" and provide metadata about the file (e.g., format version, reference used).
- **Header line**: Starts with a single "#" and lists all the fields that will appear in the body of the file.
- **Data lines**: Each data line contains information about a single variant and follows the columns specified in the header.

## 3. An Example

```plaintext
##fileformat=VCFv4.2
##fileDate=2023-11-07|11:32:51AM|CDT|-0500
##source=octopusV
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr17,length=83257441>
##contig=<ID=GL000008.2,length=209709>
##contig=<ID=GL000194.1,length=191469>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=BND,Description="Breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="The software used to identify the SV">
##INFO=<ID=RTID,Number=1,Type=String,Description="Associated ID for reciprocal translocations if available">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">
##FILTER=<ID=PASS,Description="All filters passed, variant is most likely true">
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample
GL000008.2	105803	svim.BND.1	N	N[chr13:67876737[	12	PASS	SVTYPE=TRA;END=67876737;SVLEN=.;CHR2=chr13;SUPPORT=10;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/1:1,3:.:+:12:TRA:svim.BND.1:svim.vcf:N:N[chr13:67876737[:GL000008.2_105803-chr13_67876737
GL000008.2	105804	svim.BND.2	N	]chr13:67886734]N	12	PASS	SVTYPE=TRA;END=67886734;SVLEN=.;CHR2=chr13;SUPPORT=10;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/0:2,8:.:-:12:TRA:svim.BND.2:svim.vcf:N:]chr13:67886734]N:GL000008.2_105804-chr13_67886734
GL000194.1	99821	svim.BND.3	N	N]chrX:119789373]	4	PASS	SVTYPE=TRA;END=119789373;SVLEN=.;CHR2=chrX;SUPPORT=4;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	1/1:0,4:.:+-:4:TRA:svim.BND.3:svim.vcf:N:N]chrX:119789373]:GL000194.1_99821-chrX_119789373
GL000008.2	13306	svim.DEL.1	TTCTCCCACTTTTTGATGGGGTTGTTTTTTTCTTGTAAATTTGTTTG	T	1	PASS	SVTYPE=DEL;END=13352;SVLEN=46;CHR2=GL000008.2;SUPPORT=1;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/0:1,0:46:+-:1:DEL:svim.DEL.1:svim.vcf:TTCTCCCACTTTTTGATGGGGTTGTTTTTTTCTTGTAAATTTGTTTG:T:GL000008.2_13306-GL000008.2_13352
GL000008.2	38637	svim.INS.2	G	GCTCGGGCCAAACTAGATGCTACCTTAATACACGTCTCACAA	1	PASS	SVTYPE=INS;END=38637;SVLEN=41;CHR2=GL000008.2;SUPPORT=1;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/1:1,0:41:+-:1:INS:svim.INS.2:svim.vcf:G:GCTCGGGCCAAACTAGATGCTACCTTAATACACGTCTCACAA:GL000008.2_38637-GL000008.2_38678
chr1    16725233    svim.BND.1076   N   <INV>   7   PASS    SVTYPE=INV;END=16725253;SVLEN=20;CHR2=chr1;SUPPORT=6;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+ GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO  1/0:5,3:20:+:7:INV:svim.BND.1076:svim.vcf:N:[chr1:16725253[N:chr1_16725233-chr1_16725253
chr7	117502838	svim.BND.8880	N	<DUP>	1	PASS	SVTYPE=DUP;END=131967043;SVLEN=14464205;CHR2=chr7;SUPPORT=1;SVMETHOD=octopusV;RTID=.;AF=.;STRAND=+-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/1:0,1:14464205:+-:1:DUP:svim.BND.8880:svim.vcf:N:]chr7:131967043]N:chr7_117502838-chr7_131967043
chr8	124536858	svim.BND.7772	N	N[chr17:39911421[	1	PASS	SVTYPE=TRA;END=39911421;SVLEN=.;CHR2=chr17;SUPPORT=1;SVMETHOD=octopusV;RTID=svim.BND.72;AF=.;STRAND=+-	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/0:2,2:.:+-:1:TRA:svim.BND.7772:svim.vcf:N:N[chr17:39911421[:chr8_124536858-chr17_39911421
chr17	39911421	svim.BND.72	N	N[chr8:124536858[	7	PASS	SVTYPE=TRA;END=124536858;SVLEN=.;CHR2=chr17;SUPPORT=6;SVMETHOD=octopusV;RTID=svim.BND.7772;AF=.;STRAND=+	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/1:2,4:.:+:7:TRA:svim.BND.72:svim.vcf:N:N[chr8:124536858[:chr17_39911421-chr8_124536858
```

## 4. Meta-Information Lines

Meta-information lines contain various keywords prefixed with double hash marks "##" that define how the data lines should be interpreted. Each meta-information line is a directive typically defining the file format version, the reference genome, contigs, and information specific to the format or analysis tools.

## 5. Header Line

The header line starts with a single "#" and provides the names of the columns in the data lines:

```plaintext
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample
```

## 6. Data Lines

Each data line corresponds to a structural variant and includes the following mandatory columns:

- **CHROM**: The chromosome on which the variant occurs.
- **POS**: The position of the variant on the chromosome.
- **ID**: A semi-colon-separated list of unique identifiers for the variant.
- **REF**: The reference allele at the position of the variant.
- **ALT**: The alternate allele, specifying the variant found in a sample.
- **QUAL**: A quality score for the variant call.
- **FILTER**: Indicates if the variant passed filtering.
- **INFO**: Additional information about the variant.
- **FORMAT**: The format of the sample-specific data.
- **Sample**: The sample-specific data for the variant.

## 7. INFO Field

The INFO field contains several subfields separated by semicolons, including:

- **SVTYPE**: Type of structural variant (DUP, INV, TRA, DEL, INS).
- **CHR2**: Chromosome for the end position of the SV.
- **END**: End position of the variant; set to 0 for TRA.
- **SVLEN**: Length of the SV; calculated as start-end for DEL, end-start for DUP and INV; set to '.' for TRA.
- **SUPPORT**: Number of pieces of evidence supporting the SV.
- **SVMETHOD**: Software used to identify the SV.
- **RTID**: Identifier for reciprocal translocations, if applicable.
- **AF**: Allele frequency of the SV.
- **STRAND**: Strand orientation of the SV.

Each tag in the INFO field would be used to annotate a structural variant in the SVCF file, providing a quick reference for key details of each event.

## 8. FORMAT Field

The FORMAT field includes:

- **GT**: Genotype, represented as allele indices separated by a slash (`/`) for diploid genomes or a pipe (`|`) for phased genotypes.
- **AD**: Allelic depths for the ref and alt alleles in the order listed.
- **LN**: The length of the SV event.
- **ST**: Strands for the SV event, if applicable.
- **QV**: Quality Value for the SV event.
- **TY**: The type of SV event.
- **ID**: A unique identifier for the SV event.
- **SC**: Source caller, often the VCF file name or a user-defined string.
- **REF**: The reference bases for the SV event.
- **ALT**: The alternate bases for the SV event (or `.` if none). It can also be BND pattern.
- **CO**: Coordinates for the SV event indicating the affected region.

## 9. Conventions

### SVTYPE Convention

Within our customized VCF format, the `SVTYPE` field is restricted to five canonical types of structural variants:

- `DUP`: Duplication, a segment of DNA that is duplicated on the genome.
- `INV`: Inversion, a segment of DNA that is inverted in the genome.
- `TRA`: Translocation, a segment of DNA that has been excised from one part of the genome and relocated to another position.
- `DEL`: Deletion, a segment of DNA that is missing from the genome.
- `INS`: Insertion, a segment of DNA that is added to the genome.

### Special Rules for `TRA`

For translocations (`TRA`), `SVLEN` is represented by a dot (`.`) to avoid erroneous length calculations that could potentially reduce the overall length of the reference sequence.

```plaintext
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample
chrX	12345678	svim.TRA.91011	A	<TRA>	60	PASS	SVTYPE=TRA;END=3456;SVLEN=.;SUPPORT=5;SVMETHOD=octopusV	GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO	0/1:20,30::.::60:TRA:svim.TRA.91011:example.vcf:A:<TRA>:chrX_12345678-chrY_23456789
```

### Missing Values

In situations where data is not available, a dot (`.`) is used to represent the absence of data, not `NA`. This convention aligns with the VCF specification for missing information.

### FORMAT Field: SC

The `SC` subfield within the FORMAT field stands for "source" and is by default populated with the name of the VCF file from which the SV data is derived. Users have the option to customize this subfield to facilitate tracing the origin of the SV (such as the sample, file, platform, etc.) in subsequent analyses. This is distinct from the `SVMETHOD` field, which indicates the software (typically `octopusV`) that generated the current file.

### SVLEN Calculation

The `SVLEN` field reflects the length of the structural variation and is calculated as follows:

- For `DEL`: `SVLEN` is the positive value calculated by subtracting the `POS` (start position) from the `END`, indicating a deletion from the reference genome. The length of the deletion is the number of base pairs lost.
- For `DUP` and `INV`: `SVLEN` is the positive value calculated by subtracting the `POS` (start position) from the `END`, indicating an increase in the length of the genome. The length is the number of base pairs gained.
- For `TRA`: `SVLEN` is not applicable as the length does not change in a simple translocation event. The field is represented by a period (`.`) to indicate this.

Here is the corrected representation:

- For `DEL`:

  ```
  SVLEN = END - POS
  ```

- For `DUP` and `INV`:

  ```
  SVLEN = END - POS
  ```

- For `TRA`:

  ```
  SVLEN = .
  ```

- For `INS`:
  ```
  SVLEN = Fixed number
  ```