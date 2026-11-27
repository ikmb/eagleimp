# Qref file specification

The qref format is used as a reference panel format in the EagleImp software for genotype imputation. It loads significantly faster than the standard vcf.gz or bcf variant formats.

Qref is a binary format and consists of:

### A) Header
**size:** 8 bytes

```
'Q' 'R' 'E' 'F' <v_major> <v_minor> 0 <chrom>
```

- `<v_major>` / `<v_minor>`: version number of the Qref format.
  The EagleImp software decides whether the current file is compatible with the current EagleImp version. In general, the desired major version has to match while the minor version may be smaller for backward compatibility.
- `<chrom>`: the chromosome number of this file (will be casted from byte to integer format), allowed numbers are 1 to 24 whereby 23 is used for chrX and 24 for chrY.

### B) Constants
**size:** 4 x 8 bytes

```
<Nref> <NrefHap> <Mref> <MrefMA>
```

- `<Nref>` Number of samples in the Qref
- `<NrefHap>` Number of samples in the Qref that are haploid
- `<Mref>` Number of variants in the Qref
- `<MrefMA>` Number of variants in the Qref that were multi-allelic (before splitting during conversion to Qref)

### C) Flags for haploid samples
**size:** Nref x 1 bytes or 0 bytes

```
<hapflag_0> <hapflag_1> ... <hapflag_(Nref-1)>
```

Each byte represents a flag if the corresponding sample is haploid (1) or not (0).
The flags are only present if NrefHap > 0 and NrefHap < Nref or chrom == 23 (chromosome X).
The flags are required to restore haploid samples whose haplotype data was encoded as homozygous diploid.

### D) Variant positions
**size:** Mref x 8 bytes

```
<pos_0> <pos_1> ... <pos_(Mref-1)>
```

Each 8 byte block is casted to int64 format and represents the *0-based* chromsomal position of a variant. Variants are ordered by position. Int64 type numbers are stored little-endian.

### E) Allele frequencies
**size:** Mref x 4 bytes

```
<af_0> <af_1> ... <af_(Mref-1)>
```

Each 4 byte block is casted to float (32bit) format and represents the allele frequency of the corresponding variant. The float types are stored little-endian.

### F) Alleles
**size:** 2x Mref C strings (zero delimited)

```
<ref_0> 0 <alt_0> 0 <ref_1> 0 <alt_1> 0 ... <ref_(Mref-1)> 0 <alt_(Mref-1)> 0
```

For each variant a pair of alleles (major and minor) is encoded as two C strings delimited by the zero-character. Major allele first.

### G) Variant IDs
**size:** Mref C strings

```
<id_0> 0 <id_1> 0 ... <id_(Mref-1)> 0
```

For each variant a zero-delimited C string denoting an identifier for the variant (third column in VCF format). The identifier may be empty in which case only the zero delimiter is present (empty C string).

### H) Flags for multi-allelic variants
**size:** Mref x 1 bit rounded to the next multiple of 64 bytes (512 bits)

```
byte 0: <maflag_07><maflag_06><maflag_05><maflag_04><maflag_03><maflag_02><maflag_01><maflag_00>
byte 1: <maflag_15><maflag_14><maflag_13><maflag_12><maflag_11><maflag_10><maflag_09><maflag_08>
...
```

For each variant a bit representing if the variant was split from a multi-allelic before conversion to Qref.
A byte is filled from LSB to MSB, i.e. less significant bits have lower variant indices than more significant bits.
Rounding to block sizes of 64 bytes is done for internal performance reasons.

### I) Haplotype data
**size:** Mref x (8 byte + bitvector of different lengths)

```
<encsize_0> <hapvec_0>
<encsize_1> <hapvec_1>
...
<encsize_(Mref-1)> <hapvec_(Mref-1)>
```

All haplotypes are encoded with a 0-bit representing the reference allele, a 1-bit representing the alternative allele.

The first 8 bytes (`<encsize>` as uint64, little-endian) before each haplotype vector indicate if it was runlength encoded or not.

* **Case `<encsize>` == 0:**
    The vector is **not** runlength encoded.

    `<hapvec>` is a vector of 2x Nref bits representing the haplotypes of each (diploid) sample at the corresponding variant. Note, that all haploid samples are encoded as homozygous diploid.

    Each byte in the bit vector is filled from LSB to MSB.

    The size of the bit vector is 2x Nref rounded to the next multiple of 64 bytes (512 bits). Rounding to block sizes of 64 bytes is done for internal performance reasons.

    No runlength encoding occurs if a sequence has too high variation and the runlength encoded sequence would occupy more space than the not encoded one.

    ```
    byte_0: <pat_3> <mat_3> <pat_2> <mat_2> <pat_1> <mat_1> <pat_0> <mat_0>
    byte_1: <pat_7> <mat_7> <pat_6> <mat_6> <pat_5> <mat_5> <pat_4> <mat_4>
    ...
    ```

* **Case `<encsize>` > 0:**
    The vector is runlength encoded.

    `<encsize>` determines the size of the encoded `<hapvec>` in bytes.
    `<hapvec>` is a byte vector where each byte or double-byte indicates the runlength of **alternating bits** (starting with `0`), depending on the most-significant bit (MSB) of the current byte:

    1. MSB is `0`:

        ```
        byte_0: 0 x x x x x x x
        ```
        The runlength of the current bit is `byte_0 & 0x7f`. (`&` is bit-wise AND.)

    2. MSB is `1`:

        ```
        byte_0: 1 x x x x x x x
        byte_1: y y y y y y y y
        ```
        The runlength of the current bit is `byte_1 << 7 | byte_0 & 0x7f`. (`<<` is arithmetic left-shift, `|` is bit-wise OR, `&` is bit-wise AND.)

    If the runlength for a bit is larger than the current maximum of 32767, the current bit is encoded with this maximum runlength, the next (alternated) bit is encoded with a runlength of 0, then the remaining runlength is encoded. If the remainder is still greater than 32767, the process continues until the remainder is less than 32767 where the encoding directly follows above rules.

**Examples:**

* No runlength encoding. The haplotype data is:

    ```
    0x0000000000000000 | 0x12 0x34 0x56 0x78 0x9a 0xbc 0xde 0xff 0x1c 0x33 ...
    ```
    The first 8 bytes are zero (`0x0000000000000000`), so the next 2xNref bits directly represent the haplotype data.

* Runlength encoded sequence:

    ```
    0x0000000000000008 | 0x0a 0x03 0xaa 0x02 0xff 0xff 0x00 0x20
    ```
    The size of the encoded vector is 8 bytes (`0x0000000000000008`), so the following 8 bytes indicate a runlength encoded sequence.
    All runlength encoded sequences start with the 0-bit. (To start a runlength encoded sequence with 1, the first runlength has to be 0.)

    `0x0a`: The MSB is zero, so the sequence starts with 10 subsequent zeros (`0x0a & 0x7f = 0x0a = 10`).

    `0x03`: The MSB is again zero, so the next bits are 3 subsequent ones (`0x03 & 0x7f = 0x03 = 3`).

    `0xaa 0x02`: The MSB of the current byte is one, so the next byte is also part of the runlength. So, the next bits are 298 subsequent zeros: `0x02 << 7 | (0xaa & 0x7f) = 0x0100 | 0x002a = 256 + 42 = 298`.

    `0xff 0xff`: As before, both bytes build the runlength, which results in the maximum of 32767 here (`0xff << 7 | (0xff & 0x7f) = 32767`). Thus the sequence continues with 32767 subsequent ones.

    `0x00`: The runlength is zero, so no zeros follow the 32767 ones from before.

    `0x20`: The MSB is zero, so `0x20 = 32` ones are attached to the 32767 ones from before (as the runlength for zeros was zero before).
    In fact, the last four bytes indicate an example how to encode a runlength exceeding the maximum of 32767 (here 32799).

