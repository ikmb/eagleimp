#
#    Copyright (C) 2018-2021 by Lars Wienbrandt,
#    Institute of Clinical Molecular Biology, Kiel University
#
#    This file is part of EagleImp.
#
#    EagleImp is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    EagleImp is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with EagleImp. If not, see <https://www.gnu.org/licenses/>.
#

BEGIN {
  chunkline=0
  summaryline=0
  nchunks=0
  nchunklines=0
  ntsamples=0
  ntvars=0
  nrsamples=0
  nrvars=0
  nimpout=0
  physrange=0.0
  genrange=0.0
  refaltsflip=0
  refalt=0
  sflip=0
  dropnotfound=0
  dropmatgt=0
  dropmaref=0
  dropmono=0
  dropmismatch=0
  nunphased=0
  nrefonly=0
  exmarefonly=0
  miss=0
  impconfl=0
}

# needs to be parsed before the general match
/Chunk/{ # note, that this also matches the "Chunks" field!
  if ($2 != "Chunks:") { # ignore the total-number-of-chunks field
    chunkline = FNR
  }
}

# needs to be parsed before the general match
/Summary/{
  # store the line number of the summary
  summaryline = FNR
}

# general
{
  # reset the line number of "Summary" and "Chunk" whenever a new file is read,
  # and set chromosome
  if (FNR==1) {
    chunkline = 0
    summaryline = 0
    # chromosome counter
    nchrom++
    # get the basename of the file and split into fields separated by '.'
    split(FILENAME,fnfields,"/")
    split(fnfields[length(fnfields)],fnfieldsx,".")
    chr=fnfieldsx[1]
    #gsub("chr","",chr)
    chrom[nchrom]=chr
  }

  # generally, copy the input if we are in a chunk description
  if (chunkline > 0 && summaryline == 0) {
    # exception: manipulate the chunk index
    if ($1=="Index:") {
      nchunks++ # one-based index
      chunklines[nchunklines] = "    Index: " chrom[nchrom] "/" $2
    } else chunklines[nchunklines] = $0
    nchunklines++
  }
}

#/Chromosome/ we do not parse the chromosome name from here as it may not be correctly reflected (e.g. X_PAR1 is just nominated as X)
/Target samples/{ # they should be the same in all (max) three input files!
  # we are taking the maximum though
  if (ntsamples < $3) ntsamples=$3
}
/Target variants/{
  ntvars+=$3
}
/Reference samples/{ # they should be the same in all (max) three input files!
  # we are taking the maximum though
  if (nrsamples < $3) nrsamples=$3
}
/Reference variants/{
  nrvars+=$3
}
/Common variants/{
  if (summaryline>0) {
    nphased[nchrom]=$3
  }
}
/Exclusive reference variants/{
  if (summaryline>0) {
    nimputed[nchrom]=$4
  }
}
/Imputation output variants/{
  if (summaryline>0) {
    nimpout+=$4
  }
}
/Physical distance range/{
  if (summaryline>0) {
    physrange+=$4
  }
}
/Genetic distance range/{
  if (summaryline>0) {
    genrange+=$4
  }
}
/REF ALT swap strand flips/{
  if (summaryline>0) {
    refaltsflip+=$6
  }
}
/REF ALT swaps/{
  if (summaryline>0) {
    refalt+=$4
  }
}
/Strand flips/{
  if (summaryline>0) {
    sflip+=$3
  }
}
/Dropped target only variants/{
  if (summaryline>0) {
    dropnotfound+=$5
  }
}
/Dropped multi-allelic target variants/{
  if (summaryline>0) {
    dropmatgt+=$5
  }
}
/Dropped multi-allelic reference variants/{
  if (summaryline>0) {
    dropmaref+=$5
  }
}
# /Common reference variants from multi-allelic splits/ not in webservice
/Dropped target variants monomorphic in reference/{
  if (summaryline>0) {
    dropmono+=$7
  }
}
/Dropped allele mismatched variants/{
  if (summaryline>0) {
    dropmismatch+=$5
  }
}
# /User excluded variants/ not in webservice
/Unphased variants after phasing/{
  if (summaryline>0) {
    nunphased+=$5
  }
}
/Reference-only variants/{
  if (summaryline>0) {
    nrefonly+=$3
  }
}
/Excluded reference-only multi-allelic variants/{
  if (summaryline>0) {
    exmarefonly+=$5
  }
}
# /Reference-only bi-allelic variants from multi-allelic splits/ not in webservice
# missing/unphased genotypes in phasing/imputation references are ignored here
/Missing rate in target genotypes/{
  if (summaryline>0) {
    miss+=$6*nphased[nchrom]
  }
}
/Imputation conflict rate/{
  if (summaryline>0) {
    impconfl+=$4*nimputed[nchrom]
  }
}
END {
print "---"
print "- General information:"
printf "    Chromosome: "
printf chrom[1]
for (c=2; c<=nchrom; c++) printf " " chrom[c]
print ""
print "    Target samples: " ntsamples
print "    Target variants: " ntvars
print "    Reference samples: " nrsamples
print "    Reference variants: " nrvars
print ""
print "- Chunks: " nchunks
print ""
for (n=0; n<nchunklines; n++) print chunklines[n]

print "- Summary:"
np=0; for (n=1; n<=nchrom; n++) np += nphased[n]
print "    Common variants: " np
ni=0; for (n=1; n<=nchrom; n++) ni += nimputed[n]
if (ni>0) {
  print "    Exclusive reference variants: " ni
}
if (nimpout>0)  print "    Imputation output variants: " nimpout
print "    Physical distance range: " physrange " bp"
print "    Genetic distance range: " genrange " cM"
printf "    Average SNPs per cM in target: %.0f\n", (np/genrange)+0.5
if (refaltsflip>0) print "    REF ALT swap strand flips: " refaltsflip
if (refalt>0)      print "    REF ALT swaps: " refalt
if (sflip>0) print "    Strand flips: " sflip
print "    Dropped target only variants: " dropnotfound
if (dropmatgt>0) print "    Dropped multi-allelic target variants: " dropmatgt
if (dropmaref>0) print "    Dropped multi-allelic reference variants: " dropmaref
if (dropmono>0)  print "    Dropped target variants monomorphic in reference: " dropmono
if (dropmismatch>0) print "    Dropped allele mismatched variants: " dropmismatch
if (nunphased>0) print "    Unphased variants after phasing: " nunphased
print "    Reference-only variants: " nrefonly
print "    Excluded reference-only multi-allelic variants: " exmarefonly
printf "    Missing rate in target genotypes: %.8f\n", (miss/np)
if (ni>0) printf "    Imputation conflict rate: %.8f\n", (impconfl/ni)
print ""
}
