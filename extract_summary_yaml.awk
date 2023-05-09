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
  summaryline=0
  nchrom=0
  ntsamples_min=1000000000
  ntsamples_max=0
  ntvars=0
  nrsamples_min=1000000000
  nrsamples_max=0
  nrvars=0
  nphased=0
  nimputed=0
  nimpout=0
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
  missrate_min=1.0
  missrate_max=0.0
  impconfl_min=1.0
  impconfl_max=0.0
}

{
  # reset the line number of "Summary" whenever a new file is read
  if (FNR==1) summaryline = 0
}

/Summary/{
  # store the line number of the summary
  summaryline = FNR
  # get chromosome from the filename whenever the summary line is read
  nchrom++
  # get the basename of the file and split into fields separated by '.'
  split(FILENAME,fnfields,"/")
  split(fnfields[length(fnfields)],fnfieldsx,".")
  chr=fnfieldsx[1]
  gsub("chr","",chr)
  chrom[nchrom]=chr
}
/Target samples/{
  if (ntsamples_min>$3) ntsamples_min=$3
  if (ntsamples_max<$3) ntsamples_max=$3
}
/Target variants/{
  ntvars+=$3
}
/Reference samples/{
  if (nrsamples_min>$3) nrsamples_min=$3
  if (nrsamples_max<$3) nrsamples_max=$3
}
/Reference variants/{
  nrvars+=$3
}
/Common variants/{
  if (summaryline>0) {
    nphased+=$3
  }
}
/Exclusive reference variants/{
  if (summaryline>0) {
    nimputed+=$4
  }
}
/Imputation output variants/{
  if (summaryline>0) {
    nimpout+=$4
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
  mrate=$6
  if (missrate_min>mrate) missrate_min=mrate
  if (missrate_max<mrate) missrate_max=mrate
}
/Imputation conflict rate/{
  crate=$4
  if (impconfl_min>crate) impconfl_min=crate
  if (impconfl_max<crate) impconfl_max=crate
}
END {
print "---"
print "- Summary:"
printf "    Chromosomes: "
printf chrom[1]
for (c=2; c<=nchrom; c++) printf " " chrom[c]
print ""
print "    Target samples min: " ntsamples_min
print "    Target samples max: " ntsamples_max
print "    Target variants: " ntvars
print "    Reference samples min: " nrsamples_min
print "    Reference samples max: " nrsamples_max
print "    Reference variants: " nrvars
print "    Common variants: " nphased
if (nimputed>0) print "    Exclusive reference variants: " nimputed
if (nimpout>0)  print "    Imputation output variants: " nimpout
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
print "    Missing rate in target genotypes min: " missrate_min
print "    Missing rate in target genotypes max: " missrate_max
print "    Imputation conflict rate min: " impconfl_min
print "    Imputation conflict rate max: " impconfl_max
}
