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
  exmarefonly=0
  dropnotfound=0
  dropmatgt=0
  dropmaref=0
  dropmono=0
  dropmismatch=0
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
  FS="[<>]+"; $0=$0
  if (ntsamples_min>$7) ntsamples_min=$7
  if (ntsamples_max<$7) ntsamples_max=$7
}
/Target variants/{
  FS="[<>]+"; $0=$0
  ntvars+=$7
}
/Reference samples/{
  FS="[<>]+"; $0=$0
  if (nrsamples_min>$7) nrsamples_min=$7
  if (nrsamples_max<$7) nrsamples_max=$7
}
/Reference variants/{
  FS="[<>]+"; $0=$0
  nrvars+=$7
}
/were used for phasing/{
  FS="[> ]"; $0=$0
  nphased+=$4
}
/were imputed/{
  FS="[> ]"; $0=$0
  nimputed+=$4
}
/Imputation output contains/{
  FS="[ ]+"; $0=$0
  nimpout+=$4
}
/REF\/ALT were swapped AND/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    refaltsflip+=$10
  }
}
/REF\/ALT were swapped in/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    refalt+=$6
  }
}
/Strands were flipped in/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    sflip+=$6
  }
}
/variants not found in reference/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    dropnotfound+=$3
  }
}
/multi-allelic variants in target/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    dropmatgt+=$3
  }
}
/but multi-allelic in reference/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    dropmaref+=$3
  }
}
/monomorphic in reference/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    dropmono+=$3
  }
}
/variants with allele mismatches/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    dropmismatch+=$3
  }
}
/reference-only multi-allelic variants/{
  if (summaryline>0) {
    FS="[ ]+"; $0=$0
    exmarefonly+=$3
  }
}
/missing rate in target genotypes/{
  FS="[ <]+"; $0=$0
  mrate=$(NF-1)
  if (missrate_min>mrate) missrate_min=mrate
  if (missrate_max<mrate) missrate_max=mrate
}
/No imputation conflicts/{
  if (impconfl_min>0.0) impconfl_min=0.0
  # max will be unchanged here
}
/Imputation conflict rate/{
  FS="[ <]+"; $0=$0
  crate=$(NF-1)
  if (impconfl_min>crate) impconfl_min=crate
  if (impconfl_max<crate) impconfl_max=crate
}
END {
print "<h3>SUMMARY:</h3>"
print "<table>"
printf "<tr><td>Chromosomes:</td> <td>"
printf chrom[1]
for (c=2; c<=nchrom; c++) printf ", " chrom[c]
print "</td></tr>"
printf "<tr><td>Target samples:</td> <td>" ntsamples_min
if (ntsamples_min != ntsamples_max) printf " - " ntsamples_max
print "</td></tr>"
print "<tr><td>Target variants:</td> <td>" ntvars "</td></tr>"
printf "<tr><td>Reference samples:</td> <td>" nrsamples_min
if (nrsamples_min != nrsamples_max) printf " - " nrsamples_max
print "</td></tr>"
print "<tr><td>Reference variants:</td> <td>" nrvars "</td></tr>"
print "</table>"
print "<p class='pinfo'>Variants common in target and reference used for phasing: <b>" nphased "</b></p>"
if (nimputed>0) print "<p class='pinfo'>Variants exclusively in reference that were imputed: <b>" nimputed "</b></p>"
if (nimpout>0) print "<p class='pinfo'>Variants in imputation output: <b>" nimpout "</b></p>"
print "<p class='pinfo'>"
if (refaltsflip>0) print "REF/ALT were swapped AND strands were flipped in " refaltsflip " variants.<br>"
if (refalt>0) print "REF/ALT were swapped in " refalt " variants.<br>"
if (sflip>0) print "Srands were flipped in " sflip " variants.<br>"
print "Dropped " dropnotfound " variants not found in reference.<br>"
if (dropmatgt>0) print "Dropped " dropmatgt " multi-allelic variants in target.<br>"
if (dropmaref>0) print "Dropped " dropmaref " variants bi-allelic in target but multi-allelic in reference.<br>"
if (dropmono>0) print "Dropped " dropmono " variants bi-allelic in target but monomorphic in reference.<br>"
if (dropmismatch>0) print "Dropped " dropmismatch " variants with allele mismatches.<br>"
print "Excluded " exmarefonly " reference-only multi-allelic variants from imputation.<br>"
print "</p><p class='pinfo'>"
printf "Av. missing rate in target genotypes: " missrate_min
if (missrate_max > missrate_min) printf " - " missrate_max
print "<br>"
printf "Imputation conflict rate: " impconfl_min
if (impconf_max > impconfl_min) printf " - " impconfl_max
print "</p>"
}
