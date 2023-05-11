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
print "---" # YAML header
}

{
  # in general, print every line as it is, but filter the YAML header of each input file
  # and apply the chromosome tag fetched from the file prefix to the context field
  if (FNR==1) { # set chromosome name at start of each file
    # get the basename of the file and split into fields separated by '.'
    split(FILENAME,fnfields,"/")
    split(fnfields[length(fnfields)],fnfieldsx,".")
    chr=fnfieldsx[1]
  }

  if ($1 == "Context:") {
    # modify context field
    $1=""
    print "    Context: " chr $0
  } else if ($1 != "---")
    print
}
