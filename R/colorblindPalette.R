##==============================================================================
## Source:
## http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
## http://www.somersault1824.com/wp-content/uploads/2015/02/color-blindness-palette.png
##
## Code by Tony Wong (twong@psu.edu)
##==============================================================================
#===============================================================================
# Copyright 2017 Tony Wong
#
# MESS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MESS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MESS.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================
mycol=rbind(
              c(0,0,0),
              c(0,73,73),
              c(0,146,146),
              c(255,109,182),
              c(255,182,119),
              c(73,0,146),
              c(0,109,219),
              c(182,109,255),
              c(109,182,255),
              c(182,219,255),
              c(146,0,0),
              c(146,73,0),
              c(219,109,0),
              c(36,255,36),
              c(255,255,109)
            )
mycol=mycol/max(mycol)

mycol.rgb <- rep(0,nrow(mycol))
for (i in 1:nrow(mycol)) {
    mycol.rgb[i] <- rgb(mycol[i,1],mycol[i,2],mycol[i,3])
}

##==============================================================================
## End
##==============================================================================
