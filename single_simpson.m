## Copyright (C) 2023 Antoine
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.


## Author: Antoine <Antoine@ANTONMACHINE>
## Created: 2023-04-02

function I = single_simpson(a, b, f)
  h = (b-a)/2;
  I = f(a);
  I = I + f(b);
  I = I + 4*f((a+b)/2);
  I = (h/3)*I;
endfunction
