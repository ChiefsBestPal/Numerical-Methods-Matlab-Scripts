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

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} composite_simpson2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Antoine <Antoine@ANTONMACHINE>
## Created: 2023-04-26

function I = composite_simpson (m,a,b,f)



  h=(b-a)/(2*m);
  x1=linspace(a+h,b-h,m);
  x2=linspace(a+2*h,b-2*h,m-1);
  I=h/3*(f(a)+f(b)+4*sum(f(x1))+2*sum(f(x2)));



    % h is the segment size
  # h = (b - a)/n;

  % X stores the summation of first and last segment
  #X = f(a)+f(b);

  % variables Odd and Even to store
  % summation of odd and even
  % terms respectively
  #Odd = 0;
  #Even = 0;
  #for i = 1:2:n-1
  #    xi=a+(i*h);
  #    Odd=Odd+f(xi);
  #endfor
  #for i = 2:2:n-2
  #    xi=a+(i*h);
   #   Even=Even+f(xi);
  #endfor

  % Formula to calculate numerical integration
  % using Simpsons 1/3 Rule
  #I = (h/3)*(X+4*Odd+2*Even);
endfunction
