## Copyright (C) 2023 antoi
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

## Author: antoi <antoi@DESKTOP-4HB4J2I> Antoine Cantin
## Created: 2023-03-18
function [A,a,residuals,rmse] = NormalInconsistentEqs (varargin)
    % rmse is root mean squared error computed with norm and residuals
    % a (a_hat) is model parameters
    % residuals is absolute err vector wise (litteraly 'residuals')
    % x is vector
    % y is vector
    % A is inconsistent coefficient matrix
    [x,y] = varargin{1:2};
    #x = standardize_to_rowvector(x);
    #y = standardize_to_rowvector(y);
    n = length(x);
    assert (n == length(y));
    A = [varargin{3}(x)];
    for i = 4:length(varargin)
      A = horzcat(A,varargin{i}(x));
    endfor
    a = A'*A\A'*y; #PA=LU solve of each side transposed: Normal Eqs
    residuals = (A*a)-y;
    rmse = norm(residuals,2) / sqrt(n);


endfunction
