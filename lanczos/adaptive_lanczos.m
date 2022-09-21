## Copyright (C) 2022 david
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
## @deftypefn {} {@var{retval} =} adaptive_lanczos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: david <david@DAVE>
## Created: 2022-05-23

function [Am Bm Cm] = adaptive_lanczos (A,B,C,sigma,tol)
  m=1;
  eps=1;
  Hant=zeros(n,n);;
  [n n]=size(A);
  while(eps>tol)
    [V W]= symetric_lanczos(A,B,C,sigma,m);
    Am=W(:,:,m)'*A*V(:,:,m);
    Bm=W(:,:,m)'*B;
    Cm=C*V(:,:,m);
    H=Cm*inv((m*eye(n)-Am))*Bm;
    eps9norm(H-Hant);
    Hant=H;
    m=m+1;
  endwhile
endfunction
