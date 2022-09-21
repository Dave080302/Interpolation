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
## @deftypefn {} {@var{retval} =} unsymetric_lanczos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: david <david@DAVE>
## Created: 2022-05-23

function [Vseries, Wseries] = unsymetric_lanczos (A, V, W, m)
  [n n]=size(A);
  [n p] = size(V);
  Vseries=zeros(n,p,m);
  Wseries=zeros(n,p,m);
  Vaux=zeros(n,p,m+1); %seriile de matrice W~ si V~ folosite in pseudocod
  Waux=zeros(n,p,m+1);
  [alfa,beta] = qr(W'*V);
  Vseries(:,:,1)=V*inv(beta);
  Wseries(:,:,1)=W*alfa;
  Vaux(:,:,2)=A*V(:,:,1);
  Waux(:,:,2)=A'*W(:,:,1);
  for j=1:m-1
    alfaj=Wseries(:,:,j)'*Vaux(:,:,j+1);
    Vaux(:,:,j+1)=Vaux(:,:,j+1)-Vseries(:,:,j)*alfaj;
    Waux(:,:,j+1)=Waux(:,:,j+1)-Wseries(:,:,j)*alfaj;
    [Vseries(:,:,j+1) beta] = qr(Vaux(:,:,j+1));
    [Wseries(:,:,j+1) alfa] = qr(Waux(:,:,j+1));
    [U sigma Z] = svd(Wseries(:,:,j+1)'*Vseries(:,:,j+1));
    alfa=alfa*U*(sigma^(1/2));
    beta=(sigma^(1/2))*Z'*beta;
    Vseries(:,:,j+1)=Vseries(:,:,j+1)*Z*(sigma^(-1/2));
    Wseries(:,:,j+1)=Wseries(:,:,j+1)*U*(sigma^(-1/2));
    Vaux(:,:,j+2)=A*Vseries(:,:,j+1)-Vseries(:,:,j)*alfa;
    Waux(:,:,j+2)=A'*Wseries(:,:,j+1)-Wseries(:,:,j)*beta';
  endfor
  

endfunction
