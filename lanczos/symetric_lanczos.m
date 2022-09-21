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
## @deftypefn {} {@var{retval} =} symetric_lanczos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: david <david@DAVE>
## Created: 2022-05-23

function [V,W] = symetric_lanczos (A,B,C,sigma,m)
  [n n] = size(A);
  [n p] = size(B);
  V= zeros(n,p,m+1);
  W= zeros(n,p,m+1);
  S=inv((A-sigma(:,:,1)*eye(n)))*B;
  R=(inv((A-sigma(:,:,1))))'*C';
  V(:,:,1)=S;
  W(:,:,1)=R;
  for k=1:m
     if(k<m)
      if(sigma(:,:,k+1)==inf)
      S=A*V(:,:,k);
      R=A'*W(:,:,k);
      else
      S=inv((A-sigma(:,:,k+1)*eye(n)))*V(:,:,k);
      R=(inv((A-sigma(:,:,k+1)*eye(n))))'*W(:,:,k);
      endif
      H=W(:,:,k)'*S;
      G=V(:,:,k)'*R;
      S=S-V(:,:,k)*H;
      R=R-W(:,:,k)*G;
      [V(:,:,k+1) H] = qr(S);
      [P D Q]= svd(W(:,:,k+1)'*V(:,:,k+1));
      V(:,:,k+1)=V(:,:,k+1)*Q*(D^((-1)*(1/2)));
      W(:,:,k+1)=W(:,:,k+1)*P*(D^((-1)*(1/2)));
      H=D^(1/2)*Q'*H;
      G=D^(1/2)*P'*G;
     else
       if(sigma(:,:,k+1)==inf)
        S=A*B;
        R=A'*C;
       else
        S=inv(A)*B;
        R=inv(A)'*C';
       endif
       H=W(:,:,m)'*S;
       G=V(:,:,m)'*R;
     endif
     H=W(:,:,m)'*S;
      G=V(:,:,m)'*R;
      S=S-V(:,:,m)*H;
      R=R-W(:,:,m)*G;
      [V(:,:,m+1) H] = qr(S);
      [P D Q]= svd(W(:,:,m+1)'*V(:,:,m+1));
      V(:,:,m+1)=V(:,:,m+1)*Q*(D^((-1)*(1/2)));
      W(:,:,m+1)=W(:,:,m+1)*P*(D^((-1)*(1/2)));
  endfor
endfunction
