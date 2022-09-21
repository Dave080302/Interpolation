function r = fxy(f, x, y)
    % =========================================================================
    % Aproximeaza derivata fata de x si y a lui f in punctul (x, y).
    % =========================================================================
    
    % TODO: Calculeaza derivata.
    [m n]=size(f);
    if (y<=1 || x<=1 || y>=m || x>=n)
      r=0; 
    else
      r=(f(y-1,x-1)+f(y+1,x+1)-f(y-1,x+1)-f(y+1,x-1))/4;
    endif
    
endfunction