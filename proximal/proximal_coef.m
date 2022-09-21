function a = proximal_coef(f, x1, y1, x2, y2)
    % =========================================================================
    % Calculeaza coeficientii a pentru Interpolarea Proximala intre punctele
    % (x1, y1), (x1, y2), (x2, y1) si (x2, y2).
    % =========================================================================
    
    % TODO: Calculeaza matricea A.
    A=[1 x1 y1 x1*y1; 1 x1 y2 x1*y2; 1 x2 y1 x2*y1; 1 x2 y2 x2*y2];
    % TODO: Calculeaza vectorul b.  
    b=zeros(4,1);  
    b(1)=f(y1,x1);
    b(2)=f(y2,x1);
    b(3)=f(y1,x2);
    b(4)=f(y2,x2);
    % TODO: Calculeaza coeficientii.
    a=inv(A)*b;
endfunction
