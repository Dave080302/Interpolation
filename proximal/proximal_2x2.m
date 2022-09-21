function out = proximal_2x2(f, STEP = 0.1)
    % ===================================================================================
    % Aplica Interpolare Proximala pe imaginea 2 x 2 f cu puncte intermediare echidistante.
    % f are valori cunoscute în punctele (1, 1), (1, 2), (2, 1) ?i (2, 2).
    % Parametrii:
    % - f = imaginea ce se va interpola;
    % - STEP = distan?a dintre dou? puncte succesive.
    % ===================================================================================
    
    % TODO: Defineste coordonatele x si y ale punctelor intermediare.
    x_int=1: STEP : 2;
    y_int=1: STEP : 2;
    % Se afl? num?rul de puncte.
    n = length(x_int);
    % TODO: Cele 4 puncte încadratoare vor fi aceleasi pentru toate punctele din interior.
    x1=1;
    y1=1;
    x2=1;
    y2=2;
    x3=2;
    y3=1;
    x4=2;
    y4=2;
    % TODO: Initializeaza rezultatul cu o matrice nula n x n.
    out=zeros(n,n);
    % Se parcurge fiecare pixel din imaginea finala.
    for i = 1 : n
        for j = 1 : n
        % TODO: Afla cel mai apropiat pixel din imaginea initiala.
        % TODO: Calculeaza pixelul din imaginea finala.
        d(1)=sqrt((x_int(i)-x1)^2+(y_int(j)-y1)^2) ;
        d(2)=sqrt((x_int(i)-x2)^2+(y_int(j)-y2)^2) ;
        d(3)=sqrt((x_int(i)-x3)^2+(y_int(j)-y3)^2) ;
        d(4)=sqrt((x_int(i)-x4)^2+(y_int(j)-y4)^2) ;
        min_dist=min(d);
        if(j!=ceil(n/2) && i!=ceil(n/2))
          if(min_dist==d(1))
            out(i,j)=f(x1,y1);
          elseif(min_dist==d(2))
            out(i,j)=f(x2,y2);
          elseif(min_dist==d(3))
            out(i,j)=f(x3,y3);
          else
            out(i,j)=f(x4,y4);
          endif   
         elseif(i==ceil(n/2) && j==ceil(n/2))
          out(i,j)=f(x4,y4);
         elseif(i==ceil(n/2))
             if(j<=ceil(n/2))
                out(i,j)=f(x3,y3);
             elseif(j>(n/2))
                out(i,j)=f(x4,y4);
             endif   
        elseif (j==ceil(n/2)) 
              if (i<=ceil(n/2))
                out(i,j)=f(x2,y2);
              elseif (i>(n/2))
                out(i,j)=f(x4,y4);
              endif
      
      endif   
     endfor
    endfor

endfunction