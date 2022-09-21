function out = proximal_2x2_RGB(img, STEP = 0.1)
    % ==============================================================================================
    % Aplica Interpolare Proximala pe imaginea 2 x 2 definita img cu puncte intermediare echidistante.
    % img este o imagine colorata RGB -Red, Green, Blue.
    % =============================================================================================

    % TODO: Extrage canalul rosu al imaginii.
    red=img(:,:,1);
    % TODO: Extrage canalul verde al imaginii.
    green=img(:,:,2);
    % TODO: Extrage canalul albastru al imaginii.
    blue=img(:,:,3);
    % TODO: Aplic? functia proximal pe cele 3 canale ale imaginii.
    R=proximal_2x2(red,STEP);
    G=proximal_2x2(green,STEP);
    B=proximal_2x2(blue,STEP);
    % TODO: Formeaza imaginea finala concatenând cele 3 canale de culori.
    out=cat(3,R,G,B);
endfunction
