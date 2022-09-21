function out = proximal_rotate_RGB(img, rotation_angle)
    % =========================================================================
    % Aplica Interpolarea Proximala pentru a roti o imagine RGB cu un unghi dat.
    % =========================================================================
    
    % TODO: Extrage canalul rosu al imaginii.
    red=img(:,:,1);
    % TODO: Extrage canalul verde al imaginii.
    green=img(:,:,2);
    % TODO: Extrage canalul albastru al imaginii.
    blue=img(:,:,3);
    % TODO: Aplica rotatia pe fiecare canal al imaginii.
    R=proximal_rotate(red,rotation_angle);
    G=proximal_rotate(green,rotation_angle);
    B=proximal_rotate(blue,rotation_angle);
    % TODO: Formeaza imaginea finala concatenând cele 3 canale de culori.
    out=cat(3,R,G,B);
endfunction