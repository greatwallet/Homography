function [x2,y2] = applyHomography(H,x1,y1)
    x2 = (H(1,1) * x1 + H(1,2) * y1 + H(1,3)) ./ (H(3,1) * x1 + H(3,2) * y1 + H(3,3));
    
    y2 = (H(2,1) * x1 + H(2,2) * y1 + H(2,3)) ./ (H(3,1) * x1 + H(3,2) * y1 + H(3,3));
    
    %[r,c] = size

end