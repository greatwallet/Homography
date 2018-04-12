function [H] = computeHomography(x1,y1,x2,y2)
% input 4 points 

    A = zeros(8, 8);
    for i = 1:4
        A((2 * i) - 1,:) = [x1(i), y1(i), 1, 0, 0, 0, -x1(i) * x2(i), -x2(i) * y1(i)]; 
        A(2 * i,:) = [0, 0, 0, x1(i), y1(i), 1, -x1(i) * -y2(i), -y2(i) * y1(i)];

    end
    
    B = [x2(1); y2(1); x2(2); y2(2); x2(3); y2(3); x2(4); y2(4);];
    h = A\B;
    h = [h; 1];
    %H = reshape(h, [3, 3]);
    %{
    H = [h(1), h(4), h(7);
        h(2), h(5), h(8);
        h(3), h(6), 1;];
    %}
    H = [h(1), h(2), h(3);
        h(4), h(5), h(6);
        h(7), h(8), 1;];
    
end