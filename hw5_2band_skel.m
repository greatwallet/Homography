
%
% Some sketch of how to do two-band blending for a pair of warped images.  
%
% Thise example is just for two images.  Once you get this working
% you will need to incorporate this into the code you wrote in part 1 & 2
% for blending together multiple images where you loops over all the images
%

% Let I1 and I2 be the warped images
% Let M1 and M2 be the correpsonding masks
% J = warped images
I1 = J{1};
I2 = J{2};
I3 = J{3};
M1 = rgb2gray(J{1});
M2 = rgb2gray(J{2});
M3 = rgb2gray(J{3});

%
% separate out frequency bands for the image and mask
%
% h = ones(5,5) / 25;
gFilt = fspecial('gaussian', [5,5], 1);

% low frequency band is just blurred version of the image
I1_L = imfilter(I1, gFilt);   
I2_L = imfilter(I2, gFilt); 
I3_L = imfilter(I3, gFilt);

% high frequency band is whatever is left after subtracting the low frequencies
I1_H = I1-I1_L;  
I2_H = I2-I2_L;
I3_H = I3-I3_L;

% low frequency alpha should be feathered version of M1 & M2
A1_L = imfilter(M1,gFilt);  
A2_L = imfilter(M2, gFilt); 
A3_L = imfilter(M3, gFilt);

% normalize the alpha masks to sum to 1 at every pixel location
% (we avoid dividing by zero by adding a term to the denominator 
% anyplace the sum is 0)
Asum = A1_L + A2_L + A3_L;  
A1_L = A1_L ./ (Asum + (Asum==0));
A2_L = A2_L ./ (Asum + (Asum==0));
A3_L = A3_L ./ (Asum + (Asum==0));



%for high frequencies use a very sharp alpha mask which is 
% alpha=1 for which ever image has the most weight at each 
% location
A1_H = double(A1_L > A2_L); 
A2_H = double(A2_L > A1_L);
A3_H = double(A3_L > A1_L);

% normalize the alpha masks to sum to 1 
% technically we shouldn't have to do this the way we've constructed
% A1_H and A2_H above, but just to be safe.
Asum = A1_H + A2_H + A3_H;  
A1_H = A1_H ./ (Asum + (Asum==0));
A2_H = A2_H ./ (Asum + (Asum==0));
A3_H = A3_H ./ (Asum + (Asum==0));

%
% now combine the results using alpha blending
%

J_L = A1_L .* I1_L  + A2_L .* I2_L + A3_L .* I3_L;  % low frequency band blend result
J_H = A1_H .* I1_H  + A2_H .* I2_H + A3_H .* I3_H;  % high frequency band blend result
J = J_L + J_H; % combined bands


%
% display some of the intermediate results 
%
figure(1);
subplot(3,3,1); imshow(I1_L); 
subplot(3,3,2); imshow(I2_L);
subplot(3,3,3); imshow(J_L);  title('low frequency band');

subplot(3,3,4); imshow(I1_H);
subplot(3,3,5); imshow(I2_H);
subplot(3,3,6); imshow(J_H); title('high frequency band');

subplot(3,3,4); imshow(I1);
subplot(3,3,5); imshow(I2);
subplot(3,3,6); imshow(J); title('combined');

figure(2);
imshow(J);