%Loading images
imnames = {'IMG_1347.JPG','IMG_1348.JPG','IMG_1349.JPG'};
nimages = length(imnames);
baseim = 1; %index of the central "base" imag


for i = 1:nimages
  ims{i} = im2double(imread(imnames{i}));
  ims_gray{i} = rgb2gray(ims{i});
  [h(i),w(i),~] = size(ims{i});
end

% 
% %movingPoints = cell(nimages);
% %fixedPoints = cell(nimages);
% % get corresponding points between each image and the central base image
% for i = 1:nimages
%    if (i ~= baseim)
%      %run interactive select tool to click corresponding points on base and non-base image
%      .... = cpselect(ims{baseim},ims{i},'Wait',true);
%      index = i - 1;
%      [movingPoints, fixedPoints] = cpselect(ims{baseim},ims{i},'Wait',true);
%      data{index} = [movingPoints, fixedPoints];
%      %print(movingPoints);
%      %refine the user clicks using cpcorr
%      movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,ims_gray{i},ims_gray{baseim});
%      %print(movingPointsAdjusted);
%      %dataPoints(i) = movingPointsAdjusted;
%    end
% end
% 
% %
% % verify that the points are good!
% % some example code to plot some points, you will need to modify
% % this based on how you are storing the points etc.
% %
% 
% figure;
% subplot(2,1,1); 
% imagesc(ims{1});
% hold on;
% plot(data{1,1}(1,1), data{1,1}(1,2), 'b*', data{1,1}(2,1), data{1,1}(2,2),...
%     'r*', data{1,1}(3,1), data{1,1}(3,2), 'g*', data{1,1}(4,1), data{1,1}(4,2), 'y*')
% subplot(2,1,2);
% 
% imshow(ims{2});
% hold on;
% plot(data{1,1}(1,3), data{1,1}(1,4), 'b*', data{1,1}(2,3), data{1,1}(2,4),...
%     'r*', data{1,1}(3,3), data{1,1}(3,4), 'g*', data{1,1}(4,3), data{1,1}(4,4), 'y*') 
% 
% figure;
% subplot(2,1,1); 
% imagesc(ims{1});
% hold on;
% plot(data{1,2}(1,1), data{1,2}(1,2), 'b*', data{1,2}(2,1), data{1,2}(2,2),...
%     'r*', data{1,2}(3,1), data{1,2}(3,2), 'g*', data{1,2}(4,1), data{1,2}(4,2), 'y*')
% subplot(2,1,2);
% 
% imshow(ims{3});
% hold on;
% plot(data{1,2}(1,3), data{1,2}(1,4), 'b*', data{1,2}(2,3), data{1,2}(2,4),...
%     'r*', data{1,2}(3,3), data{1,2}(3,4), 'g*', data{1,2}(4,3), data{1,2}(4,4), 'y*')
% 
% 
% % at this point it is probably a good idea to save the results of all your clicking
% % out to a file so you can easily load them in again later on
% %
% 
% %save mypts.mat x1 y1 x2 y2;
% save mypts.mat data
% 
% to reload the points: 

load mypts.mat;
%

for i = 1:nimages
    if (i ~= baseim)
        index = i - 1;
        H{i} = computeHomography(data{1,index}(:,3), data{1,index}(:,4),data{1,index}(:,1),data{1,index}(:,2));
        %H{i} = computeHomography(x1,x2,y1,y2);
    else
        H{i} = eye(3); %homography for base image is just the identity matrix
    end
end

%
% compute where corners of each warped image end up
%

for i = 1:nimages
  cx = [1 1 w(i) w(i)];  %corner coordinates based on h,w for each image
  cy = [1 h(i) 1 h(i)];
  
  [cx_warped{i},cy_warped{i}] = applyHomography(H{i},cx,cy);
end

% 
% find corners of a rectangle that contains all the warped image
%  corner points
%

totalXmin = [];
totalYmin = [];
totalXmax = [];
totalYmax = [];
for i = 1:nimages
    totalXmin = [totalXmin min(cx_warped{i})];
    totalYmin = [totalYmin min(cy_warped{i})];
    totalXmax = [totalXmax max(cx_warped{i})];
    totalYmax = [totalYmax max(cy_warped{i})];
    
    totalXmin = round(min(totalXmin));
    totalYmin = round(min(totalYmin));
    totalXmax = round(max(totalXmax));
    totalYmax = round(max(totalYmax));
end

% Use H and interp2 to perform inverse-warping of the source image to align it with the base image


[xx,yy] = meshgrid(totalXmin:totalXmax, totalYmin:totalYmax);  %range of meshgrid should be the containing rectangle
for i = 1:nimages
    [wx,wy] = applyHomography(inv(H{i}), xx, yy);
    R = interp2(ims{i}(:,:,1),wx,wy);
    G = interp2(ims{i}(:,:,2),wx,wy);
    B = interp2(ims{i}(:,:,3),wx,wy);
    J{i} = cat(3,R,G,B);
    
    mask{i} = ~isnan(R);  %interp2 puts NaNs outside the support of the warped image
    J{i}(isnan(J{i})) = 0;
end
 
%imagesc(J{1})
   
for i = 1:nimages
     %blur and clip mask{i} to get an alpha map for each image
    gFilt = fspecial('gaussian');
    alpha{i} = imfilter(mask{i}, gFilt);
    locations = find(~mask{i});
    alpha{i}(locations) = 0;
end

% scale alpha maps to sum to 1 at every pixel location
A = alpha{1};
for i = 2:nimages
    A = cat(3, A, alpha{i});
end

% normalize at ech location to not divide by zero
B = sum(A, 3);
B(B == 0) = 1;
B = repmat(B, [1, 1, nimages]);
A = A./B;

for i = 1 : nimages
        alpha{i} = A(:,:,i);
end


% finally blend together the resulting images into the final mosaic


K = zeros(size(J{1}));
for i = 1:nimages
    K = K + repmat(alpha{i}, [1, 1, 3]) .* J{i};
end

% display the result
figure(1); 
imagesc(K); axis image;

% save the result to include in your writeup
imwrite(K, 'final.jpg')
%%%%%%%%%%
% I1 = J{1};
% I2 = J{2};
% I3 = J{3};
% M1 = rgb2gray(J{1});
% M2 = rgb2gray(J{2});
% M3 = rgb2gray(J{3});
% 
% %
% % separate out frequency bands for the image and mask
% %
% % h = ones(5,5) / 25;
% gFilt = fspecial('gaussian', [5,5], 1);
% 
% % low frequency band is just blurred version of the image
% I1_L = imfilter(I1, gFilt);   
% I2_L = imfilter(I2, gFilt); 
% I3_L = imfilter(I3, gFilt);
% 
% % high frequency band is whatever is left after subtracting the low frequencies
% I1_H = I1-I1_L;  
% I2_H = I2-I2_L;
% I3_H = I3-I3_L;
% 
% % low frequency alpha should be feathered version of M1 & M2
% A1_L = imfilter(M1,gFilt, 'same');  
% A2_L = imfilter(M2, gFilt, 'same'); 
% A3_L = imfilter(M3, gFilt, 'same');
% 
% % normalize the alpha masks to sum to 1 at every pixel location
% % (we avoid dividing by zero by adding a term to the denominator 
% % anyplace the sum is 0)
% Asum = A1_L + A2_L + A3_L;  
% A1_L = A1_L ./ (Asum + (Asum==0));
% A2_L = A2_L ./ (Asum + (Asum==0));
% A3_L = A3_L ./ (Asum + (Asum==0));
% 
% 
% 
% %for high frequencies use a very sharp alpha mask which is 
% % alpha=1 for which ever image has the most weight at each 
% % location
% A1_H = double(A1_L > A2_L); 
% A2_H = double(A2_L > A1_L);
% A3_H = double(A3_L > A1_L);
% 
% % normalize the alpha masks to sum to 1 
% % technically we shouldn't have to do this the way we've constructed
% % A1_H and A2_H above, but just to be safe.
% Asum = A1_H + A2_H + A3_H;  
% A1_H = A1_H ./ (Asum + (Asum==0));
% A2_H = A2_H ./ (Asum + (Asum==0));
% A3_H = A3_H ./ (Asum + (Asum==0));
% 
% %
% % now combine the results using alpha blending
% %
% 
% J_L = A1_L .* I1_L  + A2_L .* I2_L + A3_L .* I3_L;  % low frequency band blend result
% J_H = A1_H .* I1_H  + A2_H .* I2_H + A3_H .* I3_H;  % high frequency band blend result
% J = J_L + J_H; % combined bands
% 
% 
% %
% % display some of the intermediate results 
% %
% figure(1);
% subplot(3,3,1); imshow(I1_L); 
% subplot(3,3,2); imshow(I2_L);
% subplot(3,3,3); imshow(J_L);  title('low frequency band');
% 
% subplot(3,3,4); imshow(I1_H);
% subplot(3,3,5); imshow(I2_H);
% subplot(3,3,6); imshow(J_H); title('high frequency band');
% 
% subplot(3,3,4); imshow(I1);
% subplot(3,3,5); imshow(I2);
% subplot(3,3,6); imshow(J); title('combined');
% 
% figure(2);
% imshow(J);
