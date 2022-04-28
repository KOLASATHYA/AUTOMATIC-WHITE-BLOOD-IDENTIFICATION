
clc;
clear all;
close all;
[fna,pna]=uigetfile({'*.jpg;*.bmp;*.png'},'Select file');
fanme=[pna fna];
im = double(imread(fanme));
im=imresize(im,[200 300]);
figure;imshow(uint8(im));title('Input image');
tic
%  Lab conversion and color thresholding
im=im/255;
R=im(:,:,1);
G=im(:,:,2);
B=im(:,:,3);
% Set a threshold
T = 0.008856;

[M, N,plane] = size(R);
s = M * N;
RGB = [reshape(R,1,s); reshape(G,1,s); reshape(B,1,s)];

% RGB to XYZ
MAT = [0.412453 0.357580 0.180423;
       0.212671 0.715160 0.072169;
       0.019334 0.119193 0.950227];
XYZ = MAT * RGB;

% Normalize for D65 white point
X = XYZ(1,:) / 0.950456;
Y = XYZ(2,:);
Z = XYZ(3,:) / 1.088754; 

XT = X > T;
YT = Y > T;
ZT = Z > T;

Y3 = Y.^(1/3); 

fX = XT .* X.^(1/3) + (~XT) .* (7.787 .* X + 16/116);
fY = YT .* Y3 + (~YT) .* (7.787 .* Y + 16/116);
fZ = ZT .* Z.^(1/3) + (~ZT) .* (7.787 .* Z + 16/116);

L = reshape(YT .* (116 * Y3 - 16.0) + (~YT) .* (903.3 * Y), M, N);
a = reshape(500 * (fX - fY), M, N);
b = reshape(200 * (fY - fZ), M, N);

  L = cat(3,L,a,b);
figure;imshow(L);title('L*a*b converted image'); 

% color thresholding

      out2 = colorthreshold(lab2uint8(L),'otsu');
      figure, imshow(uint8(out2));
      title('Color Thresholded image');
     OUCT=uint8(out2); 
     ReLAb=zeros(size(L));
     ReRGB=zeros(size(L));
Th=100;
for i=1:M
    for j=1:N
        if OUCT(i,j,1)<200 && OUCT(i,j,1)>90
            ReLAb(i,j,:)=L(i,j,:);
            ReRGB(i,j,:)=im(i,j,:);
        end
    end
end

   figure, imshow(ReLAb);title('L*a*b Mask image');
   figure, imshow(ReRGB);title('RGB Mask image');
   
   

%----RGB 2 Gray conversion

Grayim=rgb2gray(ReRGB);
figure, imshow(Grayim);title('Gray Mask image');
  

L=imadjust(Grayim);
H = histeq(Grayim);
%% Brighten most of the details in the image except the nucleus
R1=imadd(L,H);
%% Highlight all the objects and its borders in the image including the cell nucleus
R2 = imsubtract(R1,H);
%% Remove almost all the other blood components while retaining the nucleus with minimum affect of distortion on the nucleus part of the white blood cell.
% R3=imadd(R1,R2);
figure;
subplot(221);
imshow(Grayim);title('Gray image')
subplot(222);
imshow(H);title('Enhanced image')
subplot(223);
imshow(R1);title('Image except the nucleus')
subplot(224);
imshow(R2);title('Highlighted all objects')
toc
%check histogram
figure;
hi=imhist(R2);
hi1=hi(1:2:256);
horz=1:2:256;
bar(horz,hi1);
axis([0 255 0 1400]);
set(gca, 'xtick', 0:50:255);
set(gca, 'ytick', 0:2000:15000);
 xlabel('Gray level' );
ylabel('No of pixels' );
title('Histogram before opening the image');

%=====================

% Otsu thresholding   
   rk=otsu(R2);
  imthO=(Grayim>rk); 
figure, imshow(imthO);title('Otsu Binary Image');
  
%implements a 3-by-3 average filter
R3 = average_filter(imthO);
R31=double(im2bw(R3));
figure; imshow(R31);title('Average Filtered image');

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(R3), hy, 'replicate');
Ix = imfilter(double(R3), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure
imshow(gradmag,[]), title('Gradient magnitude')
% Finally we are ready to compute the watershed-based segmentation.
L = watershed(gradmag);
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Watershed segmentation')
%  Mark the Foreground Objects
se = strel('disk', 2);
Io = imopen(Grayim, se);
figure
imshow(Io), title('Opening (Io)')

Ie = imerode(Grayim, se);
Iobr = imreconstruct(Ie, Grayim);
figure
imshow(Iobr), title('Opening-by-reconstruction (Iobr)')

Ioc = imclose(Io, se);
figure
imshow(Ioc), title('Opening-closing (Ioc)')

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure
imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

fgm = Iobrcbr>0.1;
figure
imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')

I2 = Grayim;
I2(fgm) = 255;
figure
imshow(I2), title('Regional maxima superimposed on original image (I2)')


% se2 = strel('disk', 1);
% % fgm2 = imclose(fgm, se2);
% fgm = imerode(fgm, se2);
% This procedure tends to leave some stray isolated pixels that must be removed. You can do this using bwareaopen, which removes all blobs that have fewer than a certain number of pixels.
% getin=inputdlg('Enter the Area opening dimension(small objects-100 & Large objects-250)');
% aREAIN=str2double(getin{1,1});
aREAIN=50;
fgm4 = bwareaopen(fgm,aREAIN);
I3 = Grayim;
I3(fgm4) = 255;
figure
imshow(I3)
title('Modified regional maxima superimposed on original image (fgm4)')
figure
imshow(fgm4)
title('After Eliminating Platelets')
% 
% Compute Background Markers
% Now you need to mark the background. In the cleaned-up image, Iobrcbr, the dark pixels belong to the background, so you could start with a thresholding operation.


D = bwdist(fgm4);
DL = watershed(D);
bgm = DL == 0;
figure
imshow(bgm), title('Broken Watershed Nucleus')
%% eliminate border touching segments
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(fgm4), hy, 'replicate');
Ix = imfilter(double(fgm4), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure
imshow(gradmag,[]), title('Gradient magnitude')
% Finally we are ready to compute the watershed-based segmentation.
L = watershed(gradmag);
% Step 6: Visualize the Result
% One visualization technique is to superimpose the foreground markers, background markers, and segmented object boundaries on the original image. You can use dilation as needed to make certain aspects, such as the object boundaries, more visible. Object boundaries are located where L == 0.
I4 = Grayim;
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
figure
imshow(I4)
title('Markers and object boundaries superimposed on original image ')

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')

%% Final 
figure;
imshow(im);
hold on
boundaries = bwboundaries(fgm4);	
numberOfBoundaries = size(boundaries);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 1.5);
end
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=imthO;
aa=imresize(a,[256 256]);
u_bw_filename = aa;

b=im;
 bb=imresize(b,[256 256]);
u_GT_filename = im2bw(bb);

u_GT = [((u_GT_filename)) > 0 ];
u_bw = [((u_bw_filename)) > 0 ];

temp_obj_eval = objective_evaluation_core(u_bw, u_GT);

disp('PSNR value --');
disp(temp_obj_eval.PSNR);
