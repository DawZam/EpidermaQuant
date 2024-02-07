%% EpidermaQuant
%  Author: Dawid Zamojski
%  Email: dawid.zamojski@polsl.pl
%  ======================

clc
clear variables
close all

%% 0. IMAGE SELECTION 

% image selection for analysis
[I, path] = uigetfile('*.jpg', 'Choose one image');
str = strcat(path, I);
image = imread(str);

% a -> FLG; b -> K10; c -> Ki67 (the best marker for reference)
marker = 'c';

% select an appropriate group for analysis:
gr_badawcza_wybor = {'FLG', 'K10', 'Ki67', 'HSPA2', 'nowe_FLG'};
gr_badawcza = 'FLG';

%% 1A. INITIAL COLOR NORMALIZATION

% color normalization
normal_Reinhard = my_norm(image, marker, 'A');
normal_RGBHist = my_norm(image, marker, 'B');
normal_Macenko = my_norm(image, marker, 'C');
%normal_SCD = Norm(image, target, 'SCD');

% results
figure()
subplot(1,3,1); imshow(normal_Reinhard); title('Reinhard color normalization') 
subplot(1,3,2); imshow(normal_RGBHist); title('RGBHist color normalization')
subplot(1,3,3); imshow(normal_Macenko); title('Macenko color normalization')

% image choices based on results from color normalization
image = uint8(255 * mat2gray(normal_Reinhard)); % definitely the best
%image = normal_RGBHist; % also exposes negative samples
%image = normal_Macenko; % dedicated to H&E images

%% 1B. PLASTIC-MEMBRANE CUTTING OFF (drawassisted oraz roipoly) 

% ranges to remove plastic
if gr_badawcza(1:3) == gr_badawcza_wybor{1} % FLG
    left = 206;
    right = 225;
    tresh = 0.6;
elseif gr_badawcza(1:3) == gr_badawcza_wybor{2} % K10
    left = 203;
    right = 224.5;
    tresh = 0.55;
elseif gr_badawcza(1:4) == gr_badawcza_wybor{3} % Ki67
    left = 212;
    right = 225;
    tresh = 0.5;
elseif gr_badawcza(1:5) == gr_badawcza_wybor{4} % HSPA2
    left = 222;
    right = 228;
    tresh = 0.55;
elseif gr_badawcza(1:8) == gr_badawcza_wybor{5} % nowe_FLG
    left = 219;
    right = 223;
    tresh = 0.6;
end

% function to find ares of plastic
plastic_det = my_plastic_range(image, left, right, tresh);

% results
figure(); imshow(plastic_det); title("Effects of morphological operations on R-channel - plastic")

%% 1C. DECONVOLUTION 
% image division into channels for deconvolution
image2 = double(image);
ImgR = image2(:,:,1);
ImgG = image2(:,:,2);
ImgB = image2(:,:,3);

% deconvolution function without color normalization
[ImgR_back, ImgG_back, ImgB_back, Dye01_transmittance, Dye02_transmittance, Dye03_transmittance, LUTdye01, LUTdye02, LUTdye03, Q3x3Mat] = Colour_Deconvolution2(ImgR, ImgG, ImgB, 3, 0, 1);

% conversion to the appropriate variable
% grayscale
HEM_gray = mat2gray(Dye01_transmittance, [0, 255]);
DAB_gray = mat2gray(Dye02_transmittance, [0, 255]);
RES_gray = mat2gray(Dye03_transmittance, [0, 255]);

% real colors
HEM_real = ind2rgb(uint8(Dye01_transmittance), LUTdye01);
DAB_real = ind2rgb(uint8(Dye02_transmittance), LUTdye02);
RES_real = ind2rgb(uint8(Dye03_transmittance), LUTdye03);

% grayscale results
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(HEM_gray); title('Hematoxylin')
subplot(2,4,3); imshow(DAB_gray); title('DAB')
subplot(2,4,4); imshow(RES_gray); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(HEM_gray); title('Hematoxylin')
subplot(2,4,7); imhist(DAB_gray); title('DAB')
subplot(2,4,8); imhist(RES_gray); title('Residual')

% results real colors
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(HEM_real); title('Hematoxylin')
subplot(2,4,3); imshow(DAB_real); title('DAB')
subplot(2,4,4); imshow(RES_real); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(HEM_real); title('Hematoxylin')
subplot(2,4,7); imhist(DAB_real); title('DAB')
subplot(2,4,8); imhist(RES_real); title('Residual')

% deconvolution function 2
imageHDAB = my_deconvolution(image, 'd', 'A');

% results
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(imageHDAB(:,:,1)); title('Hematoxylin')
subplot(2,4,3); imshow(imageHDAB(:,:,2)); title('DAB')
subplot(2,4,4); imshow(imageHDAB(:,:,3)); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(imageHDAB(:,:,1)); title('Hematoxylin')
subplot(2,4,7); imhist(imageHDAB(:,:,2)); title('DAB')
subplot(2,4,8); imhist(imageHDAB(:,:,3)); title('Residual')

% choice of deconvolution method
HEM_deconv = HEM_gray; %imageHDAB(:,:,1);
DAB_deconv = DAB_gray; %imageHDAB(:,:,2);
RES_deconv = RES_gray; %imageHDAB(:,:,3);

%% 1D. INTENSITY VALUE CORRECTION
% image intensity value correction
% imadjust (stretching to 0:1; no saturation)
imadjust_HEM = imadjust(HEM_deconv); % better one
imadjust_DAB = imadjust(DAB_deconv);

% results
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(imadjust_HEM); title('Hematoxylin after imadjust')
subplot(2,4,3); imshow(imadjust_DAB); title('DAB after imadjust')
subplot(2,4,4); imshow(RES_deconv); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(imadjust_HEM); title('Hematoxylin after imadjust')
subplot(2,4,7); imhist(imadjust_DAB); title('DAB after imadjust')
subplot(2,4,8); imhist(RES_deconv); title('Residual')

% imadjust (stretching to 0:1; changing the saturation value)
imadjustTOP_HEM = imadjust(HEM_deconv, [0.01, 0.985]); % definitely the best
imadjustTOP_DAB = imadjust(DAB_deconv, [0.01, 0.985]); 

% results
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(imadjustTOP_HEM); title('Hem after imadjust with change in saturation')
subplot(2,4,3); imshow(imadjustTOP_DAB); title('DAB after imadjust with change in saturation')
subplot(2,4,4); imshow(RES_deconv); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(imadjustTOP_HEM); title('Hem after imadjust with change in saturation')
subplot(2,4,7); imhist(imadjustTOP_DAB); title('DAB after imadjust with change in saturation')
subplot(2,4,8); imhist(RES_deconv); title('Residual')

% adapthisteq - CLAHE algorithm (stretching and normalization)
adapthisteq_HEM = adapthisteq(HEM_deconv);
adapthisteq_DAB = adapthisteq(DAB_deconv);

% results
figure()
subplot(2,4,1); imshow(image); title('Original')
subplot(2,4,2); imshow(adapthisteq_HEM); title('Hematoxylin after adapthisteq')
subplot(2,4,3); imshow(adapthisteq_DAB); title('DAB after adapthisteq')
subplot(2,4,4); imshow(RES_deconv); title('Residual')

% histograms
subplot(2,4,5); imhist(image); title('Original')
subplot(2,4,6); imhist(adapthisteq_HEM); title('Hematoxylin after adapthisteq')
subplot(2,4,7); imhist(adapthisteq_DAB); title('DAB after adapthisteq')
subplot(2,4,8); imhist(RES_deconv); title('Residual')


% image selection based on deconvolution method and results from intensity value correction 
image_HEM = HEM_gray; %imadjust_HEM
image_DAB = DAB_gray; %imadjust_DAB
image_RES = RES_gray; %imageHDAB(:,:,3)

%% 2. FINDING THE BACKGROUND AREA
    
% note the correct format of the image:
% for method 'a' -> imageHDAB(:,:,1)
[FilledImg, BBox, fill] = my_background_plastic(imageHDAB(:,:,1), 'a', plastic_det);

%% 3. IMAGE ROTATION + 4. IMAGE CROPPING
    
% images to be processed - mask is first in cell array
% images can be multiple, mask determines where to cropimage_set{1} = fill;
image_set{2} = image;
image_set{3} = image_HEM;
image_set{4} = image_DAB;
image_set{5} = image_RES;

[image_new_a, kat_a] = my_rotation(image_set, BBox, 'k', 'h');
[image_new_e, kat_e] = my_rotation(image_set, BBox, 'h', 'h');

 % elimination of artifacts outside the mask area with white color
for s = 3:5
    image_new_a{s}(~image_new_a{1}) = 1;
    image_new_e{s}(~image_new_e{1}) = 1;
end

% sum of pixels results
figure()
subplot(2,3,1); imshow(image); title('Original before')
subplot(2,3,2); imshow(image_new_a{1}); title("Mask | Sum of pixels | Angle = " + num2str(round(kat_a)))
subplot(2,3,3); imshow(image_new_a{2}); title("Original after")
subplot(2,3,4); imshow(image_new_a{3}); title("Hem")
subplot(2,3,5); imshow(image_new_a{4}); title("DAB")
subplot(2,3,6); imshow(image_new_a{5}); title("Residual")

% Hough results
figure()
subplot(2,3,1); imshow(image); title('Original before')
subplot(2,3,2); imshow(image_new_e{1}); title("Mask | Hough | Angle = " + num2str(round(kat_e)))
subplot(2,3,3); imshow(image_new_e{2}); title("Original after")
subplot(2,3,4); imshow(image_new_e{3}); title("Hem")
subplot(2,3,5); imshow(image_new_e{4}); title("DAB")
subplot(2,3,6); imshow(image_new_e{5}); title("Residual")

%% 5. SPECIFIC ANALYSIS OF THE IMAGE
    
% THRESHOLD - CHECK THE AVERAGE PROPORTION OF THE DAB IMAGE
srednia_proporcja = 100*sum(sum(image_new_a{4} < 0.6))/(size(image_new_a{4}, 1)*size(image_new_a{4}, 2));

if srednia_proporcja < 0.0629
    disp("Image contains too small of DAB - probably it is a negative sample");
elseif srednia_proporcja >= 0.0629
    % K-MEANS CLUSTERING
    % 'a' - k-means clustering without GMM
    % 'b' - clustering of k-means with GMM
    cluster_method = 'a';

    % for HEM image
    [mask_HEM, cluster_HEM, cluster_before_HEM, lab_pic_seg_HEM, j_HEM] = my_k_means(image_new_a{3}, cluster_method);
    
    % for DAB image
    [mask_DAB, cluster_DAB, cluster_before_DAB, lab_pic_seg_DAB, j_DAB] = my_k_means(image_new_a{4}, cluster_method);
    
    % results (application layout)
    figure()
    for i = 1:j_HEM
        subplot(j_HEM+2,1,1); imshow(image_new_a{2}); title("Original image")
        subplot(j_HEM+2,1,2); imshow(image_new_a{3}); title("HEM image")
        subplot(j_HEM+2,1,i+2); imshow(cluster_HEM{i}, []); title("The result of the cluster " + i + " for HEM image")
    end
    
    figure()
    for i = 1:j_DAB
        subplot(j_DAB+2,1,1); imshow(image_new_a{2}); title("Original image")
        subplot(j_DAB+2,1,2); imshow(image_new_a{4}); title("DAB image")
        subplot(j_DAB+2,1,i+2); imshow(cluster_DAB{i}, []); title("The result of the cluster " + i + " for DAB image")
    end 
    
    % outline the DAB areas from the selected cluster (last in order) on the original image
    BW_bound = imbinarize(cluster_DAB{j_DAB});

    [B_bound, L_bond, N_bound] = bwboundaries(BW_bound);
    
    figure(); imshow(image_new_a{2}); title("Percentage of DAB-stained tissue = ")
    hold on;
    for k = 1:length(B_bound)
       boundary = B_bound{k};
       if(k > N_bound)
         % external edges
         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
       else
         % internal edges 
         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
       end
    end
    hold off

    % without plastic
    %[B_bound, L_bond, N_bound] = bwboundaries(image_new_a{1});
    
    %figure(); imshow(image_new_a{2}); title("Image original with applied mask outlines without plastic")
    %hold on;
    %for k = 1:length(B_bound)
       %boundary = B_bound{k};
       %if(k > N_bound)
         % external edges
         %plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
       %else
         % internal edges 
         %plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
       %end
    %end
    %hold off 
    
    % histograms
    % clustering and graphing
    j_HEM_ = 1:(j_HEM*2);
    evens_HEM = j_HEM_(mod(j_HEM_,2)==0); %paired elements
    odds_HEM = j_HEM_(mod(j_HEM_,2)==1); %unpaired elements
    
    j_DAB_ = 1:(j_DAB*2);
    evens_DAB = j_DAB_(mod(j_DAB_,2)==0); %paired elements
    odds_DAB = j_DAB_(mod(j_DAB_,2)==1); %unpaired elements
    
    figure()
    for i = 1:j_HEM
        subplot(j_HEM,2,odds_HEM(i)); imshow(cluster_HEM{i}, []); title("The result of the cluster " + i + " for HEM image")
        subplot(j_HEM,2,evens_HEM(i)); imhist(cluster_HEM{i}); title("Histogram of the cluster " + i + " for HEM image")
    end 
    
    figure()
    for i = 1:j_DAB
        subplot(j_DAB,2,odds_DAB(i)); imshow(cluster_DAB{i}, []); title("The result of the cluster " + i + " for DAB image")
        subplot(j_DAB,2,evens_DAB(i)); imhist(cluster_DAB{i}); title("Histogram of the cluster " + i + " for DAB image")
    end 
    
    
    % CALCULATION OF % DAB OCCUPANCY ON THE IMAGE RELATIVE TO THE ENTIRE AREA OF THE TISSUE SLICE
    if cluster_method == 'a'
        DAB_cluster = sum(sum(cluster_before_DAB{j_DAB} < 1));
        tissue_mask = sum(sum(fill == 1));
        DAB_area = (DAB_cluster/tissue_mask) * 100;
        display(DAB_area)
    elseif cluster_method == 'b'
        DAB_cluster = sum(sum(cluster_before_DAB{j_DAB} < 255));
        tissue_mask = sum(sum(fill == 1));
        DAB_area = (DAB_cluster/tissue_mask) * 100;
        display(DAB_area)
    end 

    % DAB SEGMENTATION
    [average_intenisty_DAB, picture_grey_255] = my_intensity_DAB(image_new_a{3}, image_new_a{4});
    display(average_intenisty_DAB)
    
    % results
    figure()
    imshow(picture_grey_255); title('DAB stained tissue after applying mask to original image')
end 