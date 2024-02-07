function [cluster_DAB, j_DAB, image_new_a, DAB_area] = DAB_detection_system(image)
    %% 1. INITIAL COLOR NORMALIZATION + DECONVOLUTION
    marker = 'c'; 
       
    % color normalization
    normal_Reinhard = my_norm(image, marker, 'A');
    
    % image choices based on results from color normalization
    image = uint8(255 * mat2gray(normal_Reinhard));
   
    % image division into channels for deconvolution
    image2 = double(image);
    ImgR = image2(:,:,1);
    ImgG = image2(:,:,2);
    ImgB = image2(:,:,3);
    
    % deconvolution function without color normalization
    [~, ~, ~, Dye01_transmittance, Dye02_transmittance, Dye03_transmittance, ~, ~, ~, ~] = Colour_Deconvolution2(ImgR, ImgG, ImgB, 3, 0, 1);
    
    % conversion to the appropriate variable
    % grayscale
    HEM_gray = mat2gray(Dye01_transmittance, [0, 255]);
    DAB_gray = mat2gray(Dye02_transmittance, [0, 255]);
    RES_gray = mat2gray(Dye03_transmittance, [0, 255]);
    
    % deconvolution function 2
    imageHDAB = my_deconvolution(image, 'd', 'A');
    
    % image selection based on deconvolution method and results from intensity value correction 
    image_HEM = HEM_gray; %imadjust_HEM
    image_DAB = DAB_gray; %imadjust_DAB
    image_RES = RES_gray; %imageHDAB(:,:,3)
    
    %% 2. FINDING THE BACKGROUND AREA
    
    % note the correct format of the image:
    % for method 'a' -> imageHDAB(:,:,1)
    [~, BBox, fill] = my_background(imageHDAB(:,:,1), 'a'); % creation of a mask 
    
    %% 3. IMAGE ROTATION + 4. IMAGE CROPPING
    
    % images to be processed - mask is first in cell array
    % images can be multiple, mask determines where to crop
    image_set{1} = fill;
    image_set{2} = image;
    image_set{3} = image_HEM;
    image_set{4} = image_DAB;
    image_set{5} = image_RES;
    
    [image_new_a, ~] = my_rotation(image_set, BBox, 'k', 'h');

    % elimination of artifacts outside the mask area with white color
    for s = 3:5
        image_new_a{s}(~image_new_a{1}) = 1;
    end

    %% 5. SPECIFIC ANALYSIS OF THE IMAGE
    
    % THRESHOLD - CHECK THE AVERAGE PROPORTION OF THE DAB IMAGE
    srednia_proporcja = 100*sum(sum(image_new_a{4} < 0.6))/(size(image_new_a{4}, 1)*size(image_new_a{4}, 2));

    if srednia_proporcja < 0.0629
        %disp("Image contains too small of DAB - probably it is a negative sample");
        DAB_area = 0.0629;
        cluster_DAB = 0; 
        j_DAB = 0; 
        image_new_a = 0;
        
    elseif srednia_proporcja >= 0.0629
        % K-MEANS CLUSTERING
        % 'a' - k-means clustering without GMM
        % 'b' - clustering of k-means with GMM
        cluster_method = 'a';

        % for DAB image
        [~, cluster_DAB, cluster_before_DAB, ~, j_DAB] = my_k_means(image_new_a{4}, cluster_method);

        % CALCULATION OF % DAB OCCUPANCY ON THE IMAGE RELATIVE TO THE ENTIRE AREA OF THE TISSUE SLICE
        if cluster_method == 'a'
            DAB_cluster = sum(sum(cluster_before_DAB{j_DAB} < 1));
            tissue_mask = sum(sum(fill == 1));
            DAB_area = (DAB_cluster/tissue_mask) * 100;
            %display(DAB_area)
        elseif cluster_method == 'b'
            DAB_cluster = sum(sum(cluster_before_DAB{j_DAB} < 255));
            tissue_mask = sum(sum(fill == 1));
            DAB_area = (DAB_cluster/tissue_mask) * 100;
            %display(DAB_area)
        end
    end
end  