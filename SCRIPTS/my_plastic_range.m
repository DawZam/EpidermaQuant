function [fill_img_reg] = my_plastic_range(image, left, right, thresh)
% Custom function for plastic detection based on a preset color range    
    
    % channel image division
    imgR = image(:,:,1);
    imgG = image(:,:,2);
    imgB = image(:,:,3);
    
    % results for RGB image channels
    %figure()
    %subplot(2,3,1); imshow(image(:,:,1)); title('R channel - slice and plastic') 
    %subplot(2,3,2); imshow(image(:,:,2)); title('G channel - slice and plastic')
    %subplot(2,3,3); imshow(image(:,:,3)); title('B channel - slice and plastic')
    %hold on

    % appropriate treshold
    ind = image(:,:,1) < left | image(:,:,1) > right; 
    imgR(ind) = 255;
    imgG(ind) = 255;
    imgB(ind) = 255;
    
    % re-create an image from each channel
    image_pl = image;
    image_pl(:,:,1) = imgR; % only this is needed further for morphological operations 
    image_pl(:,:,2) = imgG;
    image_pl(:,:,3) = imgB;
    
    % RGB results for just plastic 
    %subplot(2,3,4); imshow(image_pl(:,:,1)); title('Channel R - plastic') 
    %subplot(2,3,5); imshow(image_pl(:,:,2)); title('Channel G - plastic')
    %subplot(2,3,6); imshow(image_pl(:,:,3)); title('Channel B - plastic')
    %hold off
    
    % morphological operations on the R channel (as in the case of the mask - removal of holes + filtering only after the single largest area)
    backg_1 = image_pl(:,:,1);
    backg_1 = imbinarize(backg_1);
    
    % morphological operations
    se = strel('disk', 15);
    
    % opening
    opened = imopen(backg_1, se);
    
    % closing
    closed = imclose(opened, se);
    
    % creation of a mask
    closed = imcomplement(closed); % color inversion
    fill_img = imfill(closed, "holes");

    % regionprops to remove oval artifacts
    regions = regionprops(fill_img, "Area", "FilledImage", "BoundingBox", "Extent", "Circularity");
    stats = regionprops('table', fill_img, 'BoundingBox', "Extent");
    
    len = length([regions.Extent]);
    LenWdRatio = zeros(1, len);
    isGreater = zeros(1, len);
    for i = 1:len
        stats.LenWdRatio = stats.BoundingBox(:,3) ./ stats.BoundingBox(:,4);
        LenWdRatio(:,i) = [regions(i).Extent]; 
        isGreater(:,i) = (LenWdRatio(i) > thresh);
    end
    
    figure()
    imshow(fill_img); 
    axis on
    hold on
    % plot red rectangles around accepted objects
    arrayfun(@(i)rectangle('Position', stats.BoundingBox(i,:), 'EdgeColor','r'), find(isGreater)); 
    % plot yellow rectangles around rejected objects
    arrayfun(@(i)rectangle('Position', stats.BoundingBox(i,:), 'EdgeColor','y'), find(~isGreater));
    hold off
    
    fill_img_reg = fill_img;

    for i = 1:len
        if LenWdRatio(i) > thresh
            x_1 = regions(i).BoundingBox(1);
            x_3 = regions(i).BoundingBox(1)+regions(i).BoundingBox(3);
            x = [x_1 x_3 x_3 x_1];
            
            y_2 = regions(i).BoundingBox(2);
            y_4 = regions(i).BoundingBox(2)+regions(i).BoundingBox(4);
            y = [y_2 y_2 y_4 y_4];
            
            fill_img_reg = regionfill(double(fill_img_reg),x,y);
            fill_img_reg = logical(fill_img_reg);
        end
    end
    
    %figure();
    %imshow(fill_img_reg); 

end 