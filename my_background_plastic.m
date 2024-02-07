function [FilledImg, BBox, fill] = my_background_plastic(image, method, plastik_mask)
% Custom function to choose the appropriate algorithm for finding the background area

    if method == 'a'
        % GRAYSCALE HEMATOXYLIN CHANNEL
        backg_1 = double(image);
        backg_1 = imbinarize(backg_1);

        % morphological operations
        se = strel('disk', 15);

        % opening
        opened = imopen(backg_1, se);
    
        % closing
        closed = imclose(opened, se);
        
        % creation of a mask
        closed = imcomplement(closed); %color inversion
        fill_img = imfill(closed, "holes");

        % plastic removing
        fill_img = fill_img-(plastik_mask == 1);
       
        % for custom images where more than one tissue fragment appears
        regions = regionprops(fill_img, "Area", "FilledImage", "BoundingBox");
        [~, ind] = max([regions.Area]);
        FilledImg = [regions(ind).FilledImage];
        BBox = regions(ind).BoundingBox;
        BBox(1) = BBox(1) + 0.5;
        BBox(2) =  BBox(2) + 0.5;
        BBox(3) = BBox(3) - 1;
        BBox(4) =  BBox(4) - 1;
        cropped = imcrop(image, BBox);
        ready = logical(cropped.* FilledImg);
        fill_ = ready;
        
        % morphological reoperation
        fill_close = imclose(fill_, se);
        fill_open = imopen(fill_close, se);
        fill = imfill(fill_open, "holes");
        
    end

end