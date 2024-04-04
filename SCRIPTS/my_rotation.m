function [img_out, angle, img_out_] = my_rotation(image, BBox, method, pos) 
% Custom function to rotate an image. Image contains a cell array with a mask and consecutive images to be cropped
 
% For paramter 'method':
% 'k' method based on pixel sum
% 'h' method based on Hough transform

% For paramter 'pos':
% 'h' rotation to horizontal
% 'v' rotation to vertical
    
    if length(image) > 1
        for i = 2:length(image)
            image{i} = imcrop(image{i}, BBox);
        end
    end

     % Rotation angle calculation - sum of pixels or Hough transform
    if method == 'k'
        img = 1-image{1};
        angle = 0:0.1:180;
        maxx = zeros(length(angle), 1);
        for k = 1:length(angle)
            maxx(k) = max(sum(imrotate(img, angle(k)), 2));
        end
        [~, ind] = max(maxx);
        angle = angle(ind);
    elseif method == 'h'
        image_edge = edge(image{1},"canny");
        [~, dim2] = size(image_edge);
        [H, theta, rho] = hough(image_edge);
        image_hough_peaks  = houghpeaks(H,1,'threshold',ceil(0.95*max(H(:))));
        image_hough_lines = houghlines(image_edge,theta,rho,image_hough_peaks,"FillGap",0.8*dim2,"MinLength",30);
        angle = image_hough_lines(1).theta + 90;
    end

    % Final rotation - placement of finished image 
    if pos == 'h'
        n = length(image);
        img_out{i} = ones(1, n);
        img_out_{i} = ones(1, n);
        for i = 1:n            
            if i == 1
                img_out{i} = imrotate(image{i}, angle);
            elseif i == 2
                % convert default black background to white 
                img_out{i} = imrotate(image{i}, angle);
                img_out_{i} = ~imrotate(true(size(image{i})), angle);
                img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 255; %value for white color (scale for image 0-255)
            else 
                % convert default black background to white 
                img_out{i} = imrotate(image{i}, angle);
                img_out_{i} = ~imrotate(true(size(image{i})), angle);
                img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 1; %value for white color (scale for image 0-1)
            end
        end
    elseif pos == 'v'
        angle = angle-90;
        n = length(image);
        img_out{i} = ones(1, n);
        img_out_{i} = ones(1, n);
        for i = 1:n
            if i == 1
                img_out{i} = imrotate(image{i}, angle);
            elseif i == 2
                % convert default black background to white 
                img_out{i} = imrotate(image{i}, angle);
                img_out_{i} = ~imrotate(true(size(image{i})), angle);
                img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 255; %wartość dla koloru białego (skala dla obrazu 0-255)
            else 
                % convert default black background to white 
                img_out{i} = imrotate(image{i}, angle);
                img_out_{i} = ~imrotate(true(size(image{i})), angle);
                img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 1; %wartość dla koloru białego (skala dla obrazu 0-1)
            end
        end
    end

    regions = regionprops(img_out{1}, "Area", "FilledImage", "BoundingBox");
    [~, ind] = max([regions.Area]);
    BBox = regions(ind).BoundingBox;
    BBox(1) = BBox(1) + 0.5;
    BBox(2) =  BBox(2) + 0.5;
    BBox(3) = BBox(3) - 1;
    BBox(4) =  BBox(4) - 1;
    for i=1:length(img_out)
        img_out{i} = imcrop(img_out{i}, BBox);
    end
    
    % Check if the height is not greater than the width. If so - rotate by 90 (add 90 degrees) and repeat the rotation step 
    if method == 'k' && pos == 'h'
        [height, width] = size(img_out{i});
        if height < width
        elseif height > width
            angle = angle + 90;
            % Final rotation - laying out the finished image 
            n = length(image);
            img_out{i} = ones(1, n);
            img_out_{i} = ones(1, n);
            for i = 1:n            
                if i == 1
                    img_out{i} = imrotate(image{i}, angle);
                elseif i == 2
                    % convert default black background to white 
                    img_out{i} = imrotate(image{i}, angle);
                    img_out_{i} = ~imrotate(true(size(image{i})), angle);
                    img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 255; %value for white color (scale for image 0-255)
                else 
                    % convert default black background to white 
                    img_out{i} = imrotate(image{i}, angle);
                    img_out_{i} = ~imrotate(true(size(image{i})), angle);
                    img_out{i}(img_out_{i}&~imclearborder(img_out_{i})) = 1; %value for white color (scale for image 0-1)
                end
            end
    
            regions = regionprops(img_out{1}, "Area", "FilledImage", "BoundingBox");
            [~, ind] = max([regions.Area]);
            BBox = regions(ind).BoundingBox;
            BBox(1) = BBox(1) + 0.5;
            BBox(2) =  BBox(2) + 0.5;
            BBox(3) = BBox(3) - 1;
            BBox(4) =  BBox(4) - 1;
            for i=1:length(img_out)
                img_out{i} = imcrop(img_out{i}, BBox);
            end
        end
    end
end