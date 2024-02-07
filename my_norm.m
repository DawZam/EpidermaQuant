function [normal_image] = my_norm(image, picture, method)
% Custom function to select color normalization method

% Normalization with different methods (ready-made functions)
%normal_SCD = Norm(image, target, 'SCD');
%normal_SCDLeeds = Norm(image, target, 'SCDLeeds');
%normal_Macenko = Norm(image, target, 'Macenko');
    
    if picture == 'a'        
        % reference image for FLG marker (white background and brown pigment present)
        target = imread('TARGET_FLG.jpg');
    elseif picture == 'b'
        % reference image for marker K10 (white background and brown pigment present)
        target = imread('TARGET_K10.jpg');
    elseif picture == 'c'
        % reference image for Ki67 marker (white background and brown pigment present)
        target = imread('TARGET_Ki67.jpg');
    end 

    if method == 'A'
        % Reinhard color normalization
        normal_image = Norm(image, target, 'Reinhard');
    elseif method == 'B'
        % RGBHist color normalization
        normal_image = Norm(image, target, 'RGBHist');
    elseif method == 'C'
        % Macenko color normalization
        normal_image = Norm(image, target, 'Macenko');
    end 

end 