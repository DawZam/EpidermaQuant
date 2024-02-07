clc
clear variables
close all

% image choice
[I, path] = uigetfile('*.jpg', 'Choose one image');
str = strcat(path, I);
image = imread(str);

% function to detect DAB areas on immunohistochemical images
[cluster_DAB, j_DAB, image_new_a, DAB_area] = DAB_detection_system(image);

if cluster_DAB == 0
    % result
    display(DAB_area)
else
    % result (application)
    figure()
    for i = 1:j_DAB
        subplot(j_DAB+2,1,1); imshow(image_new_a{2}); title("Original image")
        subplot(j_DAB+2,1,2); imshow(image_new_a{4}); title("DAB image")
        subplot(j_DAB+2,1,i+2); imshow(cluster_DAB{i}, []); title("The result of the cluster " + i + " for DAB image")
    end 
    
    % outline application of DAB areas from the selected cluster (last in order) to the original image 
    BW_bound = imbinarize(cluster_DAB{j_DAB});
    
    [B_bound, L_bond, N_bound] = bwboundaries(BW_bound);
        
    figure(); imshow(image_new_a{2}); title("Original image with applied contours of DAB areas based on clustering results")
    hold on;
    for k = 1:length(B_bound)
        boundary = B_bound{k};
        if(k > N_bound)
            % external edges
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
        else
            % internal edges 
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
        end
    end
    hold off
    
    % histograms
    figure()
    for i = 1:j_DAB
        subplot(j_DAB,2,odds_DAB(i)); imshow(cluster_DAB{i}, []); title("The result of the cluster " + i + " for DAB image")
        subplot(j_DAB,2,evens_DAB(i)); imhist(cluster_DAB{i}); title("Cluster histogram " + i + " for DAB image")
    end 
    
    % results
    display(DAB_area)
end 