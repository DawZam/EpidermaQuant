function [mask, cluster, cluster_before, lab_pic_seg, j] = my_k_means(lab_pic, method)
% Custom function for k-means clustering
    
    if method == 'a'
        % progress for the original cropped image: 
    
        % assignment of variables for further analysis
        %imgRGB = double(image);
    
        % change image type to lab
        %lab_tissue = rgb2lab(imgRGB);
        %lab_pic = lab_tissue(:,:,1:3);
        lab_pic = im2single(lab_pic);
    
        % number of colors (subjective n = 3)
        [n,m,~] = size(lab_pic); 
        X = reshape(lab_pic,n*m,[]);
        [X, ~, ~] = normInp(X); % when we have GMM then without this stage?; X is pixel_GMM_ext
    
        % options for k-means
        opts = statset('Display', 'final', 'MaxIter', 1000); 
            
        % details for clustering algorithm
        clust_al = @(X,K)(kmeans(X, K, 'emptyaction', 'singleton', 'replicate', 20, 'Options', opts)); 
    
        % creation of a clustering evaluation object containing the data used to evaluate the optimal number of data clusters
        eval = evalclusters(double(X), clust_al, 'DaviesBouldin', 'KList', (1:3)); % at a later time after testing can be set to (1:7)
        
        % number of clusters
        j = eval.OptimalK; 
    
        % grouping 5 times to avoid local minima
        lab_pic_seg = imsegkmeans(lab_pic, j, 'NumAttempts', 50, 'MaxIterations', 1000); % clustering of the object in 50 clusters 
            
        % marking each pixel in the image with a pixel label
        %figure()
        %imshow(lab_pic_seg,[]); title('Image marked by cluster index')
    
        % grouping and mask creation -> select the appropriate one based on the set criterion (setting according to the intensity of white elements)
        %cluster{i} = lab_pic .* double(mask{i});
        cluster_int = zeros(j, 1);
        cluster = cell(j, 1);
        mask = cluster;
        for i = 1:j
            mask{i} = lab_pic_seg == i;
            ind = double(mask{i}) == 0;
            cluster{i} = lab_pic;
            cluster{i}(ind) = 1; % or 255 depending on the lab_pic scale 
            cluster_int(i) = mean(mean(cluster{i}(mask{i}))); 
        end 
        
        % cluster sorting
        [~, temp] = sort(cluster_int, 1, "descend");
        cluster = cluster(temp);
        cluster_before = cluster; 
    
        %figure()
        %imshow(cluster{1}, []); title('The result of the selected mask')
    end 

    if method == 'b'
        % get image and scale to [0-255]
        lab_pic = lab_pic*255;
        
        % get number of pixels per color value excluding max value(white on a screen)
        x = 0:254;
        y_orig = histcounts(lab_pic(:), [x-0.5,254.5]);
        
        % smooth data
        y = smoothdata(y_orig, 'gaussian', 9);
        
        % run GMM
        [pp_est, mu_est, sig_est] = gaussian_mixture_vector(x, y, 20, false, true);
        
        % calculate conditional probabilities
        KS = length(pp_est);
        prob = zeros(length(x), KS);
        for i = 1:1:KS
            prob(:,i) = pp_est(i)*normpdf(x, mu_est(i), sig_est(i));
        end
        tot_prob = sum(prob, 2);
        for i = 1:1:KS
            prob(:,i) = prob(:,i)./tot_prob;
        end
        
        % find pixels with at least one count
        ind = find(y_orig>0);
        
        % create matrix with conditional probabilities from GMM instead of pixel values
        pixel_GMM_ext = [];
        pixel_orig = [];
        for a = 1:length(ind)
            pixel_GMM_ext = [pixel_GMM_ext; repmat(prob(ind(a),:), y_orig(ind(a)), 1)];
            pixel_orig = [pixel_orig; x(ind(a)) * ones(y_orig(ind(a)), 1)];
        end
        
        % transform cond. prob. to get better distribution of variables
        pixel_GMM_ext_mod = pixel_GMM_ext;
        for a = 1:size(pixel_GMM_ext_mod, 2)
            pixel_GMM_ext_mod(:,a) = log10(pixel_GMM_ext_mod(:,a));
            ind = isinf(pixel_GMM_ext_mod(:,a));
            pixel_GMM_ext_mod(ind,a) = min(pixel_GMM_ext_mod(~ind,a));
            pixel_GMM_ext_mod(:,a) = pixel_GMM_ext_mod(:,a) - min(pixel_GMM_ext_mod(:,a));
            pixel_GMM_ext_mod(:,a) = pixel_GMM_ext_mod(:,a)/max(pixel_GMM_ext_mod(:,a));
        end
        
        % normalization 
        X = pixel_GMM_ext_mod; 
    
        % options for k-means
        opts = statset('Display', 'final', 'MaxIter', 1000); 
            
        % details for clustering algorithm
        clust_al = @(X,K)(kmeans(X, K, 'emptyaction', 'singleton', 'replicate', 20, 'Options', opts)); 
    
        % creation of a clustering evaluation object containing the data used to evaluate the optimal number of data clusters
        eval = evalclusters(double(X), clust_al, 'DaviesBouldin', 'KList', (1:3)); % at a later time after testing can be set to (1:7)
        
        % number of clusters
        j = eval.OptimalK; 
    
        % grouping 5 times to avoid local minima
        lab_pic_seg = kmeans(X, j, 'emptyaction', 'singleton', 'replicate', 50, 'Options', opts); % clustering of the object in 50 clusters 
        
        % clustering result - moving clusters on the matrix
        cluster_result = lab_pic_seg; 
        unik_orig = unique(pixel_orig);

        matrix_results_klaster = zeros(size(lab_pic));
        for a = 1:length(unik_orig)
            index_pixel = find(pixel_orig == unik_orig(a)); % Find indices of cluster numbers in pixel_orig
            [ind_1, ind_2] = find(round(lab_pic) == unik_orig(a)); % Find where to write clustering results in matrix_results
            matrix_results_klaster(sub2ind(size(matrix_results_klaster), ind_1, ind_2)) = cluster_result(index_pixel(1)); % Write clustering results to the matrix at the corresponding positions
        end

        % grouping and creating a mask -> select the appropriate one based on the given criterion (setting according to the intensity of white elements)
        cluster_int = zeros(j, 1);
        cluster = cell(j, 1);
        mask = cluster;
        for i = 1:j
            mask{i} = matrix_results_klaster == i;
            ind = double(mask{i}) == 0;
            cluster{i} = lab_pic;
            cluster{i}(ind) = 255; % or 1 depending on the lab_pic scale 
            cluster_int(i) = mean(mean(cluster{i}(mask{i})));
        end 
        
        % cluster sorting
        [~, temp] = sort(cluster_int, 1, "descend");
        cluster_before = cluster(temp);
        
        % 'display range' setting - to make DAB structures visible on individual clusters
        cluster_ = cell(j, 1);
        for i = 1:j
            cluster_{i} = mat2gray(cluster_before{i});
        end
        cluster = cluster_;
    end
end 