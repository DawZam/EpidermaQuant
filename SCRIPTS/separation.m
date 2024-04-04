function imageOut = separation(imageRGB, MatrixHDAB)
    % dekonvolution 
    imageOut = reshape(-log(imageRGB),[],3) * MatrixHDAB;
    imageOut = reshape(imageOut, size(imageRGB));

    % normalization
    imageOut = normalization(imageOut);
    
end