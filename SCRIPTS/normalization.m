function PictureOut = normalization(imageIn)

	PictureOut = imageIn;
	for i=1:size(imageIn,3)
		Picture = PictureOut(:,:,i);

        % min max nirmalization
		PictureOut(:,:,i) = (PictureOut(:,:,i)-min(Picture(:)))/(max(Picture(:)-min(Picture(:))));
		
        % invert image
        PictureOut(:,:,i) = 1 - PictureOut(:,:,i);

	end

end