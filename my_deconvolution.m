function [imageHDAB] = my_deconvolution(image, method, normali)
% Custom function to select the appropriate set of vectors necessary for the deconvolution process

    if method == 'a'
        %  https://github.com/jnkather/ColorDeconvolutionMatlab
        He = [ 0.6500286; 0.704031; 0.2860126 ];
        DAB = [ 0.26814753; 0.57031375; 0.77642715 ];
        Res = [ 0.7110272; 0.42318153; 0.5615672 ];

    elseif method == 'b'
        %  set of standard values for stain vectors (from python scikit; https://github.com/scikit-image/scikit-image/blob/main/skimage/color/colorconv.py)
        %  He = [0.651; 0.701; 0.29];
        %  DAB = [0.216; 0.801; 0.558]; %Eosin
        %  Res = [0.316; -0.598; 0.737];

        % from python scikit 
        He = [0.65; 0.70; 0.29];
        DAB = [0.27; 0.57; 0.78];
        Res = [0.07; 0.99; 0.11]; %Eosin

    elseif method == 'c'
        %  default Color Deconvolution Matrix proposed in Ruifork and Johnston
        %  He = [0.644211; 0.716556; 0.266844];
        %  DAB = [0.092789; 0.954111; 0.283111]; %Eosin
        %  Res = [-0.051733909968; -0.157623032505; 0.548160286737];

        % default Color Deconvolution Matrix proposed in Ruifork and Johnston
        He = [0.18; 0.20; 0.08]; 
        DAB = [0.10; 0.21; 0.29];
        Res = [0.01; 0.13; 0.01]; %Eosin

    elseif method == 'd'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.651; 0.701; 0.290]; 
        DAB = [0.269; 0.568; 0.778];
        Res = [0.633; -0.713; 0.302]; 
       
    elseif method == 'e'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.712; 0.673; 0.201]; 
        DAB = [0.477; 0.669; 0.570];
        Res = [0.583; -0.726; 0.364];
 
    elseif method == 'f'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.701; 0.676; 0.230]; 
        DAB = [0.522; 0.700; 0.487];
        Res = [0.542; -0.713; 0.444];
             
    elseif method == 'g'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.661; 0.620; 0.423]; 
        DAB = [0.329; 0.646; 0.689];
        Res = [0.369; -0.759; 0.536];
               
    elseif method == 'h'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.641; 0.646; 0.416]; 
        DAB = [0.498; 0.616; 0.611];
        Res = [0.572; -0.762; 0.302];
                 
    elseif method == 'i'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.710; 0.661; 0.241]; 
        DAB = [0.456; 0.642; 0.617];
        Res = [0.573; -0.742; 0.348];
                         
    elseif method == 'j'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.737; 0.631; 0.241]; 
        DAB = [0.524; 0.636; 0.566];
        Res = [0.535; -0.763; 0.362];

     elseif method == 'k'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.725; 0.586; 0.363]; 
        DAB = [0.481; 0.622; 0.618];
        Res = [0.391; -0.783; 0.484];

     elseif method == 'l'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.741; 0.567; 0.360]; 
        DAB = [0.463; 0.626; 0.627];
        Res = [0.340; -0.779; 0.527];

     elseif method == 'm'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.715; 0.648; 0.264]; 
        DAB = [0.478; 0.623; 0.619];
        Res = [0.567; -0.758; 0.324];

     elseif method == 'n'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.759; 0.610;  0.229]; 
        DAB = [0.468; 0.636; 0.614];
        Res = [0.488; -0.765; 0.420];  

     elseif method == 'o'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.726; 0.661; 0.188]; 
        DAB = [0.482; 0.646; 0.591];
        Res = [0.588; -0.740; 0.329];

     elseif method == 'p'
        %  https://github.com/qupath/qupath/wiki/Preprocessing
        %  https://www.imagescientist.com/brightfield-color-deconvolution
        He = [0.741; 0.606; 0.290]; 
        DAB = [0.445; 0.647; 0.619];
        Res = [0.432; -0.761; 0.485]; 

    elseif method == 'r'
        % funkcja Colour_Deconvolution2 from Fiji
        % https://github.com/landinig/IJ-Colour_Deconvolution2
        He = [1.20008420810656; -0.689387725766535;	0.636214232961718]; 
        DAB = [0.685203230666564; 0.0235969196559078; -0.710026796260320];
        Res = [-0.917768568624354; 1.50870556052618; 0.301816829168353]; 

    end

    if normali == 'A'
        % merging stain vectors with deconvolution matrix
        HDABtoRGB = [He/norm(He) DAB/norm(DAB) Res/norm(Res)]';
        RGBtoHDAB = inv(HDABtoRGB);
    
        % separate stains = performing color deconvolution
        imageHDAB = SeparateStains(image, RGBtoHDAB);

    elseif normali == 'B'
        imgRGB = double(image);

        % normalization
        p11 = He(1)/sqrt(He(1)^2+He(2)^2+He(3)^2); p12 = He(2)/sqrt(He(1)^2+He(2)^2+He(3)^2); p13 = He(3)/sqrt(He(1)^2+He(2)^2+He(3)^2);
        p21 = DAB(1)/sqrt(He(1)^2+He(2)^2+He(3)^2); p22 = DAB(2)/sqrt(He(1)^2+He(2)^2+He(3)^2); p23 = DAB(3)/sqrt(He(1)^2+He(2)^2+He(3)^2);
        p31 = Res(1)/sqrt(He(1)^2+He(2)^2+He(3)^2); p32 = Res(2)/sqrt(He(1)^2+He(2)^2+He(3)^2); p33 = Res(3)/sqrt(He(1)^2+He(2)^2+He(3)^2);

        % creation of a matrix from the values obtained
        HDAB = inv([p11 p12 p13; p21 p22 p23; p31 p32 p33]);

        % deconvolution
        imageHDAB = separation(imgRGB, HDAB);

    end 

end 