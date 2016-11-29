function objectRecover = runAlgorithm(dirName, LEDgap, LEDheight, arraysize, wavelength, NA, spsize)
    structFiles = dir(strcat(dirName, '/*.png'));
    fileNames = {structFiles.name};
    dimensions = size(imread(strcat(dirName, '/', char(fileNames(1)))));
    xlength = dimensions(1);
    ylength=  dimensions(1);
    imSeqLowRes = zeros(xlength, ylength, length(fileNames));
    for i = 1: length(fileNames)
        sampleImg = imread(strcat(dirName, '/', char(fileNames(i))));
        sampleImg(:,:,4) = [];
        imSeqLowRes(:, :, i) = rgb2gray(imresize(sampleImg, [1024 1024]));        
        %imshow(imSeqLowRes(:,:,i), []);
        %title(i);
        %figure;
    end    
    objectRecover = reconstructImage(LEDgap, LEDheight, arraysize, wavelength, NA, spsize, imSeqLowRes);
end