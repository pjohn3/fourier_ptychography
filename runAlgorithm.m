function objectRecover = runAlgorithm(dirName, LEDgap, LEDheight, arraysize, wavelength, NA, spsize)
    structFiles = dir(strcat(dirName, '/*.png'));
    fileNames = {structFiles.name};
    dimensions = size(imread(strcat(dirName, '/', char(fileNames(1)))));
    xlength = dimensions(1);
    ylength=  dimensions(2);
    imSeqLowRes = zeros(xlength, ylength, length(fileNames));    
    for i = 1: length(fileNames)
        sampleImg = imread(strcat(dirName, '/', char(fileNames(i))));
        sampleImg(:,:,4) = [];
        imSeqLowRes(:, :, i) = rgb2gray(sampleImg);      
    end    
    objectRecover = reconstructImage(LEDgap, LEDheight, arraysize, wavelength, NA, spsize, imSeqLowRes);
    
    figure;
    title('Recovered magnitude');    
    imshow(abs(objectRecover), []);
    
    figure;
    title('Recovered phase');    
    imshow(angle(objectRecover), []);
    
    figure;
    title('Recovered FT');    
    imshow(log(fftshift(fft2(objectRecover))), []);
end