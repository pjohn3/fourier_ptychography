function objectRecover = reconstructImage(LEDgap, LEDheight, arraysize, wavelength, NA, spsize, imSeqLowRes)

    % "Set up" the LED grid. Where is each LED located in x, y plane?    
    % Consider a row-major flattened version of the LED matrix.
    % xlocation represents the x-location of the nth element in the flattened array.
    % ylocation represents the y-location of the nth element in the flattened array.

    xlocation = zeros(1, arraysize^2);
    ylocation = zeros(1, arraysize^2);
    for i=1:arraysize % from topleft to bottomright
        xlocation(1, 1 + arraysize *(i - 1):arraysize + arraysize*(i - 1)) = (-(arraysize - 1)/2:1:(arraysize - 1) / 2) * LEDgap;
        ylocation(1, 1 + arraysize *(i - 1):arraysize + arraysize*(i - 1)) =((arraysize - 1)/2 - (i - 1)) * LEDgap;
    end;
        
    % Compute the wavevector of the light emitted by each LED in the x and y directions.
    k0 = 2 * pi / wavelength;
    kx_relative = -sin(atan(xlocation/LEDheight));
    ky_relative = -sin(atan(ylocation/LEDheight));
    
    % Fs = 2*pi/spsize must be > 2 * max frequency in signal. So kmax for camera = pi/spsize.    
    kmax = pi/spsize;

    % Angular wavenumber components
    kx = k0 * kx_relative;
    ky = k0 * ky_relative;
    cutoffFrequency = NA * k0;

    % If we capture regions of higher spatial frequency with our LED arrays, 
    % our final reconstructed image must sample the spatial domain at a high enough 
    % frequency to incorporate the ones in our image. Solve this issue by adjusting 
    % the output image's pixel size (which sets a limit on our highest represented frequency) 
    % to be small enough to accommodate the highest spatial frequency represented by any of the 
    % low res images.

    % To match Nyquist constraints, dynamically compute psize here.
    % Fs = 2*pi/psize must be > 2 * max frequency in signal.

    psize = min(pi/max([max(kx), max(ky)]), spsize/1.5);

    % Generate the low-pass filtered images.
    [m1, n1] = size(imSeqLowRes(:,:,1));
    m = ceil((spsize / psize) * m1);
    n = ceil((spsize / psize) * n1);
    [kxm, kym] = meshgrid(-kmax:kmax / ((n1 - 1)/2):kmax, -kmax:kmax/((m1 - 1) / 2):kmax);
    CTF = ((kxm.^2+kym.^2)<cutoffFrequency^2); % Coherent Transfer function. 
    
    % With the above psize explanation, we ensure that an arbitrary kx and ky
    % value will scale to an appropriate index in range of (1, length of matrix)
    % corresponding to the center of the NA constraint. Call this scaled index
    % kxc, kyc.

    % But is that really what we want?

    % What if an image's kxc and kyc is somewhere near the edge of the image?
    % Then the circle (NA constraint) in the Fourier space might go out of bounds!
    % So, modify computations to scale kx, ky to be closer to the center so that
    % the NA constraint lies within the boundaries of the matrix. The
    % implication here is that kxc, kyc should fall within NA from the left, and NA
    % from the right of the image. Hence the subtraction of 2*NA (and an
    % additional 1, just in case) from the original range of values within which
    % kxc and kyc could fall.

    % Spatial frequency resolution (# radians per pixel in the final output image) 
    % = Angular Fs/N
    % = (2pi/psize)/N
    dkx = 2 * pi / (psize * (n - 2 * size(CTF, 1) - 1));
    dky = 2 * pi / (psize * (m - 2 * size(CTF, 2) - 1));
    
    % recover the high resolution image
    seq = gseq(arraysize); % define the order of recovery, we start from the center (the 113th image) to the edge of the spectrum (the 225th image)
    objectRecover = ones(m, n); % initial guess of the object
    objectRecoverFT = fftshift(fft2(objectRecover));
    loop = 5;

    for tt=1:loop
        for i3=1:arraysize^2
            i2=seq(i3);
            kxc = round((n+1)/2+kx(1,i2)/dkx);
            kyc = round((m+1)/2+ky(1,i2)/dky); 
            kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2); 
            kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
            lowResFT = (m1/m)^2 * objectRecoverFT(kyl:kyh,kxl:kxh).*CTF; 
            im_lowRes = ifft2(ifftshift(lowResFT));
            im_lowRes = (m/m1)^2 *imresize(imSeqLowRes(:,:,i2), [m1, n1]).*exp(1i.*angle( im_lowRes)); 
            lowResFT=fftshift(fft2(im_lowRes)).*CTF;
            objectRecoverFT(kyl:kyh,kxl:kxh)=(1-CTF).*objectRecoverFT(kyl:kyh,kxl:kxh) + lowResFT;
        end;
    end;
    objectRecover=ifft2(ifftshift(objectRecoverFT));    
end