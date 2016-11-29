function imSeqLowRes = illuminateImages(LEDgap, LEDheight, arraysize, wavelength, NA, spsize, object)
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

    psize = pi/max([max(kx), max(ky)]);

    % Generate the low-pass filtered images.
    [m, n] = size(object); % image size of the high resolution object.
    m1 = floor(m / (spsize / psize));
    n1 = floor(n / (spsize / psize)); % image size of the final output.
    [kxm, kym] = meshgrid(-kmax:kmax / ((n1 - 1)/2):kmax, -kmax:kmax/((n1 - 1) / 2):kmax);
    CTF = ((kxm.^2+kym.^2)<cutoffFrequency^2); % Coherent Transfer function.
    imSeqLowRes = zeros(m1, n1, arraysize^2); % output low-res image sequence    
    
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

    objectFT = fftshift(fft2(object));
    for tt = 1:arraysize^2
        kxc = round((n+1) / 2 + kx(1, tt)/dkx);
        kyc = round((m+1) / 2 + ky(1, tt)/dky);
        kyl = round(kyc - (m1-1)/2);kyh = round(kyc+(m1-1)/2);
        kxl = round(kxc - (n1-1)/2);kxh = round(kxc+(n1-1)/2);
        imSeqLowFT = (m1/m)^2 * objectFT(kyl:kyh, kxl:kxh).*CTF;
        imSeqLowRes(:,:,tt) = abs(ifft2(ifftshift(imSeqLowFT)));
    end;

end