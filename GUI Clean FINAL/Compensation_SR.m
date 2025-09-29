%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:-->        CompensationSR                                 
%                                                                                                               
%                                                                                  
% Description:-->  This file is the contains the functions that are needed for the 
%                  OPEN_SR_QPI_APP. Some extra functions that are mainly for the
%                  Vortex?Legendre compensation are inlcuded in separate
%                  files. 
%                                                                                  
% Authors:-->      Sofia Obando-Vasquez(1), Ana Doblas(1) 
%                  (1) Department of Electrical and Computer Engineering, University 
%                      of Massachusetts Dartmouth
% Email:-->        sobandovasquez@umassd.edu                                         
% Date:-->         9/29/2025                                                      
% version 1.0 (2025)                                                               
 
classdef Compensation_SR
    methods(Static)

        function [holo,M,N,m,n] = holo_read(filename)
            % Inputs:
            %   filename  -> MxNxC numeric array 
            % Outputs:
            %   holo      -> NxM matrix with the first channel of the image.
            %                Size: N rows (height), M columns (width).
            %   M         -> Scalar, number of columns (image width).
            %   N         -> Scalar, number of rows (image height).
            %   m, n      -> NxM coordinate grids centered at zero, suitable for frequency-
            %                domain work and spatial indexing.
            holo = filename(:,:,1);
            [N,M] = size(holo);
            [n,m] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
        end


        function [holo_filter,HOLO_FT,circmask] = spatialFilter(holo,M,N, region)
            % Inputs:
            %   holo    -> NxMxC numeric array with a hologram.
            %   M       -> Scalar, number of columns (image width).
            %   N       -> Scalar, number of rows (image height).
            %   region  -> Scalar flag selecting the spectrum quadrant for the +1 order seed:   
            % Outputs:
            %   holo_filter -> NxM complex array; inverse-FFT of the masked spectrum.
            %                  This is the spatial-domain hologram after spectral filtering.
            %   HOLO_FT     -> NxM complex array; frequency-domain hologram after applying:
            %                    (i) quadrant rectangular pre-mask, and
            %                    (ii) adaptive circular binary mask. 
            %   circmask    -> NxM logical/double binary mask (0/1) created based on an adaptive threshold.
          
            fft_holo_1 = fftshift(fft2(fftshift(holo(:,:,1))));       
            % Initialize a filter with zeros
            filter = zeros(N,M);
            % Create a filter mask for the desired region
            if region==1
                filter(1:round(N/2-(N*0.15)),round(M/2+(M*0.15)):M) = 1; % 1nd quadrant
            elseif region==2
                filter(1:round(N/2-(N*0.15)),1:round(M/2-(M*0.15))) = 1;  % 2nd quadrant
            elseif region==3
                filter(round(N/2+(N*0.15)):N,1:round(M/2-(M*0.15))) = 1; % 3nd quadrant
            else
                filter(round(N/2+(N*0.15)):N,round(M/2+(M*0.15)):M) = 1; % 4nd quadrant
            end

            ft_filtered_holo_1 = fft_holo_1 .* filter;
            filtered_spect1 = log(abs(ft_filtered_holo_1).^2);
            filtered_spect1(isinf(filtered_spect1)) = 0;
            maskSize = 60;

            if region==1
                mean_corner = mean2(filtered_spect1(1:maskSize,M-maskSize:M));
                std_corner = std2(filtered_spect1(1:maskSize,M-maskSize:M));
            elseif region==2
                mean_corner = mean2(filtered_spect1(1:maskSize,1:maskSize));
                std_corner = std2(filtered_spect1(1:maskSize,1:maskSize));
            elseif region==3
                mean_corner = mean2(filtered_spect1(N-maskSize:N,1:maskSize));
                std_corner = std2(filtered_spect1(N-maskSize:N,1:maskSize));
            else
                mean_corner = mean2(filtered_spect1(N-maskSize:N,M-maskSize:M));
                std_corner = std2(filtered_spect1(N-maskSize:N,M-maskSize:M));
            end

            TH1 = mean_corner + 0.8*std_corner;
            circ = Compensation_SR.circular_filter(TH1, filtered_spect1,N);
            Filter_Spectrum_CirckMask = ft_filtered_holo_1.*circ;
            holo_filter = IFT(Filter_Spectrum_CirckMask);
            HOLO_FT = Filter_Spectrum_CirckMask;
            circmask = circ; 
        end

        function circ = circular_filter(TH, Filtered_spectrum, N)
            % Inputs:
            %   TH                -> Scalar threshold used to binarize the (log-)spectrum.
            %   Filtered_spectrum -> Numeric 2D array from which the mask is derived.
            %   N                 -> Scalar image size used to build the mask. The output mask
            %                        is N-by-N.
            % Outputs:
            %   circ -> N-by-N logical/double binary mask.
            
            BW = imbinarize(Filtered_spectrum, TH);
            BW = bwareafilt(BW, 1);
            BW_filled = imfill(BW, 'holes');
            edges = edge(BW_filled, 'Canny');
            [y,x] = find(edges);
            A = [-2*x, -2*y, ones(size(x))];
            b = -(x.^2 + y.^2);
            params = A\b;
            center(1) = params(1);
            center(2) = params(2);
            radius = ceil(sqrt(center(1)^2 + center(2)^2 - params(3)));
            circ = ones(N);
            for r = 1:N
                for p = 1:N
                    if sqrt((r-center(2))^2+(p-center(1))^2)>radius
                        circ(r,p) = 0;
                    end
                end
            end

        end

        function [fx_max, fy_max] = findMax(HOLO_FT, M,N) 
            % Only for Structure illumination
            % Inputs:
            %   HOLO_FT -> NxM complex (or real) 2-D spectrum of the hologram.
            %   M       -> Scalar, number of columns of the spectrum (image width).
            %   N       -> Scalar, number of rows of the spectrum (image height).
            % Outputs:
            %   fx_max  -> 1x2 vector with the x/column indices (in pixels, MATLAB 1-based)
            %              of the two highest peaks in |HOLO_FT|.
            %              fx_max = [x_firstPeak, x_secondPeak]
            %
            %   fy_max  -> 1x2 vector with the y/row indices (in pixels, MATLAB 1-based)
            %              of the two highest peaks in |HOLO_FT|.
            %              fy_max = [y_firstPeak, y_secondPeak] 
                    
            maxValue_1 = max(max(abs(HOLO_FT)));
            [fy_max_1 ,fx_max_1] = find(abs(HOLO_FT) == maxValue_1);
            mask = ones(M,N);
            mask(fy_max_1 - 10:fy_max_1 + 10,fx_max_1 - 10:fx_max_1 + 10)= 0;
            fx_max_L = fx_max_1;
            fy_max_L = fy_max_1;
            HOLO_FT = HOLO_FT .* mask;

            maxValue_1 = max(max(abs(HOLO_FT)));
            [fy_max_1, fx_max_1] = find(abs(HOLO_FT) == maxValue_1);
            fx_max_D = fx_max_1(1);
            fy_max_D = fy_max_1(1);

            fy_max = [fx_max_L,fx_max_D];
            fx_max = [fy_max_L,fy_max_D];
        end 

        function Gplus_demod = minimaztion_heuristic(holo_FT, fx_max, fy_max, askFn)
            % Inputs:
            %   holo_FT -> 2-D complex spectrum (size: Ny x Nx), 
            %   fx_max  -> Scalar or [1xK] integer column index/indices (1-based) of the
            %              +1 order peak(s) used by the cost function for seeding.
            %   fy_max  -> Scalar or [1xK] integer row index/indices (1-based) matching fx_max.
            %   askFn   -> (optional) function handle for an interactive yes/no prompt:
            %                
            % Outputs:
            %   Gplus_demod -> 2-D complex array (size ~ size(holo_FT,1) x size(holo_FT,2)),
            %                  the demodulated +1 order field. This is the result after the final
            %                  accepted iteration (i.e., when askFn returns false).
            if nargin < 4 || isempty(askFn)
                % fallback if called from scripts: use the standalone function if you have it
                askFn = @(img,prompt) askRepeatWithPreview(img, prompt);
            end

            repeat_minimization = true;
            while repeat_minimization
                % test_theta = randi([0 360])*pi/180;
                lb = 0;
                ub = 360*pi/180;
                options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, ...
                    'MaxIterations', 500, 'FunctionTolerance', 1e-8);
                cf_particleswarm = @(test_theta) Compensation_SR.costFunction_SIDHM(test_theta, holo_FT, fx_max, fy_max);

                [theta, ~] = particleswarm(@(params) Compensation_SR.costFunction_SIDHM(params, holo_FT, fx_max, fy_max), 1, lb, ub, options);
                [Gdemod] = Compensation_SR.demComp2SIDHM(theta, holo_FT);
                Gplus_demod = Gdemod(:,: ,1);
                imgPreview = log(abs(Gplus_demod).^2);
                repeat_minimization = askFn(imgPreview, 'Repeat the minimization?');

            end
            close all
        end

        function Gminus_demod = minimaztion_heuriticII(holo_FT, fx_max, fy_max, askFn)
            % Inputs:
            %   holo_FT -> 2-D complex spectrum (Ny x Nx). 
            %   fx_max  -> Scalar or vector of column indices (1-based) for +/-1 order seeds
            %              used by the cost function (must align with fy_max).
            %   fy_max  -> Scalar or vector of row indices (1-based) for +/-1 order seeds
            %              used by the cost function (must align with fx_max).
            %   askFn   -> (optional) function handle to query whether to repeat:
            %                repeat = askFn(imgPreview, promptString)
            %              
            % Output:
            %   Gminus_demod -> 2-D complex array (Ny x Nx): demodulated –1 order field,
            %                   taken as the second slice of demodulation output:
            %                     Gdemod(:,:,2)
            %                   This corresponds to the final accepted iteration (last loop
            %                   before askFn returns false).

            if nargin < 4 || isempty(askFn)
                % fallback if called from scripts: use the standalone function if you have it
                askFn = @(img,prompt) askRepeatWithPreview(img, prompt);
            end
            repeat_minimization = true;
            while repeat_minimization
                test_theta = randi([0 360])*pi/180;
                lb = 0;
                ub = 360*pi/180;
                options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, ...
                    'MaxIterations', 500, 'FunctionTolerance', 1e-8);
                cf_particleswarm = Compensation_SR.costFunction_SIDHMII(test_theta,holo_FT,fx_max,fy_max);

                [theta, ~] = particleswarm(@(params) Compensation_SR.costFunction_SIDHMII(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
                [Gdemod] =  Compensation_SR.demComp2SIDHM(theta, holo_FT);
                Gminus_demod = Gdemod(:,:,2);
                imgPreview = log(abs(Gminus_demod).^2);
                repeat_minimization = askFn(imgPreview, 'Repeat the minimization?');
            end
            close all
        end

        function [cf] = costFunction_SIDHM(theta, FTHolo, fx_max,fy_max)
            % Inputs:
            %   theta    -> Scalar (in radians). Demodulation angle to be evaluated.
            %   FTHolo   -> 2-D complex array (Ny x Nx). Fourier-domain hologram to demodulate.
            %   fx_max   -> [1x2] vector of column indices (x, 1-based) for the two peaks.
            %   fy_max   -> [1x2] vector of row    indices (y, 1-based) for the two peaks.
            %
            % Output:
            %   cf       -> Scalar in [0,1]. Normalized magnitude ratio at peak #1 over the
            %               sum at peaks #1 and #2 after demodulation with theta.

            [Dtemp] = Compensation_SR.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,1);
            cf =  abs(Dplus(fx_max(1),fy_max(1))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end

        function [cf] = costFunction_SIDHMII(theta, FTHolo, fx_max,fy_max)
            % Inputs:
            %   theta   -> Scalar (radians). Demodulation angle to evaluate.
            %   FTHolo  -> 2-D complex array (Ny x Nx). Fourier-domain hologram to demodulate.
            %   fx_max  -> [1x2] vector with column indices (x, 1-based) of candidate peaks.
            %   fy_max  -> [1x2] vector with row    indices (y, 1-based) of candidate peaks.
            %
            % Output:
            %   cf      -> Scalar in [0,1]. Normalized magnitude at the second location over
            %              the sum of magnitudes at both locations, after demodulation with
            %              the given theta.
            
            [Dtemp] = Compensation_SR.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,2);
            cf =  abs(Dplus(fx_max(2),fy_max(2))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end

        function [D] = demComp2SIDHM(theta, H)

            % Inputs:
            %   theta -> Scalar (radians). Relative phase shift between the two SI/DHM channels.
            %   H     -> X-by-Y-by-no complex array. Frequency-domain data with no channels.
            %            Expected no = 2, where:
            %              H(:,:,1) = first measurement (phase 0)
            %              H(:,:,2) = second measurement (phase theta)
            %
            % Outputs:
            %   D     -> X-by-Y-by-no complex array. Demodulated components obtained by
            %            linearly combining the two input channels:
            %              D(:,:,1) = "+1" component (first demodulated field)
            %              D(:,:,2) = "–1" component (second demodulated field)

            [X, Y, no] = size(H);
            D = zeros(X,Y,no);
            M = 1/2*[exp(1i*0) exp(-1i*0);exp(1i*theta) exp(-1i*theta)];
            Minv = pinv(M);
            D(:,:,1) = Minv(1,1).*H(:,:,1) + Minv(1,2).*H(:,:,2);
            D(:,:,2) = Minv(2,1).*H(:,:,1) + Minv(2,2).*H(:,:,2);
        end

        function gplus2 = Vortex_for_SR(Gplus_demod, M,N,k, m,n, lambda, dxy, max_order)
            % Source:
            %   Implemented following the method described in: 
            %   Ortega, K., Restrepo, R., Padilla-Vivanco, A., Castaneda, R., 
            %   Doblas, A., & Trujillo, C. (2025). 
            %   Intricate quantitative phase imaging via vortex-Legendre 
            %   high-order phase compensation.  
            %   Optics and Lasers in Engineering, 195, 109318.
            %
            % Inputs:
            %   Gplus_demod -> NxM complex array. Demodulated +1 order spectrum (frequency domain).
            %   M, N        -> Scalars. Image width (cols) and height (rows), respectively.
            %   k           -> Scalar. Wavenumber (2*pi/lambda) or provided externally.
            %   m, n        -> NxM coordinate grids (rows=y and cols=x, centered), typically from meshgrid.
            %   lambda      -> Scalar. Wavelength (same units as dxy).
            %   dxy         -> Scalar. Pixel size in object/image plane (same units as lambda).
            %   max_order   -> Scalar integer >= 2. Highest Legendre order used for aberration fitting.
            %
            % Output:
            %   gplus2      -> NxM complex array. Compensated +1 order field after:
            %                   (i) vortex-based tilt removal
            %                   (ii) Legendre-based phase aberration compensation
            %                   (includes piston correction if enabled).
            
            fx_0 = M/2;
            fy_0 = N/2;
            hologram = IFT(Gplus_demod);
            % Algorithm parameters to change
            medFilt = 3; % Change between 3 or 1
            NoPistonCompensation = false; % "False" to enable piston compensation.
            UsePCA = true; % "True" to enable compensation using PCA.
            limit= 500/2; % Spatial frequency support region region for calculating the Legendre coefficients.
            holo_filtered = hologram;

            [fxOverMax,fyOverMax] = Compensation_SR.maxVortex(hologram);
            % Vortex subpixel frequencies calculation
            ft_holo = log(abs(fftshift(fft2(fftshift(hologram)))).^2);
            field = medfilt2(ft_holo, [medFilt, medFilt], 'symmetric');
            [positions] = vortexCompensation(field, fxOverMax, fyOverMax);

            % Tilt compensation aberration using optical vortex
            [ref_wave] = Compensation_SR.reference_wave...
                (M,N,m,n,lambda,dxy,positions(1),positions(2),k,fx_0,fy_0);
            field_compensate = ref_wave.*holo_filtered;

            [~, Legendre_Coefficients] = ...
                LegendreCompensation(field_compensate, limit, NoPistonCompensation, UsePCA);

            % Vortex + piston
            if ~NoPistonCompensation
                PistonCorrtection = ones(size(field_compensate)).*Legendre_Coefficients(1);
                vortexPistonCorrection = field_compensate.*(exp(1i.*PistonCorrtection));
            end
            gridSize = size(angle(field_compensate),1);
            [X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), ...
                -1:(2 / gridSize):(1 - 2 / gridSize));
            dA=(2 / gridSize) ^ 2;

            order = 2:max_order;
            [polynomials] = squareLegendrefitting(order, X, Y);
            Legendres = reshape(polynomials, [size(polynomials, 1)*size(polynomials, 2) ...
                size(polynomials, 3)]);
            zProds = Legendres.'* Legendres * dA;
            zNorm = bsxfun(@rdivide, Legendres, sqrt(diag(zProds).'));
            Legendres = (bsxfun(@times, (ones(size(order))').', zNorm));
            Legendres_norm_const =sum(Legendres.^2,1)*dA;

            WavefrontReconstructed_Vect = sum(repmat(Legendre_Coefficients(2:size(order,2)+1)./ ...
                sqrt(Legendres_norm_const(1:size(order,2))),[size(Legendres,1) 1]).*Legendres(:,1:size(order,2)),2);

            WavefrontReconstructed = reshape(WavefrontReconstructed_Vect, ...
                [size(polynomials, 1) size(polynomials, 2)]);

            compensatedHologram = (abs(field_compensate).*exp(1i.* angle(field_compensate))./ ...
                exp((1i).*WavefrontReconstructed));

            if ~NoPistonCompensation
                PistonCorrtection = ones(size(compensatedHologram)).*Legendre_Coefficients(1);
                vortexPistonCorrection = compensatedHologram.*(exp(1i.*PistonCorrtection));
            end
            gplus2 = vortexPistonCorrection;
        end

        function [fx_max,fy_max] = maxVortex(holo)
            % Inputs:
            %   holo   -> 2-D real or complex array (Ny x Nx). Spatial-domain hologram.
            %             The function computes ft_holo = fftshift(fft2(fftshift(holo))).
            % Outputs:
            %   fx_max -> Scalar column index (x, 1-based) of the maximum |FT(holo)|.
            %   fy_max -> Scalar row    index (y, 1-based) of the maximum |FT(holo)|.
            
            ft_holo = fftshift(fft2(fftshift(holo)));
            [~, linearIndex] = max(abs(ft_holo), [], 'all', 'linear');
            [fy_max, fx_max] = ind2sub(size(ft_holo), linearIndex);
        end

        function [ref_wave] = reference_wave(M,N,m,n,lambda,dxy,fx_max,fy_max,k,fx_0,fy_0)
            % Inputs:
            %   M, N     -> Scalars. Image width (cols) and height (rows), respectively.
            %   m, n     -> N-by-M coordinate grids (rows=y and cols=x), typically centered, e.g.:
            %                 [n, m] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
            %   lambda   -> Scalar. Illumination wavelength (same length units as dxy).
            %   dxy      -> Scalar. Pixel pitch/size in the image/object plane (same units as lambda).
            %   fx_max   -> Scalar. Column index (x, 1-based) of the measured +1 order peak.
            %   fy_max   -> Scalar. Row    index (y, 1-based) of the measured +1 order peak.
            %   k        -> Scalar. Wavenumber, k = 2*pi/lambda.
            %   fx_0     -> Scalar. Column index of the FFT center (usually M/2 after fftshift).
            %   fy_0     -> Scalar. Row    index of the FFT center (usually N/2 after fftshift).
            %
            % Output:
            %   ref_wave -> N-by-M complex array. Plane-wave reference.

            theta_x = asin((fx_0 - fx_max) * lambda / (M * dxy));
            theta_y = asin((fy_0 - fy_max) * lambda / (N * dxy));
            ref_wave = exp(1i * k * (sin(theta_x) * n * dxy + sin(theta_y) * m * dxy));
        end

        function [BW_shifted] = shiftMaskByMax(I1, I2, BW) 
            % Inputs:
            %   I1 -> 2-D real or complex array. Source image (if complex, abs(.) is used).
            %   I2 -> 2-D real or complex array. Target image (if complex, abs(.) is used).
            %   BW -> 2-D logical (or numeric 0/1) mask of the same size as I1/I2. 
            % Output:
            %   BW_shifted -> 2-D logical mask, same size as BW. Result of shifting BW by
            %                 (dx, dy), where (dx, dy) is the column/row displacement that
            %                 maps the global maximum of I1 onto the global maximum of I2.
           
            % Make sure we work with real images
            if ~isreal(I1), I1 = abs(I1); end
            if ~isreal(I2), I2 = abs(I2); end

            % 1) Locate maxima
            [~, k1] = max(I1(:)); [y1, x1] = ind2sub(size(I1), k1);
            [~, k2] = max(I2(:)); [y2, x2] = ind2sub(size(I2), k2);

            % 2) Pixel shift needed to move I1 → I2
            dx = x2 - x1;                % columns
            dy = y2 - y1;                % rows
            
            % 3) Apply same integer shift to the binary mask (no wrap-around)
            [M,N] = size(BW);
            BW_shifted = false(M,N);

            [r,c] = find(BW);            % coordinates of 1's
            r2 = r + dy;  c2 = c + dx;   % shifted coordinates

            % keep only points that stay inside the image
            keep = r2>=1 & r2<=M & c2>=1 & c2<=N;
            idx = sub2ind([M,N], r2(keep), c2(keep));
            BW_shifted(idx) = true;
        end

        function Ft_norm = normalization_FT(Given_FT)
            % Input:
            %   Given_FT -> 2-D numeric array. Complex.
            % Output:
            %   Ft_norm  -> 2-D numeric array (same size as Given_FT) scaled to [0, 1].

            Ft_Gminus = Given_FT;
            Ft_norm = (Ft_Gminus - min(min(Ft_Gminus)))/(max(max(Ft_Gminus)) - min(min(Ft_Gminus)));
        end

        function rad_ratio = CalcEnhandcement (CircSum,CricMaskCent1)  
            % Inputs:
            %   CircSum        -> 2-D numeric array. Image/map whose positive values define
            %                     the "flower" region. Internally binarized as F = (CircSum > 0).
            %   CricMaskCent1  -> 2-D numeric/logical array. Binary mask of a single circle
            %                     centered at the same frame. Internally binarized as C = (CricMaskCent1 > 0).
            % Output:
            %   rad_ratio -> Scalar. Enhancement ratio defined as: r_f / r_c

            F = CircSum>0;                    
            C = CricMaskCent1>0;                  
            % --- area (sub-pixel estimate) ---
            Ac = bwarea(C);
            % --- single-circle radius from its area ---
            r_c = sqrt(Ac/pi);
            % --- circumscribed radius of the flower ---
            ctr = regionprops(F,'Centroid');  ctr = ctr.Centroid;   % [x y]
            BWb = bwperim(F);                 % boundary pixels of the flower
            [y,x] = find(BWb);
            r_f = max(hypot(x-ctr(1), y-ctr(2)));   % farthest boundary point
            % --- metrics ---
            rad_ratio  = r_f / r_c;
        end

    end
end
