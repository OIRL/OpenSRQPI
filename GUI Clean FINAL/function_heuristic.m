%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: functions_evaluation                                                  %
%                                                                              %
% The script contains all implemented function for SI_DHM_heuristic            %
%                                                                              %
% Authors: Raul Castaneda, Sofia Obando, Carlos Trujillo, Rene Restrepo,       %
%           Ana Doblas.                                                        %
% Applied Optics Group EAFIT univeristy                                        %
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef function_heuristic
    methods(Static)



        function [cf] = costFunction_SIDHM(theta, FTHolo, fx_max,fy_max)

            % Function to apply The cost function to the first demodulation
            % Inputs: theta - angle for the demodulation
            %         FTHolo - Fourier Transform of the Hologram
            %         fx_max,fy_max - coordinates of the maximum within the
            %         +1 diffraction order
            % Output: cf - cost function for the first demodulation

            % cf = 0;
      
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,1);
            cf =  abs(Dplus(fx_max(1),fy_max(1))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [cf] = costFunction_SIDHMII(theta, FTHolo, fx_max,fy_max)

            % Function to apply The cost function to the second demodulation
            % Inputs: theta - angle for the demodulation
            %         FTHolo - Fourier Transform of the Hologram
            %         fx_max,fy_max - coordinates of the maximum within the
            %         +1 diffraction order
            % Output: cf - cost function for the second demodulation

            % cf = 0;
        
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,2);
            cf =  abs(Dplus(fx_max(2),fy_max(2))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [D] = demComp2SIDHM(theta, H)

            % Function to demodulate de periodic pattern of SI-DHM
            % Inputs: theta - pair of angles for demodulate
            %         H - Matrix with the two holograms, dimensions: (N,M,2)
            % Output: D - Metrix with the two demdulation G+ and G- dimensions: (N,M,2)

            [X, Y, no] = size(H);
            D = zeros(X,Y,no);
            M = 1/2*[exp(1i*0) exp(-1i*0);exp(1i*theta) exp(-1i*theta)];
            Minv = pinv(M);

            D(:,:,1) = Minv(1,1).*H(:,:,1) + Minv(1,2).*H(:,:,2);
            D(:,:,2) = Minv(2,1).*H(:,:,1) + Minv(2,2).*H(:,:,2);
        end


        function  [ref] = phase_rec(filename, dx, dy, lambda)
            % Function to compensate the phase image
            % Inputs: Hologram - Hologram to compensate
            %         dx, dy  - pixel size in um
            %         lambda - wavelenght in um
            %   
            % Output: ref - numerical reference wave to compensate the hologram

            % Main function to perform phase reconstruction
            holo = filename;

            % Get the size of the hologram
            [N,M] = size(holo);
            % Create a meshgrid for the hologram
            [m,n] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);

            % Calculate the Fourier Transform of the hologram and shift the zero-frequency component to the center
            ft_holo = fftshift(fft2(fftshift(holo)));
            %imagesc(abs(ft_holo).^.2), title 'entrada del phase rec after demod'

            [~,idx] = max(ft_holo(:));
            % Get the maximum values of fx and fy
            [fy_max,fx_max] = ind2sub([N,M],idx);
            
            % Find the maximum value in the filtered spectrum
            % Define wavenumber
            k = 2 * pi / lambda;
            % Calculate the center frequencies for fx and fy
            fx_0 = M/2;
            fy_0 = N/2;

            % Calculate the angles for the compensation wave
            theta_x = asin((fx_0 - fx_max +1) * lambda / (M * dx));
            theta_y = asin((fy_0 - fy_max +1) * lambda / (N * dy));
            % Calculate the reference wave
            ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));
        end

    end
end
