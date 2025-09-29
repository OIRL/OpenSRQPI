%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  FT
%
% Purpose:
%   Compute the centered Fourier transform of the input array:
%     - 2-D input  -> uses fft2 with fftshift pre/post
%     - 3-D input  -> uses fftn with fftshift pre/post (applies to all dims)
%
% Input:
%   A -> Numeric array (real or complex), either:
%        * 2-D (Ny x Nx): a single image/map
%        * 3-D (Ny x Nx x Nz): a volume or multi-channel stack
%
% Output:
%   F -> Same size as A, complex array containing the Fourier transform of A,
%        centered via fftshift (zero frequency at the middle).
%
% Notes:
%   - For 2-D: F = fftshift( fft2( fftshift(A) ) )
%     For 3-D: F = fftshift( fftn( fftshift(A) ) )
%   - Input can be real or complex; output is generally complex.
%   - Ensure A is properly windowed/padded beforehand if needed to control
%     spectral leakage or achieve desired frequency sampling.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F = FT(A)
if eq(length(size(A)),3)
    F = fftshift( fftn( fftshift( A ) ) );
else
    F = fftshift( fft2( fftshift( A ) ) );
end