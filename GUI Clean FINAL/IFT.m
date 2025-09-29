%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  IFT
%
% Purpose:
%   Compute the centered inverse Fourier transform of the input array:
%     - 2-D input  -> ifft2 with fftshift pre/post
%     - 3-D input  -> ifftn with fftshift pre/post (applies to all dims)
%
% Input:
%   A -> Numeric array (real or complex), either:
%        * 2-D (Ny x Nx): frequency-domain image/map
%        * 3-D (Ny x Nx x Nz): frequency-domain volume or multi-channel stack
%
% Output:
%   F -> Complex array, same size as A. The inverse Fourier transform of A with
%        zero frequency centered (via fftshift before and after the transform).
%
% Notes:
%   - 2-D path: F = fftshift( ifft2( fftshift(A) ) )
%     3-D path: F = fftshift( ifftn( fftshift(A) ) )
%   - The result is generally complex; use abs(.), angle(.), or real(.) as needed.
%   - MATLAB's ifft2/iffftn include 1/(N) scaling internally; no extra scaling here.
%   - Ensure A is fftshift-centered to match the centering assumed by this routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F = IFT(A)
if length(size(A)) == 3
    F = fftshift( ifftn( fftshift( A ) ) );
else
    F = fftshift( ifft2( fftshift( A ) ) );
end