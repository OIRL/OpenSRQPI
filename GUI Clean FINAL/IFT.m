function F = IFT(A)
if length(size(A)) == 3
    F = fftshift( ifftn( fftshift( A ) ) );
else
    F = fftshift( ifft2( fftshift( A ) ) );
end