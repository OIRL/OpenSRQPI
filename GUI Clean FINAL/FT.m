function F = FT(A)
if eq(length(size(A)),3)
    F = fftshift( fftn( fftshift( A ) ) );
else
    F = fftshift( fft2( fftshift( A ) ) );
end