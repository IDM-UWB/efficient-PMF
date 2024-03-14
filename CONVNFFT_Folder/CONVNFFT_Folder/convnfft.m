function A = convnfft(A, B,Npa)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
%
%   C = CONVNFFT(A, B)
%   OPTIONS is structure with following optional fields
%       - 'GPU', boolean. If GPU is TRUE Jacket/GPU FFT engine will be used
%       By default GPU is FALSE.
%       - 'Power2Flag', boolean. If it is TRUE, use FFT with length rounded
%       to the next power-two. It is faster but requires more memory.
%       Default value is TRUE.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
%
% See also conv, conv2, convn.
%
%   Author: Bruno Luong <brunoluong@yahoo.com>

nd = ndims(A);%max(ndims(A),ndims(B));
% work on all dimensions by default
dims = 1:nd;
dims = reshape(dims, 1, []); % row (needed for for-loop index)
% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
ifun = @(m,n) ceil((n-1)/2)+(1:m);
% faster FFT if the dimension is power of 2
lfftfun = @(l) 2^nextpow2(l);
% Do the FFT
subs(1:ndims(A)) = {':'};
for dim=dims
    m = Npa;%size(A,dim);
    n = Npa;%size(B,dim);
    % compute the FFT length
    l = lfftfun(m+n-1);
    A = fft(A,l,dim);
%     B = fft(B,l,dim);
    subs{dim} = ifun(m,n);
end
% inplace product to save 1/3 of the memory
inplaceprod(A,B);
for dim=dims
    A = ifft(A,[],dim);
end
% Make sure the result is real
A = real(A(subs{:}));

end % convnfft

