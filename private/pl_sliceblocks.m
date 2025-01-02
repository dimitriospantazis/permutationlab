function [sndx,ssize,nslices] = pl_sliceblocks(n,slicesize)
%
% Slices a dimension of size n into slices of size 'slicesize'

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.

nslices = ceil(n/slicesize);

for s = 1:nslices-1
    sndx{s} = (s-1)*slicesize+1 : slicesize*s;
    ssize(s) = slicesize;
end
sndx{nslices}=(nslices-1)*slicesize+1:n;
ssize(nslices) = length(sndx{nslices});