function [data_cell] = pl_mat2cell(data)
%
% Converts an N dimensional array (n x m1 x m2 x ...) into a cell (n x 1) 
% of N-1 dimensional arrays (m1x m2 x ...)
% Data in this form will be processed faster in permutationlab
%

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.


n = size(data,1);

for i = 1:n
    data_cell{i} = shiftdim(data(i,:,:),1);
end






