function Cell_slice = pl_cell_slice(Cell,ndx)
%
% Slices cell elements keeping the indices 'ndx'

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.

N = length(Cell);

%compute mean
for i = 1:N
    Cell_slice{i} = Cell{i}(ndx);
end


