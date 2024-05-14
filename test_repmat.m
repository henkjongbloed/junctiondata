% Example matrix
N = 3; M = 4;
val = randn(N, M);  % Replace this with your actual matrix
p=5;
% Example weights vector
weights = 1:p;  % Replace this with your actual weights vector

% Number of copies in the third dimension
num_copies = numel(weights);

% Replicate the matrix along the third dimension
replicated_array = repmat(val, [1, 1, num_copies]);

% Multiply each slice by the corresponding weight
weighted_array = bsxfun(@times, replicated_array, permute(weights, [1, 3, 2]));

% Display the size of the weighted array
disp(size(weighted_array));