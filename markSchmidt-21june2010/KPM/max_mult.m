function M = max_mult(A, B)
%% Just like the matrix multipliation of A and B, except we max instead of sum
M = squeeze(max(bsxfun(@times, A, permute(B, [3 1 2])), [], 2)); 
end