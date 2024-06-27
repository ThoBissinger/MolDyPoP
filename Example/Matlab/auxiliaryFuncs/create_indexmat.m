function indexmat = create_indexmat(TCF_length,bin_length)
%CREATE_INDEXMAT Creates index matrix for given TCF_length and bin_length
    indexmat = zeros(TCF_length,bin_length);
    for i = 0:TCF_length-1
        indexmat(i + 1,:) = 1 + i * bin_length : (i + 1) * bin_length;
    end
end