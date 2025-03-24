function W = build_W_from_H(H)
    % Determine the number of terms by counting the fields in the structure
    fields = fieldnames(H); % Get the field names of the structure
    num_terms = length(fields) / 2;  % Each term has (H_k, G_k)
    lambda = randn(num_terms, 1); % Random coefficients

    % Initialize W (size based on the dimensions of the Hermitian matrices)
    first_field = fields{1}; % Get the first field name
    matrix_size = size(H.(first_field)); % Size of the Hermitian matrix
    W = zeros(matrix_size(1) * matrix_size(1), matrix_size(2) * matrix_size(2)); % Initialize W

    % Loop through each pair of Hermitian matrices (H_k, G_k)
    for k = 1:num_terms
        H_k = H.(sprintf('H%d', 2*k - 1));  % Extract H_k from H1, H3, ...
        G_k = H.(sprintf('H%d', 2*k));      % Extract G_k from H2, H4, ...
        
        % Add the contribution to W
        W = W + lambda(k) * kron(H_k, G_k);
    end
end
