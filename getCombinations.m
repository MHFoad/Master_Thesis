function output = getCombinations(input)
N_input = length(input);
input_ind = 1:N_input;
input_sizes = cellfun(@numel,input);

output = cell(1,N_input);
for output_ind=1:N_input
    input_var = input{input_ind(output_ind)};
    dim_var = ones(1,N_input);
    dim_var(output_ind) = numel(input_var);
    input_var = reshape(input_var,dim_var);
    dim_var = input_sizes;
    dim_var(output_ind) = 1;
    output{output_ind} = repmat(input_var,dim_var);
end