
a = 0;
b = 0;
n_state = 2;

for i=1:1000
    epsilon = rand(n_state, n_state);

    rowsum = sum(epsilon,2);
    epsilon = bsxfun(@rdivide, epsilon, rowsum);

    if any(sum(epsilon,2) ~= ones(n_state,1))
        a = a + 1;
    end
    
    if any(sum(epsilon,2) - ones(n_state,1) ~= zeros(n_state,1))
        b = b+1;
    end
end

fprintf('a= %i \n', a)
fprintf('b= %i \n', b)