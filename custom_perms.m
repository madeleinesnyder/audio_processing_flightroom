function P = custom_perms(vec)
    n = numel(vec);
    if n == 1
        P = vec;
    else
        P = [];
        % Recursively get permutations of the rest of elements
        for i = 1:n
            temp = vec;
            temp(i) = [];
            Ptemp = custom_perms(temp);
            
            % Append the removed element to the permutations of the n-1 subarray
            P = [P; vec(i) * ones(size(Ptemp, 1), 1), Ptemp];
        end
    end
end
