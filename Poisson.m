classdef Poisson < handle
    
    methods(Static)
        function [I_y I_x] = get_gradients(I)
            % Get the vertical (derivative-row) and horizontal (derivative-column) gradients
            % of an image. (or multiple images)
            
            I_y = diff(I,1,1);
            I_y = padarray(I_y, [1 0], 0, 'pre');
            I_x = diff(I,1,2);
            I_x = padarray(I_x, [0 1], 0, 'pre');
        end
        
        function [I] = solve(t_y, t_x, mask, t_y_weights, t_x_weights)
            % Solve for the image which best matches the target vertical gradients
            % t_y and horizontal gradients t_x, e.g. the one which minimizes sum of squares
            % of the residual
            % 
            %    sum of (I[i,j] - I[i-1,j] - t_y[i,j])^2 + (I[i,j] - I[i,j-1] - t_x[i,j])^2
            % 
            % Only considers the target gradients lying entirely within the mask.
            % The first row of t_y and the first column of t_x are ignored. 
            % Optionally, you may pass in an array with the weights corresponding to each target
            % gradient. The solution is unique up to a constant added to each of the pixels.
            
            if nargin < 5
                t_x_weights = [];
            end
            if nargin < 4
                t_y_weights = [];
            end
        
            if isempty(t_y_weights)
                t_y_weights = ones(size(t_y));
            end
            if isempty(t_x_weights)
                t_x_weights = ones(size(t_x));
            end
            [M, N] = size(mask);
            numbers = get_numbers(mask);
            A = get_A(mask, t_y_weights, t_x_weights);
            b = get_b(t_y, t_x, mask, t_y_weights, t_x_weights);
            
            % The solution is unique up to a constant added to each of the pixels.
            % We can reset the range of x to start from 0.
            warning off;
            x = A\b;
            x = x - min(x);
            warning on;
            
            % test x within a precision
            if max(abs(A*x-b)) > 1e-10
                fprintf('x = A\\b is not accurate.\n');
            end
            
            I = zeros(size(mask));
            for i = 1:M
                for j = 1:N
                    if mask(i,j)
                        I(i,j) = x(numbers(i,j));
                    end
                end
            end
        end
        
        function [I] = solve_L1(t_y, t_x, mask)
            % Same as solve(), except using an L1 penalty rather than least squares.
            
            EPSILON = 0.0001;
            
            % We minimize the L1-norm of the residual
            %
            %    sum of |r_i|
            %    r = Ax - b
            %
            % by alternately minimizing the variational upper bound
            %
            %    |r_i| <= a_i * r_i^2 + 1 / (4 * a_i)
            %
            % with respect to x and a. When r is fixed, this bound is tight for a = 1 / (2 * r).
            % When a is fixed, we optimize for x by solving a weighted
            % least-squares problem.
            
            I = Poisson.solve(t_y, t_x, mask);
            for i = 1:20
                [I_y I_x] = Poisson.get_gradients(I);
                t_y_err = mask .* abs(I_y - t_y);
                t_x_err = mask .* abs(I_x - t_x);
                
                t_y_weights = 1. ./ (2. * max(t_y_err, EPSILON));
                t_x_weights = 1. ./ (2. * max(t_x_err, EPSILON));
                
                try
                    I = Poisson.solve(t_y, t_x, mask, t_y_weights, t_x_weights);
                catch ME
                    % Occasionally the solver fails when the weights get very large
                    % or small. In that case, we just return the previous iteration's
                    % estimate, which is hopefully close enough.
                    return;
                end
            end
        end
        
    end % method(Static)
    
end

function [numbers] = get_numbers(mask)
    [M N] = size(mask);
    numbers = zeros(size(mask));
    count = 0;
    for i = 1:M
        for j = 1:N
            if mask(i,j)
                count = count + 1;
                numbers(i,j) = count;
            end
        end
    end
end

function [b] = get_b(t_y, t_x, mask, t_y_weights, t_x_weights)
    [M N] = size(mask);
    t_y = t_y(2:end,:);
    t_y_weights = t_y_weights(2:end,:);
    t_x = t_x(:,2:end);
    t_x_weights = t_x_weights(:,2:end);
    numbers = get_numbers(mask);
    K = max(numbers(:));
    b = zeros(K,1);

    % horizontal derivatives
    for i = 1:M
        for j = 1:N-1
            if mask(i,j) && mask(i,j+1)
                n1 = numbers(i,j);
                n2 = numbers(i,j+1);

                % row (i,j): -x_{i,j+1} + x_{i,j} + t
                b(n1) = b(n1) - t_x(i,j) * t_x_weights(i,j);

                % row (i,j+1): x_{i,j+1} - x_{i,j} - t
                b(n2) = b(n2) + t_x(i,j) * t_x_weights(i,j);
            end
        end
    end

    % vertical derivatives
    for i = 1:M-1
        for j = 1:N
            if mask(i,j) && mask(i+1,j)
                n1 = numbers(i,j);
                n2 = numbers(i+1,j);

                % row (i,j): -x_{i+1,j} + x_{i,j} + t
                b(n1) = b(n1) - t_y(i,j) * t_y_weights(i,j);

                % row (i,j+1): x_{i+1,j} - x_{i,j} - t
                b(n2) = b(n2) + t_y(i,j) * t_y_weights(i,j);
            end
        end
    end
end
        
function [A] = get_A(mask, t_y_weights, t_x_weights)
    [M N] = size(mask);
    numbers = get_numbers(mask);
    K = max(numbers(:));

    t_y_weights = t_y_weights(2:end,:);
    t_x_weights = t_x_weights(:,2:end);

    % herizontal derivatives
    count = 0;
    for i = 1:M
        for j = 1:N-1
            if mask(i,j) && mask(i,j+1)
                count = count + 1;
            end
        end
    end
    data = zeros(4*count,1);
    row  = zeros(4*count,1);
    col  = zeros(4*count,1);
    count = 0;
    for i = 1:M
        for j = 1:N-1
            if mask(i,j) && mask(i,j+1)
                n1 = numbers(i,j);
                n2 = numbers(i,j+1);

                % row (i,j): -x_{i,j+1} + x_{i,j} + t
                row(4*count+1) = n1;
                col(4*count+1) = n2;
                data(4*count+1) = -t_x_weights(i,j);

                row(4*count+2) = n1;
                col(4*count+2) = n1;
                data(4*count+2) = t_x_weights(i,j);

                % row (i, j+1): x_{i,j+1} - x_{i,j} - t
                row(4*count+3) = n2;
                col(4*count+3) = n2;
                data(4*count+3) = t_x_weights(i,j);

                row(4*count+4) = n2;
                col(4*count+4) = n1;
                data(4*count+4) = -t_x_weights(i,j);
                
                count = count + 1;
            end
        end
    end

    data1 = data;
    row1  = row;
    col1  = col;

    % vertical derivatives
    count = 0;
    for i = 1:M-1
        for j = 1:N
            if mask(i,j) && mask(i+1,j)
                count = count + 1;
            end
        end
    end
    data = zeros(4*count,1);
    row  = zeros(4*count,1);
    col  = zeros(4*count,1);
    count = 0;
    for i = 1:M-1
        for j = 1:N
            if mask(i,j) && mask(i+1,j)
                n1 = numbers(i,j);
                n2 = numbers(i+1,j);

                % row (i,j): -x_{i+1,j} + x_{i,j} + t
                row(4*count+1) = n1;
                col(4*count+1) = n2;
                data(4*count+1) = -t_y_weights(i,j);

                row(4*count+2) = n1;
                col(4*count+2) = n1;
                data(4*count+2) = t_y_weights(i,j);

                % row (i, j+1): x_{i+1,j} - x_{i,j} - t
                row(4*count+3) = n2;
                col(4*count+3) = n2;
                data(4*count+3) = t_y_weights(i,j);

                row(4*count+4) = n2;
                col(4*count+4) = n1;
                data(4*count+4) = -t_y_weights(i,j);
                
                count = count + 1;
            end
        end
    end

    data2 = data;
    row2  = row;
    col2  = col;

    data = [data1; data2];
    row  = [row1; row2];
    col  = [col1; col2];

    A = sparse(row, col, data, K, K);
end
