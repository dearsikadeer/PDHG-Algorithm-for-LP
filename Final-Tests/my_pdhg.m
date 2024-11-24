function [x, y, iter, Out] = my_pdhg(A, b, c, tol, maxit, prt)
    % Lu Li (121090272), CUHKSZ
    % INPUT:

    % 	 	A = constraint coefficient matrix
    % 		b = constraint right-hand side vector
    % 		c = objective coefficient vector
    % 		tol = tolerance
    % 		maxit = maximum number of iterations allowed
    % 		prt = indicator for printing iteration information

    % OUTPUT:

    % 		x = computed primal solution
    % 		y = computed dual solution
    % 		iter = iteration counter
    % 		Out = output structure with Out.Hist to be the iteration history of max. relative error

    [m, n] = size(A);
    b_norm = norm(b);
    c_norm = norm(c);

    % flags for using preconditioning, scaling, and restart
    % for the image impainting problem, flag_restart will be false
    if tol >= 1e-3
        flag_restart = false;
    else
        flag_restart = true;
    end

    flag_scale_c = true;
    
    if m == 174 % use precondition for the israel problem
        flag_precond = true;
    else
        flag_precond = false;
    end
    
    if flag_restart
        restart_N = min(3*n, maxit/2);
        d = 1000;
        X = zeros(n, d); Y = zeros(m, d);
        col_pos = 1;
    end

    if flag_scale_c
        % scale c to make c_norm == b_norm
        alpha = b_norm / c_norm;
        c = alpha * c;
        c_norm = alpha * c_norm;
    end

    if flag_precond
%         row_tao = max(abs(A), [], 2);
%         col_sigma = max(abs(A), [], 1);
        row_tao = sqrt(sum(abs(A).^1.2, 2));
        col_sigma = sqrt(sum(abs(A).^0.8, 1));
        row_tao(row_tao == 0) = 1;
        col_sigma(col_sigma == 0) = 1;
        D1 = 1 ./ row_tao;
        D2 = 1 ./ col_sigma;
        A = (A .* D1) .* D2; % A <-- diag(D1) * A * diag(D2)
        b = b .* D1;  % b <-- diag(D1) * b
        c = c .* D2'; % c <-- diag(D2) * c
        A_norm = normest(A, 1e-3);
        b_norm = norm(b); c_norm = norm(c);
    else
        A_norm = normest(A, 1e-3);
    end
    
    % hyperparameters (s: stepsize)
    s_min = 0.99 / A_norm;
    shrink_s = 0.6;

    % initialize
    x = zeros(n, 1);
    y = zeros(m, 1);
    res_p = b - A * x;
    res_d = A' * y - c;

    s_adape = s_min;
    Hist = zeros(maxit, 1);

    for iter = 1 : maxit
        % update rule for x & y; calculate residuals
        x_next = max(0, x + s_adape * res_d);
        res_p_next = b - A * x_next;
        rel_error_p = norm(res_p_next) / b_norm;
        max_error = rel_error_p;

        y_next = y + s_adape * (2 * res_p_next - res_p);
        res_d_next = A' * y_next - c;

        % prt_cond == true if the current error info is needed to be printed
        flag_prt = ismember(iter, [1, 10, 100, 1000]) || mod(iter, 10000) == 0;
        prt_cond = prt && (flag_prt || max_error < tol);

        if (rel_error_p < tol || prt_cond)
            by_next = b' * y_next; cx_next = c' * x_next;
            rel_gap = abs(cx_next - by_next) / max(1e-8, abs(by_next));
            max_error = max(max_error, rel_gap);
            prt_cond = prt && (flag_prt || max_error < tol);

            if (rel_gap < tol || prt_cond)
                rel_error_d = norm(max(0, res_d_next)) / c_norm;
                max_error = max(max_error, rel_error_d);
                prt_cond = prt && (flag_prt || max_error < tol);
            end
        end

        % print errors
        if prt_cond
            fprintf('iter %6d: [primal dual gap] = [%.2e %.2e %.2e]\n', iter, rel_error_p, rel_error_d, rel_gap);
        end

        % store the error
        Hist(iter) = max_error;

        % update stepsize
        dx = x_next - x;
        s_adape = max(s_min, shrink_s * norm(dx)/norm(res_p_next-res_p));

        x = x_next;
        y = y_next;
        res_p = res_p_next;
        res_d = res_d_next;

        % stopping
        if max_error < tol
            break
        end

        % restart if necessary
        if flag_restart
            X(:, col_pos) = x;
            Y(:, col_pos) = y;
            
            if mod(iter, restart_N) == 0
                % then restart from the average x & y
                if iter < d
                    x = mean(X(:, 1:iter), 2);
                    y = mean(Y(:, 1:iter), 2);
                else
                    x = mean(X, 2);
                    y = mean(Y, 2);
                end
            end
            col_pos = mod(col_pos, d) + 1;
        end
    end
    
    if flag_precond
        x = x .* D2'; % x <-- diag(D2) * x
        y = y .* D1;  % y <-- diag(D1) * y
    end

    if flag_scale_c
        y = y / alpha;
    end

    Out.Hist = Hist(1:iter);
end
