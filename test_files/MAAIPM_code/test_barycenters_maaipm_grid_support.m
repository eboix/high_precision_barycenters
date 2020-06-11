function c = test_barycenters_maaipm_grid_support(n,k,eps_inv,fprefix)
    % n is support size of distribution
    % k is number of distributions
    % eps_inv is side length of discretization grid

    % Based on code from:
    % Ge, D., Wang, H., Xiong, Z., & Ye, Y. (2019).
    % Interior-point Methods Strike Back: Solving the Wasserstein Barycenter Problem.
    % arXiv preprint arXiv:1905.12895.

    fname = ['../experiment_data/', fprefix, '_n', num2str(n), 'k', num2str(k), '.txt'];
    fileID = fopen(fname,'r');
    A = fscanf(fileID, '%f');
    fclose(fileID);
    A = reshape(A, 3, n*k);

    prob_instance.stride = repmat(n, [1,k]);
    prob_instance.supp = A(1:2,:);
    prob_instance.w = A(3,:);


    [X, Y] = meshgrid(linspace(-1,1,eps_inv), linspace(-1,1,eps_inv));
    num_support_points = length(X(:));
    supports.supp = [X(:)'; Y(:)'];
    supports.w = repmat(1 / num_support_points, [1, num_support_points]);
    db = { prob_instance };
    c0 = { supports };


    %% Compute Wasserstein Barycenter (Pre-specified support MAAIPM)
    options.method='fixed_maaipm'; % {'ibp','gurobi','badmm', 'admm', 'ibp'}
    options.ipmouttolog = 1;
    options.ipmtol_primal_dual_gap = 1e-7;
    if num_support_points > 2*n %SLRM/DLRM
        options.largem = 1;
    else
        options.largem = 0;
    end

    [c, OT]=Wasserstein_Barycenter(db, c0, options);

    sum(c{:}.w)
end
