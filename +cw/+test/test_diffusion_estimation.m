function test_diffusion_estimation
    delta_t = 1;
    sigma = 1;
    
    steps = sigma*randn(2,40);
    
    delta_rs = sqrt( steps(1,:) .^ 2 + steps(2,:) .^2 );
    [data_cdf,x_locations] = ecdf(delta_rs);
    
    mu_max = max(steps,[],2);
    mu_min = min(steps,[],2);
    sigma_max = max(mu_max-mu_min);
    
    D_guess = std2(steps)^2/(2*delta_t)
    D_true = sigma^2/(2*delta_t)
    
    num_states = 1;
    num_iterations = 1;
    
    beta0 = [ones(1,num_states) / num_states, D_guess*ones(1,num_states)];

    problem = createOptimProblem('lsqcurvefit','x0',beta0,'objective',@(b,x) cdf_modelfun(b,x,delta_t),...
    'lb',zeros(1,num_states*2),'ub',[ones(1,num_states),sigma_max*ones(1,num_states)],'xdata',x_locations,'ydata',data_cdf);

    ms = MultiStart('PlotFcns',@gsplotbestf);
    [beta,errormulti] = run(ms,problem,num_iterations);

    num_D = numel(beta) / 2;

    weights = beta(1:num_D);
    weights = weights / sum(weights)
    Ds = beta(num_D+1:end)

    [Dsort, Is] = sort(Ds,'descend');

    residuals(num_states,:) = data_cdf-cdf_modelfun(beta,x_locations,delta_t);
    
end

function y = cdf_modelfun(b,x,delta_t)
    %%% Coefficients
    % weights followed by diffusion coefficients

    num_D = numel(b) / 2;
    
    weights = b(1:num_D);
    weights = weights / sum(weights);
    Ds = b(num_D+1:end);
    
    exponentials = zeros(numel(x),num_D);
    
    for D_idx = 1:num_D
        
        exponentials(:,D_idx) = weights(D_idx) .* exp(-x.^2 ./ (4*Ds(D_idx)*delta_t));
    end
    
    y = 1 - sum(exponentials,2);

end