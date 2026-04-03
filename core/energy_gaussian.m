function [energy] = energy_gaussian( params, cfg )
%ENERGY_GAUSSIAN Calculate parameter energy for gaussian likelihood model
% 
%  [energy] = energy_generic( params, cfg )

% check prior probability on parameters
logprior = cfg.logpdf_prior_fcn(params);
if isinf(logprior)
    energy = cfg.big_energy;
    return;
end
energy = -logprior;

% equillibrate, if required
if isfield(cfg, 'equilibrate_fcn')
    [err,state,~] = cfg.equilibrate_fcn( params );
    if (err)
        energy = cfg.big_energy;
        return;
    end
    init = state(end,:);
else
    init = [];
end

% simulate experiments
parworkers = 0;
if isfield(cfg, 'parallel') && cfg.parallel
    parworkers = Inf; % Use all available workers
    if isfield(cfg, 'maxlabs')
        parworkers = min([cfg.maxlabs, parworkers]);
    end
end

expt_energy = zeros(cfg.nexpt, 1);
expt_time = zeros(cfg.nexpt, 1);

parfor (d = 1:cfg.nexpt, parworkers)

    t_expt = tic;

    % simulate experiment (default initial conditions)
    [err,~,obsv] = cfg.data{d}.protocol_fcn( cfg.data{d}.time, init, params );
    if (err)
        expt_energy(d) = Inf;
        expt_time(d) = toc(t_expt);
        continue; % parfor does not support return
    end

    % normalize obsv
    [obsv] = norm_obsv( obsv, params, cfg );

    % heuristic penalties (do this before transformating obsv)
    penalty_val = 0;
    if isfield(cfg.data{d}, 'heuristic_penalty_fcn')
        penalty_val = cfg.data{d}.heuristic_penalty_fcn(obsv, params);
        if isinf(penalty_val)
            expt_energy(d) = Inf;
            expt_time(d) = toc(t_expt);
            continue;
        end
    end

    % if necessary, transform simulated trajectory for computing fitness
    if isfield(cfg, 'transform_sim_for_fit_fcn')
        obsv = cfg.transform_sim_for_fit_fcn(obsv,params);
    end

    % calculate log-likelihood as weighted sum of square errors
    loglike = nansum(nansum( -cfg.data{d}.weight .* (obsv - cfg.data{d}.mean).^2 ./ (2*(cfg.data{d}.stdev).^2) ));

    % accumulate experiment energy
    expt_energy(d) = penalty_val - loglike;

    expt_time(d) = toc(t_expt);

end

% check for any failures
if any(isinf(expt_energy))
    energy = cfg.big_energy;
    return;
end

% update energy with valid experiment results
energy = energy + sum(expt_energy);

% penalize for slow integrations
dt = sum(expt_time);
energy = energy + cfg.timepenalty*dt^2;

% all done
return;
