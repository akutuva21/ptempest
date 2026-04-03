function [beta] = adapt_beta( beta, swap_acceptance, swap_idx, cfg )
% adapt beta
fprintf(1,'------------------------------------------------------\n');
fprintf(1,' Adapting chain temperatures . . .\n');
% based on acceptance rate of all chains
recent_swaps = swap_acceptance(:,(swap_idx-cfg.adapt_beta_interval+1):swap_idx);
acceptances = sum( recent_swaps == 1, 2 );
attempts = sum( recent_swaps ~= 0, 2 );
swap_acceptance_rate = acceptances ./ max(attempts, 1);
old_beta = beta;
% chain 1 has fixed temperature. adjust chains in order of coolest to hottest.
%   whenever a temperature changes, we proportionally increase the temperature of the hotter chains as well
for chain_idx = 2 : cfg.nchains
    % compute adaption factor
    adaption_factor = (swap_acceptance_rate(chain_idx-1)/cfg.optimal_swap_acceptance)^cfg.adapt_beta_rate;
    adaption_factor = min( [max([cfg.min_adaption_factor, adaption_factor]), cfg.max_adaption_factor] );
    % make sure beta doesn't get too close to its cooler neighbor
    adaption_factor = max( [adaption_factor, 2*beta(chain_idx)/(beta(chain_idx) + beta(chain_idx-1))] );
    % multiply this chain and higher temperature chains by adaption_factor
    for chain_idx2 = chain_idx : cfg.nchains
        beta(chain_idx2) = beta(chain_idx2)/adaption_factor;
    end
    fprintf(1,'  chain %d: recent acceptance=%-8.3g old temperature=%-6.4f new temperature=%-6.4f\n', ...
                chain_idx, swap_acceptance_rate(chain_idx-1), 1/old_beta(chain_idx), 1/beta(chain_idx) );
end

