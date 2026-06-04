function N_sim = size_Nsim_MC(sigma_fcn, cache_key, sig_inputs, opts)
%SIZE_NSIM_MC  Monte-Carlo path count for a target accuracy, cached on disk.
%
%   N_sim = SIZE_NSIM_MC(SIGMA_FCN, CACHE_KEY, SIG_INPUTS) returns the number
%   of MC paths needed so that the 95% confidence-interval half-width of the
%   price stays below the market-maker tolerance (Ex.4 hint: a market maker
%   quotes with a 10 bps bid-ask spread, so the MC noise must be smaller):
%
%       1.96 * sigma / sqrt(N) <= target   =>   N = (1.96*sigma/target)^2
%
%   This matches exactly the CI that normfit returns for the price (muCI), so
%   sizing and reporting use the same quantity.
%
%   DECOUPLING + CACHING.  The std-dev sigma of the (discounted) payoff is
%   estimated by a cheap PILOT, run as SIGMA_FCN(N_pilot). To avoid repeating
%   that pilot (which also rebuilds the FFT CDF) on every pricing call, the
%   resulting N is cached ON DISK, keyed by CACHE_KEY together with SIG_INPUTS,
%   a numeric vector fingerprinting everything the pilot depends on (model
%   params, scale factors, strikes, discount factors, grid). If any of these
%   change (e.g. after recalibration) the key no longer matches and the pilot
%   is rerun automatically. So the heavy pilot runs ONCE per exotic /
%   forward-start option and is reused on subsequent prices.
%
% INPUTS
%   sigma_fcn  : handle @(N_pilot) -> sigma_est, the sample std of the
%                discounted payoff (may be a 1xN_K row across strikes).
%   cache_key  : char id of the option, e.g. 'CoC_GL', 'fwdstart_AB'.
%   sig_inputs : numeric vector fingerprinting the pilot inputs.
%   opts       : (optional) struct, fields (defaults in brackets)
%                  .target_error [10e-4 = 10 bps]
%                  .N_pilot      [1e4]   pilot path count
%                  .cap          [1e7]   hard ceiling on N_sim
%                  .recompute    [false] ignore cache, force a fresh pilot
%                  .verbose      [true]
%                  .cache_file   [<thisdir>/Nsim_cache.mat]
%
% OUTPUT
%   N_sim      : required number of MC paths, in [1, cap].
%
% See also NTH_OUT, NORMFIT.

    if nargin < 4 || isempty(opts), opts = struct(); end
    def = struct('spread_bps', 10, 'ref', 1, 'target_error', [], ...
                 'N_pilot', 1e4, 'cap', 5e7, 'recompute', false, 'verbose', true, ...
                 'cache_file', fullfile(fileparts(mfilename('fullpath')), 'Nsim_cache.mat'));
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end
    end

    % Market-maker tolerance (Ex.4). Two readings of "10 bps":
    %   * absolute on the option price  -> ref = 1      -> target = 10e-4
    %   * 10 bps of the underlying       -> ref = F      -> target = 10e-4 * F
    % Pass opts.ref = forward for the latter, or opts.target_error directly to
    % override (e.g. a vol-based spread).
    if isempty(opts.target_error)
        opts.target_error = opts.spread_bps * 1e-4 * opts.ref;
    end

    % The target enters N, so fold it into the fingerprint: changing the spread
    % convention (ref/bps) invalidates a stale cached N automatically.
    sig_inputs = [sig_inputs(:); opts.target_error];

    % --- cache lookup --------------------------------------------------------
    cache = struct('key', {}, 'inputs', {}, 'N', {});
    if isfile(opts.cache_file)
        S = load(opts.cache_file, 'cache');
        if isfield(S, 'cache'), cache = S.cache; end
    end
    hit = 0;
    if ~opts.recompute
        for i = 1:numel(cache)
            if strcmp(cache(i).key, cache_key) ...
                    && numel(cache(i).inputs) == numel(sig_inputs) ...
                    && isequaln(round(cache(i).inputs, 10), round(sig_inputs, 10))
                hit = i; break;
            end
        end
    end
    if hit > 0
        N_sim = cache(hit).N;
        if opts.verbose
            fprintf('[size_Nsim_MC] %-14s cache hit  -> N_sim = %d\n', cache_key, N_sim);
        end
        return;
    end

    % --- pilot: estimate sigma and size N ------------------------------------
    % Save/restore the RNG state around the pilot so that the production pricing
    % downstream sees the SAME random stream whether the pilot ran (cache miss)
    % or was skipped (cache hit). This keeps prices reproducible across runs.
    rng_state = rng;
    sigma_est = sigma_fcn(opts.N_pilot);
    rng(rng_state);
    sigma_max = max(sigma_est(:));        % size to the worst strike
    N_raw     = ceil((1.96 * sigma_max / opts.target_error)^2);
    N_sim     = min(max(N_raw, 1), opts.cap);
    capped    = N_raw > opts.cap;

    if opts.verbose
        fprintf('--- PILOT (%s) ---\n', cache_key);
        fprintf('  Pilot paths:   %d\n', opts.N_pilot);
        fprintf('  Std Dev (max): %.6f\n', sigma_max);
        if opts.ref == 1
            fprintf('  Target error:  %.6f  (%.0f bps, absolute on price)\n', ...
                opts.target_error, opts.spread_bps);
        else
            fprintf('  Target error:  %.6f  (%.0f bps x ref %.4f)\n', ...
                opts.target_error, opts.spread_bps, opts.ref);
        end
        if capped
            fprintf('  Required N:    %d  -> CAPPED to %d\n', N_raw, N_sim);
        else
            fprintf('  Required N:    %d\n', N_sim);
        end
        fprintf('------------------------\n');
    end
    if capped
        warning('size_Nsim_MC:capBinding', ...
            ['Required N_sim = %d exceeds cap %d for "%s"; achieved 95%% CI ' ...
             'half-width will be ~%.6f (> %.6f target).'], ...
            N_raw, opts.cap, cache_key, 1.96*sigma_max/sqrt(opts.cap), ...
            opts.target_error);
    end

    % --- store ---------------------------------------------------------------
    rec = struct('key', cache_key, 'inputs', sig_inputs, 'N', N_sim);
    if isempty(cache), cache = rec; else, cache(end+1) = rec; end %#ok<AGROW>
    try
        save(opts.cache_file, 'cache');
    catch ME
        warning('size_Nsim_MC:cacheSaveFailed', 'Could not write cache: %s', ME.message);
    end
end
