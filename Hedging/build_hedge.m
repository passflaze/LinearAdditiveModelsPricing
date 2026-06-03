function [w, residuals] = build_hedge(port_greeks, instr_greeks, target_greeks)
% BUILD_HEDGE  Optimal static-hedge weights for a generic basket of
%              instruments against a set of target Greeks.
%
%   Solves  A*w = -b  in the least-squares sense, where b are the portfolio
%   Greeks to neutralize and A(i,j) is the i-th target Greek of the j-th
%   instrument. With M = N (e.g. 3 vanillas vs Delta/Gamma/Vega) and a
%   well-conditioned A the hedge is exact.
%
% INPUTS:
%   port_greeks   - (struct) aggregate portfolio Greeks (.delta/.gamma/.vega).
%   instr_greeks  - (N x 1 struct array) Greeks of the available instruments,
%                   each with .kind .K .mat .price .delta .gamma .vega.
%   target_greeks - (cell) Greeks to hedge, e.g. {'Delta','Gamma','Vega'}.
%
% OUTPUTS:
%   w         - (N x 1) optimal signed quantities, one per instrument.
%   residuals - (M x 1) residual risk on the targeted Greeks.

    M = numel(target_greeks);
    N = numel(instr_greeks);

    % 1. Target vector (portfolio Greeks to neutralize)
    b = zeros(M, 1);
    for i = 1:M
        b(i) = port_greeks.(lower(target_greeks{i}));
    end

    % 2. Instruments matrix
    A = zeros(M, N);
    for j = 1:N
        for i = 1:M
            A(i, j) = instr_greeks(j).(lower(target_greeks{i}));
        end
    end

    % 3. Rank / conditioning checks
    if rank(A) < min(M, N)
        warning('BUILD_HEDGE:RankDeficient', ...
            ['Rank-deficient hedge matrix: instruments are linearly dependent ', ...
             'on the targeted Greeks (e.g. same-strike call/put share gamma & vega). ', ...
             'Use distinct strikes/maturities to span Delta-Gamma-Vega.']);
    end

    % Conditioning: a poorly conditioned A still solves exactly but produces
    % huge offsetting notionals and an unstable hedge (typical of wings too
    % deep OTM). Surface it so the basket strikes can be tightened.
    kappaA = cond(A);
    fprintf('Hedge Matrix Condition Number: %.2e\n', kappaA);
    if kappaA > 1e6
        warning('BUILD_HEDGE:IllConditioned', ...
            ['Ill-conditioned hedge matrix (cond = %.2e): expect large ', ...
             'offsetting positions and unstable rebalancing. Move strikes ', ...
             'closer to ATM (reduce the wing offset).'], kappaA);
    end

    % 4. Solve (exact or least squares)
    w = - (A \ b);
    residuals = b + A * w;

    % 5. Console output
    fprintf('\n--- Static Hedging Results ---\n');
    for j = 1:N
        if strcmpi(instr_greeks(j).kind, 'future')
            fprintf('  qty[%d] %-6s        T=%d : %+.4e\n', j, instr_greeks(j).kind, instr_greeks(j).mat, w(j));
        else
            fprintf('  qty[%d] %-6s K=%-8.4g T=%d : %+.4e\n', j, instr_greeks(j).kind, instr_greeks(j).K, instr_greeks(j).mat, w(j));
        end
    end
    fprintf('  Residual Risk -> ');
    for i = 1:M
        fprintf('%s: %+.2e', target_greeks{i}, residuals(i));
        if i < M, fprintf(' | '); end
    end
    fprintf('\n  cond(A) = %.2e\n', kappaA);
end
