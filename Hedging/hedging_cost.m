function cost = hedging_cost(trade, F, costRule)
% HEDGING_COST  Bid-ask transaction cost of a (re)hedging trade.
%
%   Cost convention (agreed): the half-spread is charged on the NOTIONAL,
%   i.e. on the forward/underlying level F, for every unit traded.
%
%       cost = fut_bp * |nF| * F        (underlying future, 1 bp)
%            + opt_bp * (|nC| + |nP|) * F   (each PV option, 4 bp)
%
%   Pass the TRADE (change in positions) so that at rebalancing only the
%   incremental quantity is charged:  trade = positions_k - positions_{k-1}.
%   At inception pass the full initial positions.
%
% INPUTS:
%   trade    - (struct) .nC .nP .nF  signed quantities actually traded
%   F        - (scalar) underlying/forward level at the trade date
%   costRule - (struct) .fut_bp (default 1) .opt_bp (default 4)
%
% OUTPUT:
%   cost     - (scalar) total bid-ask cost in $ (>= 0).

    if nargin < 3 || isempty(costRule)
        costRule = struct('fut_bp', 1, 'opt_bp', 4);
    end

    bp = 1e-4;
    cost =  costRule.fut_bp * bp * abs(trade.nF)                  * F ...
          + costRule.opt_bp * bp * (abs(trade.nC) + abs(trade.nP)) * F;
end
