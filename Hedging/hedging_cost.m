function costs = hedging_cost(trade_w, instr_greeks, costRule)
% HEDGING_COST  Transaction costs (slippage) and hedging premium for a
%               generic basket of instruments.
%
%   Cost convention:
%   1. SLIPPAGE: half-spread charged on the traded PV — opt_bp on option
%      premiums, fut_bp on the futures notional — for every unit traded.
%   2. PREMIUM : actual cash flow to acquire the options (futures: none).
%
% INPUTS:
%   trade_w      - (N x 1) signed quantities actually traded.
%   instr_greeks - (N x 1 struct array) with .kind and .price per instrument.
%   costRule     - (struct) .fut_bp (default 1) .opt_bp (default 4)
%
% OUTPUT:
%   costs        - (struct) .slippage, .premium, .total in $

    if nargin < 3 || isempty(costRule)
        costRule = struct('fut_bp', 1, 'opt_bp', 4);
    end

    bp = 1e-4;
    N  = numel(instr_greeks);

    cost_slippage = 0;
    cost_premium  = 0;
    for j = 1:N
        px = instr_greeks(j).price;
        if strcmpi(instr_greeks(j).kind, 'future')
            % Futures: slippage on notional, no upfront premium.
            cost_slippage = cost_slippage + costRule.fut_bp * bp * abs(trade_w(j)) * px;
        else
            cost_slippage = cost_slippage + costRule.opt_bp * bp * abs(trade_w(j)) * px;
            cost_premium  = cost_premium  + trade_w(j) * px;
        end
    end

    costs = struct();
    costs.slippage = cost_slippage;
    costs.premium  = cost_premium;
    costs.total    = cost_slippage + cost_premium;
end
