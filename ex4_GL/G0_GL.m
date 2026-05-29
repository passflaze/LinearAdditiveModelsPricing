function G0 = G0_GL(alpha, beta)
% G0_GL  Prezzo NORMALIZZATO di una call ATM-forward sotto il modello GL.
%
% Una call ATM-forward ha strike K = F(t0,T), quindi moneyness normalizzato
%   chi = (K - F)/(sigma_T*sqrt(T)) = 0.
% La funzione prezzo normalizzata e'  G(chi) = E[(z - chi)^+]  con z ~ GL
% standardizzata (CF = cf_GL). Valutata in chi = 0:
%
%       G0 = G(0) = E[z^+] = \int_0^inf x * pdf_GL(x) dx
%
% Il prezzo (scontato, dimensionale) di una call ATM-forward e':
%       C_ATM(T) = B(0,T) * sigma_T * sqrt(T) * G0
% dove sigma_T e' la SCALA NORMALIZZATA sigma_t = sigma_ATM/I_0
% (NON sigma_ATM: vedi nota di coerenza nel main).
%
% Inputs
%   alpha, beta : parametri di forma GL
% Output
%   G0 : E[z^+], numero puro (dipende solo dalla forma, non da T ne' dal forward)

    integrand = @(x) pdf_GL(alpha, beta, x) .* x;
    G0 = quadgk(integrand, 0, inf);
end