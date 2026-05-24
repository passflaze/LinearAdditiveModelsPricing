function [ccdf] = ccdf_AB_FFT(eta, k, T1, T2, sigma_T1, sigma_T2)
    % Conditional complementary CDF (survival function) of the AB log-price
    % increment between reset dates T1 and T2, computed via Lewis-FFT.
    %
    % Returns 1 - F_{T2|T1}(x) on the log-moneyness grid zk.
    % The CDF is obtained as F(x) = 1 - ccdf(x).
    %
    % sigma_T1, sigma_T2 = entries of sigma_t from calibrateAB, i.e. sigmaATM/I_0 (Eq.15).
    % Calling convention:
    %   [~, ~, sigma_t, ~] = calibrateAB(...);
    %   ccdf = ccdf_AB_FFT(eta, kAB, yf(iT1), yf(iT2), sigma_t(iT1), sigma_t(iT2));
    %
    % Reference: Baviera-Massaria (2026).
    % Kernel: Digital  D(x)/B0 = (e^{ax}/2pi) * int phi(xi+ia)/(i xi - a) e^{-i xi x} dxi


    phi_T1 = charateristic_function_AB(T1, k, eta, sigma_T1);
    phi_T2 = charateristic_function_AB(T2, k, eta, sigma_T2);
    phi_cond = @(u) phi_T2(u) ./ phi_T1(u);    % conditional CF T2 | T1

    % FFT grid 
    % Digital kernel decays as 1/|xi|, slower than the Call (1/|xi|^2),
    % so we use a finer dz / larger N than call_AB_FFT.

    M  = 15;
    dz = 0.0025;
    N  = 2^M;
    dx = 2*pi / (N*dz);

    z1 = -dz * (N-1) / 2;
    x1 = -dx * (N-1) / 2;
    j  = 0:N-1;
    zk = (z1 + dz*j).';          
    xk = x1 + dx*j;              % integration grid

    % Lewis contour 
    % Strip f_t = (-p+ * I_0, p- * I_0), with p+- = +- eta + sqrt(eta^2 + 1/k).
    % For the Digital, the contour MUST satisfy a in (-p+, 0) strictly:
    % a = 0 is the pole 1/(i xi - a).
    I_0    = I0(0, k, eta);
    p_plus = eta + sqrt(eta^2 + 1/k);
    a      = max(-0.5, -0.45 * p_plus * I_0);

    % Lewis integrand for the Digital 
    int_kernel = @(csi) phi_cond(csi + 1i*a) ./ (1i*csi - a);
    fk = arrayfun(int_kernel, xk) .* exp(-1i * z1 * dx .* j);

    f_hat = exp(a*zk.') .* dx .* exp(-1i * x1 * zk.') .* fft(fk);

    ccdf = real(f_hat(:) / (2*pi));      % survival function = 1 - F(x)

    %cdf = 1 - ccdf;                      % CDF F(x) = 1 - ccdf(x)

    % Sanity check on boundary behaviour 
    % Expect ccdf -> 1 as x -> -inf and ccdf -> 0 as x -> +inf.
    % If the tails are off by more than a few percent, increase M / refine dz / pull a more negative.
    tol = 5e-2;
    if abs(ccdf(1) - 1) > tol || abs(ccdf(end)) > tol
        warning('ccdf_AB_FFT:tails', ...
            'CCDF boundary off: ccdf(-inf)=%.4f (want 1), ccdf(+inf)=%.4f (want 0). Check grid / contour a.', ...
            ccdf(1), ccdf(end));
    end
end
