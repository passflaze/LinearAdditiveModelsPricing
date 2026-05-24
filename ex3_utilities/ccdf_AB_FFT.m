function cdf = ccdf_AB_FFT(eta, k, T1, T2, sigma_T1, sigma_T2, x)
    % Conditional complementary CDF (survival function) of the AB log-price
    % increment between reset dates T1 and T2, computed via Lewis-FFT.
    %
    % Returns 1 - F_{T2|T1}(x) 
    % The CDF is obtained as F(x) = 1 - ccdf(x).
    %
    % sigma_T1, sigma_T2 = entries of sigma_t from calibrateAB, i.e. sigmaATM/I_0 (Eq.15).
    % Calling convention:
    %   [~, ~, sigma_t, ~] = calibrateAB(...);
    %   ccdf = ccdf_AB_FFT(eta, kAB, yf(iT1), yf(iT2), sigma_t(iT1), sigma_t(iT2));
    %.  x is the vector of moneyness points where to evaluate the CCDF
    %
    % Reference: Baviera-Massaria (2026).
    % Kernel: Digital  D(x)/B0 = (e^{ax}/2pi) * int phi(xi+ia)/(i xi - a) e^{-i xi x} dxi


    phi_T1 = charateristic_function_AB(T1, k, eta, sigma_T1);
    phi_T2 = charateristic_function_AB(T2, k, eta, sigma_T2);
    phi_cond = @(u) phi_T2(u) ./ phi_T1(u);    % conditional CF T2 | T1 due to additivity

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
    zk = z1 + dz*j;              
    xk = x1 + dx*j;

    % Lewis contour 
    % Strip f_t = (-p+ * I_0, p- * I_0), with p+- = +- eta + sqrt(eta^2 + 1/k).
    % For the Digital, the contour MUST satisfy a in (-p+, 0) strictly:
    % a = 0 is the pole 1/(i xi - a).
    I_0    = I0(0, k, eta);
    p_plus = eta + sqrt(eta^2 + 1/k);
    a=-0.02;

    % Lewis integrand for the Digital.
    % At extreme |xi|, both phi_T2 and phi_T1 decay to machine zero -> 0/0 = NaN.
    % These tail contributions are numerically zero; replace NaN/Inf with 0.
    int_kernel = @(csi) phi_cond(csi + 1i*a) ./ (1i*csi - a);
    fk_raw = arrayfun(int_kernel, xk);
    fk_raw(~isfinite(fk_raw)) = 0;
    fk = fk_raw .* exp(-1i * z1 * dx .* j);

    f_hat =  dx .* exp(-1i * x1 * zk) .* fft(fk);

    f_hat = interp1(zk, f_hat, x, 'spline');

    ccdf = real(f_hat(:) / (2*pi)).*exp(a*x);      % survival function = 1 - F(x)

    cdf = 1 - ccdf;                      % CDF F(x) = 1 - ccdf(x)

end
