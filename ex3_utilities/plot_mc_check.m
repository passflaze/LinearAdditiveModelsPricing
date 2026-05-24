function plot_mc_check(Z, x_grid, cdf, T1, T2)
% Compare MC simulated increments against the FFT-reconstructed CDF/PDF.
%
% Inputs:
%   Z      : Nsim x 1 simulated increments (from sample_from_cdf)
%   x_grid : grid passed to ccdf_AB_FFT
%   cdf    : CDF values on x_grid (from ccdf_AB_FFT)
%   T1, T2 : reset dates (years)

in_bounds = (cdf >= 0) & (cdf <= 1);
x_v = x_grid(find(in_bounds,1,'first') : find(in_bounds,1,'last'));
c_v = cdf(find(in_bounds,1,'first')    : find(in_bounds,1,'last'));

pdf_approx = diff(c_v) ./ diff(x_v);
x_mid      = (x_v(1:end-1) + x_v(2:end)) / 2;

figure;

subplot(1,2,1)
histogram(Z, 80, 'Normalization', 'pdf', 'FaceAlpha', 0.5)
hold on
plot(x_mid, pdf_approx, 'r-', 'LineWidth', 1.5)
xlabel('x = X_{T_2} - X_{T_1} [$]'); ylabel('pdf')
legend('MC', 'teorica (FFT)')
title(sprintf('Incremento: T_1=%.2fy \\to T_2=%.2fy', T1, T2))
xlim([-30 30]); grid on

subplot(1,2,2)
[f_emp, x_emp] = ecdf(Z);
plot(x_emp, f_emp, 'b-', 'LineWidth', 1.2); hold on
plot(x_v, c_v, 'r--', 'LineWidth', 1.5)
xlabel('x [$]'); ylabel('F(x)')
legend('MC empirica', 'FFT teorica')
xlim([-30 30]); grid on
title('CDF: MC vs modello')
end
