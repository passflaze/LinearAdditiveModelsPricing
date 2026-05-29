function plot_mc_check(Z, x_grid, cdf, T1, T2)
% Compare MC simulated increments against the FFT-reconstructed CDF/PDF.
%
% Inputs:
%   Z      : Nsim x 1 simulated increments (from sample_from_cdf)
%   x_grid : grid passed to ccdf_AB_FFT
%   cdf    : CDF values on x_grid (from ccdf_AB_FFT)
%   T1, T2 : reset dates (years)

in_bounds = (cdf > 0) & (cdf < 1);
x_v = x_grid(find(in_bounds,1,'first') : find(in_bounds,1,'last'));
c_v = cdf(find(in_bounds,1,'first')    : find(in_bounds,1,'last'));

pdf_approx = diff(c_v) ./ diff(x_v);
x_mid      = (x_v(1:end-1) + x_v(2:end)) / 2;

% Fine bins with an edge anchored exactly on x=0, so the sharp leptokurtic
% peak is not smeared across the median. Resolution ~0.25$ matches the pdf
% curvature scale (a fixed 80-bin grid was too coarse and dropped the peak).
bw    = 0.25;
edges = -30:bw:30;          % symmetric, includes 0 as an edge

figure;

% (1) pdf, linear scale: fine bins resolve the sharp peak
subplot(1,3,1)
histogram(Z, edges, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
hold on
plot(x_mid, pdf_approx, 'r-', 'LineWidth', 1.5)
xlabel('x = X_{T_2} - X_{T_1} [$]'); ylabel('pdf')
legend('MC', 'teorica (FFT)')
title(sprintf('pdf: T_1=%.2fy \\to T_2=%.2fy', T1, T2))
xlim([-30 30]); grid on

% (2) pdf, log-y: this is where the heavy AB tails are visible (on a linear
% scale they look like a flat zero)
subplot(1,3,2)
histogram(Z, edges, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
hold on
plot(x_mid, pdf_approx, 'r-', 'LineWidth', 1.5)
set(gca, 'YScale', 'log')
xlabel('x [$]'); ylabel('pdf (log)')
legend('MC', 'teorica (FFT)')
title('pdf: code (log-y)')
xlim([-30 30]); grid on

% (3) CDF: the binning-free check -- overlap here is the rigorous validation
subplot(1,3,3)
[f_emp, x_emp] = ecdf(Z);
plot(x_emp, f_emp, 'b-', 'LineWidth', 1.2); hold on
plot(x_v, c_v, 'r--', 'LineWidth', 1.5)
xlabel('x [$]'); ylabel('F(x)')
legend('MC empirica', 'FFT teorica', 'Location', 'southeast')
xlim([-30 30]); grid on
title('CDF: MC vs modello')
end
