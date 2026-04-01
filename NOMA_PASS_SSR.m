%% ============================================================
%%  NOMA-Assisted Downlink Pinching-Antenna Systems (PASS)
%%  Secrecy Sum Rate Maximization via Alternating Optimization
%%
%%  Replicates: WCL2026-0016
%%  "Secrecy Sum Rate for NOMA-Assisted Downlink
%%   Pinching-Antenna Systems"
%%  Tian et al., IEEE Wireless Communications Letters, 2026
%%
%%  Figures reproduced:
%%    Fig. 1 – Convergence of the AO algorithm
%%    Fig. 2 – SSR vs transmit power (three schemes)
%%    Fig. 3 – PSO vs 1-D exhaustive search
%%
%%  Requires: CVX toolbox  http://cvxr.com/cvx/
%%            MATLAB R2019b or later
%% ============================================================

%% ============================================================
%%  CVX SETUP FOR MATLAB ONLINE
%%
%%  MATLAB Online does not ship CVX by default. Follow these
%%  steps before running this script:
%%
%%  Step 1 – Download CVX from http://cvxr.com/cvx/download/
%%           Choose the "Standard" (free) bundle.
%%           Upload the unzipped folder to MATLAB Drive so it
%%           appears at, e.g.:   /MATLAB Drive/cvx/
%%
%%  Step 2 – Uncomment and edit the line below so the path
%%           matches the location you uploaded CVX to, then
%%           run that line once to initialise CVX for the
%%           current MATLAB Online session:
%%
%%       run('/MATLAB Drive/cvx/cvx_setup.m')
%%
%%  Alternatively you can paste that single line into the
%%  MATLAB Online Command Window and press Enter.
%%
%%  Step 3 – The setup only needs to be done ONCE per session.
%%           After that, just run this script normally.
%%
%%  Quick one-liner to copy into the Command Window:
%%       run('/MATLAB Drive/cvx/cvx_setup.m')
%% ============================================================

% ---------- Uncomment the line below and adjust the path -----
% run('/MATLAB Drive/cvx/cvx_setup.m');
% -------------------------------------------------------------

clear; clc; close all;
rng(42);            % fixed seed for reproducibility

%% ======================== PARAMETERS ========================

c          = 3e8;                               % speed of light (m/s)
fc         = 28e9;                              % carrier frequency (28 GHz)
sigma2_dBm = -97;                               % noise power (dBm)
sigma2     = 10^((sigma2_dBm - 30) / 10);       % noise power (Watts)

d  = 3;     % waveguide height above ground (m)
Dx = 60;    % service-area x-span (m)
Dy = 20;    % service-area y-span (m)
L  = Dx;    % waveguide length equals x-span (m)

N = 4;      % number of legitimate users
K = 2;      % number of eavesdroppers

% Minimum achievable-rate requirement per user (bits/s/Hz) – Eq. (13c)
Rth = 0.3 * ones(1, N);

% Alternating optimisation parameters
MAX_AO  = 30;   % outer AO iterations  (paper uses 30 in Fig. 2/3)
MAX_SCA = 20;   % SCA inner iterations
PSO_S   = 30;   % PSO swarm size
PSO_T   = 50;   % PSO maximum iterations

% Pre-factor in the free-space path-loss channel gain
% h = G / distance^2    [Eqs. (1)-(2)]
G = c^2 / (16 * pi^2 * fc^2);

%% =================== RANDOM TOPOLOGY =======================
% Legitimate users: uniform in [0,Dx] x [-Dy/2, Dy/2], z = 0
xu = Dx * rand(1, N);
yu = Dy * (rand(1, N) - 0.5);

% Eavesdroppers: same service area
xe = Dx * rand(1, K);
ye = Dy * (rand(1, K) - 0.5);

fprintf('=== NOMA-PASS Secrecy Sum Rate Simulation ===\n');
fprintf('N=%d users, K=%d eavesdroppers\n', N, K);
fprintf('d=%.1f m, L=%.0f m, fc=%.0f GHz, sigma2=%.0f dBm\n', ...
        d, L, fc/1e9, sigma2_dBm);

%% ==================== FIGURE 1: CONVERGENCE ================
fprintf('\n[Fig. 1] Convergence of the AO algorithm\n');

P_fig1_dBm = [10, 20, 30];          % transmit power levels
SSR_conv   = zeros(length(P_fig1_dBm), MAX_AO);

for pp = 1:length(P_fig1_dBm)
    P_W = dBm2W(P_fig1_dBm(pp));

    % ---- initialise ----
    xAnt  = L / 2;               % PA starts at waveguide mid-point
    alpha = ones(1, N) / N;      % equal power allocation
    xu_s  = xu;  yu_s = yu;      % will be re-sorted each iteration

    for r = 1:MAX_AO
        % Step 4 of Algorithm 1: update channel gains
        h_all  = chanGain(xAnt, xu_s, yu_s, G, d);  % 1xN
        he_all = chanGain(xAnt, xe,   ye,   G, d);  % 1xK

        % Sort users in ascending channel-gain order: h1 <= h2 <= ... <= hN
        [h_s, idx_s] = sort(h_all, 'ascend');
        xu_s = xu_s(idx_s);
        yu_s = yu_s(idx_s);

        he_w = max(he_all);     % worst (strongest) eavesdropper channel

        % Step 5 (Algorithm 1): power allocation via SCA/CVX  [Prob. P3]
        alpha = SCA_power(h_s, he_w, P_W, sigma2, N, Rth, MAX_SCA);

        % Step 6 (Algorithm 1): PA position via PSO            [Prob. P4]
        xAnt = PSO_position(xu_s, yu_s, xe, ye, ...
                             alpha, P_W, sigma2, G, d, L, PSO_S, PSO_T);

        % Evaluate and record SSR
        h_new  = chanGain(xAnt, xu_s, yu_s, G, d);
        he_new = chanGain(xAnt, xe,   ye,   G, d);
        SSR_conv(pp, r) = totalSSR(h_new, max(he_new), alpha, P_W, sigma2, N);
    end

    fprintf('  P = %2d dBm  =>  converged SSR = %.3f bits/s/Hz\n', ...
            P_fig1_dBm(pp), SSR_conv(pp, end));
end

figure(1); clf;
line_styles = {'b-o', 'r-s', 'g-^'};
for pp = 1:length(P_fig1_dBm)
    plot(1:MAX_AO, SSR_conv(pp,:), line_styles{pp}, ...
         'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', 1:2:MAX_AO);
    hold on;
end
xlabel('Iterations',                         'FontSize', 12);
ylabel('Total Secrecy Sum Rate (bits/s/Hz)', 'FontSize', 12);
title('Fig. 1: Convergence of AO Algorithm (N=4, K=2)', 'FontSize', 11);
legend(arrayfun(@(p) sprintf('P = %d dBm', p), P_fig1_dBm, ...
       'UniformOutput', false), 'Location', 'southeast', 'FontSize', 10);
grid on;  xlim([1, MAX_AO]);

%% ========== FIGURE 2: SSR vs TRANSMIT POWER (THREE SCHEMES) ===========
fprintf('\n[Fig. 2] SSR vs transmit power\n');

P_sweep = 0:5:30;       % dBm
nPts    = length(P_sweep);

SSR_proposed = zeros(1, nPts);
SSR_fixedAnt = zeros(1, nPts);
SSR_noopt    = zeros(1, nPts);

for pidx = 1:nPts
    P_W = dBm2W(P_sweep(pidx));

    %--- Scheme 1: Proposed joint CVX + PSO (AO) ---
    xAnt_j  = L / 2;
    alpha_j = ones(1, N) / N;
    xu_j    = xu;   yu_j = yu;

    for r = 1:MAX_AO
        h_j  = chanGain(xAnt_j, xu_j, yu_j, G, d);
        he_j = chanGain(xAnt_j, xe,   ye,   G, d);

        [h_js, idx_j] = sort(h_j, 'ascend');
        xu_j = xu_j(idx_j);   yu_j = yu_j(idx_j);
        he_wj = max(he_j);

        alpha_j = SCA_power(h_js, he_wj, P_W, sigma2, N, Rth, MAX_SCA);
        xAnt_j  = PSO_position(xu_j, yu_j, xe, ye, ...
                               alpha_j, P_W, sigma2, G, d, L, PSO_S, PSO_T);
    end
    h_j  = chanGain(xAnt_j, xu_j, yu_j, G, d);
    he_j = chanGain(xAnt_j, xe,   ye,   G, d);
    SSR_proposed(pidx) = totalSSR(h_j, max(he_j), alpha_j, P_W, sigma2, N);

    %--- Scheme 2: Fixed antenna at waveguide centre + CVX ---
    xAnt_f = L / 2;
    h_f    = chanGain(xAnt_f, xu, yu, G, d);
    he_f   = chanGain(xAnt_f, xe, ye, G, d);
    [h_fs, ~] = sort(h_f, 'ascend');
    alpha_f   = SCA_power(h_fs, max(he_f), P_W, sigma2, N, Rth, MAX_SCA);
    SSR_fixedAnt(pidx) = totalSSR(h_fs, max(he_f), alpha_f, P_W, sigma2, N);

    %--- Scheme 3: No optimisation (fixed PA at L/3, equal power) ---
    xAnt_n = L / 3;
    h_n    = chanGain(xAnt_n, xu, yu, G, d);
    he_n   = chanGain(xAnt_n, xe, ye, G, d);
    [h_ns, ~] = sort(h_n, 'ascend');
    SSR_noopt(pidx) = totalSSR(h_ns, max(he_n), ones(1,N)/N, P_W, sigma2, N);

    fprintf('  P=%2d dBm  Proposed=%.2f  Fixed=%.2f  NoOpt=%.2f\n', ...
            P_sweep(pidx), SSR_proposed(pidx), SSR_fixedAnt(pidx), SSR_noopt(pidx));
end

figure(2); clf;
plot(P_sweep, SSR_proposed, 'b-s',  'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(P_sweep, SSR_fixedAnt, 'r--^', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(P_sweep, SSR_noopt,    'k:o',  'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Transmit Power (dBm)',                'FontSize', 12);
ylabel('Total Secrecy Sum Rate (bits/s/Hz)', 'FontSize', 12);
title('Fig. 2: SSR vs Transmit Power',        'FontSize', 11);
legend({'Proposed: CVX+PSO', 'Fixed Antenna + CVX', 'No Optimisation'}, ...
       'Location', 'northwest', 'FontSize', 10);
grid on;

%% ============== FIGURE 3: PSO vs 1-D EXHAUSTIVE SEARCH =================
fprintf('\n[Fig. 3] PSO vs 1-D exhaustive search\n');

P_sweep3  = 0:5:30;
nPts3     = length(P_sweep3);
SSR_PSO_f = zeros(1, nPts3);
SSR_1D_f  = zeros(1, nPts3);

for pidx = 1:nPts3
    P_W3 = dBm2W(P_sweep3(pidx));

    % Obtain a reference power allocation at the centre position
    xAnt_0 = L / 2;
    h_0    = chanGain(xAnt_0, xu, yu, G, d);
    he_0   = chanGain(xAnt_0, xe, ye, G, d);
    [h_0s, idx_0] = sort(h_0, 'ascend');
    xu_0   = xu(idx_0);   yu_0 = yu(idx_0);
    alpha_0 = SCA_power(h_0s, max(he_0), P_W3, sigma2, N, Rth, MAX_SCA);

    %--- PSO positioning ---
    xAnt_p = PSO_position(xu_0, yu_0, xe, ye, ...
                          alpha_0, P_W3, sigma2, G, d, L, PSO_S, PSO_T);
    h_p    = chanGain(xAnt_p, xu_0, yu_0, G, d);
    he_p   = chanGain(xAnt_p, xe,   ye,   G, d);
    SSR_PSO_f(pidx) = totalSSR(h_p, max(he_p), alpha_0, P_W3, sigma2, N);

    %--- 1-D exhaustive search over 500 candidate positions ---
    xGrid    = linspace(0, L, 500);
    ssr_grid = zeros(1, 500);
    for jj = 1:500
        h_g  = chanGain(xGrid(jj), xu_0, yu_0, G, d);
        he_g = chanGain(xGrid(jj), xe,   ye,   G, d);
        ssr_grid(jj) = totalSSR(h_g, max(he_g), alpha_0, P_W3, sigma2, N);
    end
    SSR_1D_f(pidx) = max(ssr_grid);

    fprintf('  P=%2d dBm  PSO=%.2f  1D-Search=%.2f\n', ...
            P_sweep3(pidx), SSR_PSO_f(pidx), SSR_1D_f(pidx));
end

figure(3); clf;
plot(P_sweep3, SSR_PSO_f, 'b-s',  'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(P_sweep3, SSR_1D_f,  'r--o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Transmit Power (dBm)',                'FontSize', 12);
ylabel('Total Secrecy Sum Rate (bits/s/Hz)', 'FontSize', 12);
title('Fig. 3: Antenna Positioning – PSO vs 1-D Exhaustive Search', 'FontSize', 11);
legend({'PSO', '1-D Exhaustive Search'}, 'Location', 'northwest', 'FontSize', 10);
grid on;

fprintf('\nSimulation complete.\n');

%% ============================================================
%%                     LOCAL FUNCTIONS
%% ============================================================

%% -- chanGain --------------------------------------------------
%  Free-space LoS channel power gain, Eqs. (1)-(2).
%  h_i = c^2/(16*pi^2*fc^2) / ( (xAnt-x_i)^2 + y_i^2 + d^2 )
function h = chanGain(xAnt, x_nodes, y_nodes, G, d)
    h = G ./ ((xAnt - x_nodes).^2 + y_nodes.^2 + d^2);
end

%% -- dBm2W -----------------------------------------------------
function W = dBm2W(dBm_val)
    W = 10^((dBm_val - 30) / 10);
end

%% -- totalSSR --------------------------------------------------
%  Total secrecy sum rate, Eq. (10).
%  h    (1xN): channel gains sorted ascending  h_1 <= ... <= h_N
%  he_w (scalar): strongest eavesdropper channel gain (worst case)
%  alpha(1xN): power coefficients aligned with h
function total = totalSSR(h, he_w, alpha, P, sigma2, N)
    total = 0;
    for i = 1:N
        % Legitimate user SINR (Eq. 6)
        % U_i = sum_{j>i} alpha_j  (residual interference after SIC)
        Ui      = sum(alpha(i+1:end));
        wi      = sigma2 / (h(i) * P);         % sigma^2/(P*h_i)
        gamma_i = alpha(i) / (Ui + wi);
        Ri      = log2(1 + gamma_i);           % Eq. (7)

        % Eavesdropper SINR for user i (Eq. 8, worst case)
        % Eavesdropper has no SIC knowledge -> all other users are interference
        Ue      = 1 - alpha(i);                % sum_{j!=i} alpha_j
        we      = sigma2 / (he_w * P);         % sigma^2/(P*h_{e,worst})
        gamma_e = alpha(i) / (Ue + we);
        Re      = log2(1 + gamma_e);           % Eq. (9)

        total = total + max(0, Ri - Re);       % Eq. (10)
    end
end

%% -- SCA_power -------------------------------------------------
%  Power allocation via Successive Convex Approximation (CVX).
%  Solves Problem P3, Eq. (17).
%
%  h   (1xN) : channel gains sorted ASCENDING
%  he_w      : max eavesdropper channel gain
%  Rth (1xN) : minimum-rate requirements (bits/s/Hz)
function alpha_out = SCA_power(h, he_w, P, sigma2, N, Rth, max_iter)

    % Noise-to-signal ratios (constants for given xAnt and P)
    w   = sigma2 ./ (h * P);       % 1xN  w_i = sigma^2/(P*h_i)
    we  = sigma2  / (he_w * P);    % scalar  w_e = sigma^2/(P*h_{e,worst})
    gth = (2.^Rth - 1)';           % Nx1  gamma^th_i = 2^{R^th_i} - 1

    % Interference matrix: M(i,j)=1 for j>i, so U = M*alpha gives
    % U_i = sum_{j>i} alpha_j for each user i.
    M = triu(ones(N, N), 1);

    % Initialise with equal power
    alpha_k = ones(N, 1) / N;

    for iter = 1:max_iter
        % Linearisation point for the SCA lower bound (Eq. 16):
        %   r_tilde_i = U_i^{prev} + w_i
        U_prev  = M * alpha_k;
        r_tilde = max(U_prev + w(:), 1e-15);   % Nx1, avoid log(0)

        cvx_begin quiet
            variable a(N, 1)

            % U_a(i) = sum_{j>i} a_j  (affine in a)
            U_a = M * a;

            % Objective: sum_i S^lb_i (Eq. 17a)
            % S^lb_i = log2(U_i+a_i+w_i) + log2(1-a_i+w_e)
            %          - (U_i+w_i)/(log2 * r_tilde_i) + constants
            % All three terms: concave + concave + affine -> concave.
            % Constants (-log2(r_tilde_i), -log2(1+w_e)) omitted as they
            % do not affect the argmax.
            obj = 0;
            for ii = 1:N
                obj = obj ...
                    + log(U_a(ii) + a(ii) + w(ii)) / log(2) ...
                    + log(1 - a(ii) + we)           / log(2) ...
                    - (U_a(ii) + w(ii)) / (log(2) * r_tilde(ii));
            end

            maximize(obj)

            subject to
                a >= 0                                  % Eq. (17b)
                a <= 1                                  % Eq. (17b)
                sum(a) == 1                             % Eq. (17c)
                for ii = 1:N
                    % QoS: alpha_i >= gamma^th_i * (U_i + w_i)   Eq. (17d)
                    a(ii) >= gth(ii) * (U_a(ii) + w(ii))
                end
        cvx_end

        if contains(lower(cvx_status), 'solved')
            alpha_k = double(a);
        else
            warning('CVX: %s at SCA iter %d – keeping previous solution', ...
                    cvx_status, iter);
            break;
        end
    end

    alpha_out = alpha_k';   % return as 1xN
end

%% -- PSO_position ----------------------------------------------
%  Particle Swarm Optimisation for the PA position (Problem P4).
%  Maximises totalSSR over xAnt in [0, L].
%
%  xu_s, yu_s: user coordinates aligned with the current alpha ordering
function xAnt_opt = PSO_position(xu_s, yu_s, xe, ye, ...
                                  alpha, P, sigma2, G, d, L, S, T)

    N = length(xu_s);

    % Objective: SSR as a function of the scalar xAnt
    obj = @(xA) ssr_xAnt(xA, xu_s, yu_s, xe, ye, ...
                          alpha, P, sigma2, N, G, d);

    % PSO hyper-parameters (standard values)
    w_in = 0.7;   c1 = 1.5;   c2 = 1.5;

    % Initialise swarm uniformly in [0, L]
    pos = L * rand(S, 1);
    vel = (L / 4) * (2 * rand(S, 1) - 1);

    fit = arrayfun(obj, pos);

    pbest_pos = pos;
    pbest_fit = fit;
    [gbest_fit, gbest_idx] = max(fit);
    gbest_pos = pos(gbest_idx);

    for t = 1:T
        r1 = rand(S, 1);
        r2 = rand(S, 1);

        vel = w_in * vel ...
            + c1 .* r1 .* (pbest_pos - pos) ...
            + c2 .* r2 .* (gbest_pos - pos);

        pos = pos + vel;
        pos = max(0, min(L, pos));      % enforce xAnt in [0, L]

        fit = arrayfun(obj, pos);

        % Update personal bests
        improved = fit > pbest_fit;
        pbest_pos(improved) = pos(improved);
        pbest_fit(improved) = fit(improved);

        % Update global best
        [curr_best, curr_idx] = max(pbest_fit);
        if curr_best > gbest_fit
            gbest_fit = curr_best;
            gbest_pos = pbest_pos(curr_idx);
        end
    end

    xAnt_opt = gbest_pos;
end

%% -- ssr_xAnt --------------------------------------------------
%  SSR for a candidate PA position (PSO objective function).
function val = ssr_xAnt(xAnt, xu_s, yu_s, xe, ye, ...
                         alpha, P, sigma2, N, G, d)
    h  = G ./ ((xAnt - xu_s).^2 + yu_s.^2 + d^2);
    he = G ./ ((xAnt - xe  ).^2 + ye  .^2 + d^2);
    he_w = max(he);
    val  = 0;
    for i = 1:N
        Ui      = sum(alpha(i+1:end));
        wi      = sigma2 / (h(i) * P);
        gamma_i = alpha(i) / (Ui + wi);
        Ri      = log2(1 + gamma_i);

        Ue      = 1 - alpha(i);
        we      = sigma2 / (he_w * P);
        gamma_e = alpha(i) / (Ue + we);
        Re      = log2(1 + gamma_e);

        val = val + max(0, Ri - Re);
    end
end
