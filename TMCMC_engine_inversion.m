% TMCMC_engine_inversion.m
%
% 基于 TMCMC（Transitional Markov Chain Monte Carlo）的
% 航空发动机热力循环修正参数贝叶斯后验反演
%
% 敏感参数（待反演）: eta_k, eta_t, eta_T, eta_m
% 观测量            : 比推力 R_ud [N·s/kg], 比油耗 C_ud [kg/(N·h)]
%
% 参考文献：
%   Ching, J. & Chen, Y. (2007). Transitional Markov Chain Monte Carlo
%   Method for Bayesian Model Updating, Model Class Selection, and
%   Model Averaging. J. Eng. Mech., 133(7), 816-832.
%
% 运行方式：直接在 MATLAB 中打开本文件并运行（F5）

clear; clc; close all;

%% =====================================================================
%% 1. 参数定义
%% =====================================================================

param_labels = {'eta\_k', 'eta\_t', 'eta\_T', 'eta\_m'};   % LaTeX 标签（绘图）
param_names  = {'eta_k',  'eta_t',  'eta_T',  'eta_m'};    % 普通标签（打印）
n_params     = 4;

% 均匀先验上下界
lb = [0.84,  0.86,  0.980, 0.970];
ub = [0.860, 0.920, 0.995, 0.990];

% 工况参数
cond.T_H      = 288;     % 大气总温 (K)
cond.M_flight = 0.0;     % 飞行马赫数（地面静止）
cond.m        = 10.0;    % 涵道比
cond.pi_k     = 33.0;    % 压气机压比
cond.T_g      = 1700.0;  % 涡轮前总温 (K)

% 真值（仅用于生成虚拟观测、验证反演效果）
theta_true = [0.860, 0.880, 0.988, 0.980];

% 固定参数（不参与反演）
%   lambda, eta_v, eta_tv, eta_c1, eta_c2,
%   sigma_cc, sigma_kan, sigma_kask, sigma_ks
theta_fixed = [1.03, 0.87, 0.92, 0.95, 0.94, 1.00, 0.98, 0.99, 0.96];

assert(all(theta_true >= lb) && all(theta_true <= ub), 'theta_true 超出先验范围');

fprintf('===== 工况参数 =====\n');
fprintf('  T_H=%g K  M_flight=%g  m=%g  pi_k=%g  T_g=%g K\n\n', ...
    cond.T_H, cond.M_flight, cond.m, cond.pi_k, cond.T_g);

%% =====================================================================
%% 2. 前向模型验证
%% =====================================================================

fprintf('===== 前向模型验证 =====\n');
[y_true, ~] = engine_forward(theta_true, theta_fixed, cond);
if ~all(isfinite(y_true))
    error('前向模型在 theta_true 处输出非有限值，请检查参数设置');
end
fprintf('  R_ud = %.4f  [N·s/kg]\n',       y_true(1));
fprintf('  C_ud = %.6f [kg/(N·h)]\n\n', y_true(2));

%% =====================================================================
%% 3. 生成虚拟观测数据（噪声水平 1%）
%% =====================================================================

data = generate_virtual_data(theta_true, theta_fixed, cond, 0.01, 0.01);

fprintf('===== 虚拟观测数据 =====\n');
fprintf('  R: 真值=%.4f, 观测=%.4f, sigma_R=%.4f\n', ...
    data.R_true, data.R_obs, data.sigma_R);
fprintf('  C: 真值=%.6f, 观测=%.6f, sigma_C=%.6f\n\n', ...
    data.C_true, data.C_obs, data.sigma_C);

%% =====================================================================
%% 4. TMCMC 算法配置
%% =====================================================================

opts.N          = 2000;  % 粒子数（普通电脑建议 1000~3000）
opts.COV_target = 1.0;   % 目标权重变异系数（控制阶段数）
opts.N_MCMC     = 3;     % 每粒子每阶段 MCMC 步数
opts.scale      = 0.20;  % 提议协方差缩放因子（AM 风格）
opts.max_stages = 60;    % 最大阶段上限（安全阀）

fprintf('===== TMCMC 配置 =====\n');
fprintf('  粒子数 N=%d, COV_target=%.1f, N_MCMC=%d, scale=%.2f\n\n', ...
    opts.N, opts.COV_target, opts.N_MCMC, opts.scale);

%% =====================================================================
%% 5. 运行 TMCMC
%% =====================================================================

fprintf('===== 运行 TMCMC =====\n');
results = run_tmcmc(data, theta_fixed, cond, lb, ub, opts);
fprintf('\nTMCMC 完成，共 %d 个过渡阶段\n\n', results.n_stages);

%% =====================================================================
%% 6. 后验统计结果
%% =====================================================================

fprintf('后验统计结果\n');
fprintf('%-12s %8s %10s %10s %10s %10s\n', ...
    '参数','真值','后验均值','MAP','CI95_低','CI95_高');
fprintf('%s\n', repmat('-',1,66));
for i = 1:n_params
    fprintf('%-12s %8.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        param_names{i}, theta_true(i), results.theta_mean(i), ...
        results.theta_map(i), results.theta_ci95(i,1), results.theta_ci95(i,2));
end
fprintf('\n');

%% =====================================================================
%% 7. 后验预测对比
%% =====================================================================

[y_mn, ~] = engine_forward(results.theta_mean, theta_fixed, cond);
[y_mp, ~] = engine_forward(results.theta_map,  theta_fixed, cond);

fprintf('后验预测对比\n');
fprintf('  %-18s %14s %16s\n', '','R_ud [N·s/kg]','C_ud [kg/(N·h)]');
fprintf('  %-18s %14.4f %16.6f\n','真值',       data.R_true, data.C_true);
fprintf('  %-18s %14.4f %16.6f\n','观测值',     data.R_obs,  data.C_obs);
fprintf('  %-18s %14.4f %16.6f\n','后验均值预测',y_mn(1),     y_mn(2));
fprintf('  %-18s %14.4f %16.6f\n','MAP 预测',   y_mp(1),     y_mp(2));
fprintf('\n');
fprintf('  后验均值相对误差: R=%.3f%%, C=%.3f%%\n', ...
    abs(y_mn(1)-data.R_true)/data.R_true*100, ...
    abs(y_mn(2)-data.C_true)/data.C_true*100);
fprintf('  MAP    相对误差: R=%.3f%%, C=%.3f%%\n\n', ...
    abs(y_mp(1)-data.R_true)/data.R_true*100, ...
    abs(y_mp(2)-data.C_true)/data.C_true*100);

%% =====================================================================
%% 8. 绘图
%% =====================================================================

plot_results(results, theta_true, lb, ub, param_labels);


%% =====================================================================
%%                         局部函数定义
%% =====================================================================

%% --- 分段绝热指数 ---
function kT = piecewise_kT(T_g)
    if     T_g > 1600,                 kT = 1.25;
    elseif T_g > 1400 && T_g <= 1600, kT = 1.30;
    else,                              kT = 1.33;
    end
end

%% --- 分段燃气常数 ---
function RT = piecewise_RT(T_g)
    if     T_g > 1600,                 RT = 288.6;
    elseif T_g > 1400 && T_g <= 1600, RT = 288.0;
    else,                              RT = 287.6;
    end
end

%% --- 冷却引气系数 ---
function d = delta_cooling(T_g)
    d = max(0.0, min(0.15, 0.02 + (T_g - 1200)/100*0.02));
end

%% --- 发动机热力循环前向模型 ---
function [y, aux] = engine_forward(theta_s, theta_f, cond)
    % theta_s: [eta_k, eta_t, eta_T, eta_m]
    % theta_f: [lambda, eta_v, eta_tv, eta_c1, eta_c2,
    %           sigma_cc, sigma_kan, sigma_kask, sigma_ks]
    eta_k=theta_s(1); eta_t=theta_s(2); eta_T=theta_s(3); eta_m=theta_s(4);
    lambda=theta_f(1); eta_v=theta_f(2); eta_tv=theta_f(3);
    eta_c1=theta_f(4); eta_c2=theta_f(5);
    sigma_cc=theta_f(6); sigma_kan=theta_f(7);
    sigma_kask=theta_f(8); sigma_ks=theta_f(9);

    T_H=cond.T_H; M_f=cond.M_flight; m=cond.m;
    pi_k=cond.pi_k; T_g=cond.T_g;

    y=[NaN,NaN]; aux=struct();
    try
        k=1.4; R_a=287.3;
        a_s = sqrt(k*R_a*T_H);
        if ~isfinite(a_s)||a_s<=0, return; end
        V_f = a_s*M_f;

        kT=piecewise_kT(T_g); RT_g=piecewise_RT(T_g); d=delta_cooling(T_g);

        inner = 1 + V_f^2/(2*(k/(k-1))*R_a*T_H);
        if inner<=0, return; end
        tau_v = inner^(k/(k-1));
        T_B   = T_H*(inner^k);
        if ~isfinite(T_B)||T_B<=0, return; end

        pi_r = pi_k^((k-1)/k);
        if ~isfinite(pi_r)||pi_r<1, return; end
        T_k = T_B*(1+(pi_r-1)/eta_k);
        if ~isfinite(T_k)||T_k<=0, return; end

        g_T = 3e-5*T_g - 2.69e-5*T_k - 0.003;
        if ~isfinite(g_T)||g_T<=0, return; end

        W_k  = (k/(k-1))*R_a*T_B*(pi_r-1);
        H_g  = (kT/(kT-1))*RT_g*T_g;
        if abs(H_g)<1e-6, return; end
        num_lh = 1 - W_k/(H_g*eta_k);
        den_lh = 1 - W_k/(H_g*eta_k*eta_t);
        if abs(den_lh)<1e-10, return; end
        lh = num_lh/den_lh;
        if ~isfinite(lh), return; end

        sigma_bx = sigma_cc*sigma_kan;
        e_T  = (kT-1)/kT;
        pr_d = tau_v*sigma_bx*pi_k*sigma_kask*sigma_ks;
        if pr_d<=0, return; end
        exp_t = (1/pr_d)^e_T;
        if ~isfinite(exp_t), return; end

        trm1 = (kT/(kT-1))*RT_g*T_g*(1-exp_t);
        den2 = (1+g_T)*eta_k*eta_T*eta_t*eta_m*(1-d);
        if abs(den2)<1e-10, return; end
        L_sv = lh*(trm1 - W_k/den2);
        if ~isfinite(L_sv)||L_sv<=0, return; end

        num_x = 1 + m*V_f^2/(2*L_sv*eta_tv*eta_v*eta_c2);
        den_x = 1 + (m*eta_tv*eta_v*eta_c2)/(eta_c1*lambda);
        if abs(den_x)<1e-10, return; end
        x_pc = num_x/den_x;
        if ~isfinite(x_pc)||x_pc<=0||x_pc>=1, return; end

        sq1 = 2*eta_c1*lambda*x_pc*L_sv;
        if sq1<0, return; end
        V_j1 = (1+g_T)*sqrt(sq1)-V_f;
        sq2  = 2*(1-x_pc)/m*L_sv*eta_tv*eta_v*eta_c2+V_f^2;
        if sq2<0, return; end
        V_j2 = sqrt(sq2)-V_f;

        R_ud = V_j1/(1+m) + m*V_j2/(1+m);
        if ~isfinite(R_ud)||R_ud<=0, return; end
        C_ud = 3600*g_T*(1-d)/(R_ud*(1+m));
        if ~isfinite(C_ud)||C_ud<=0, return; end

        y = [R_ud, C_ud];
        aux.T_B=T_B; aux.T_k=T_k; aux.g_T=g_T;
        aux.L_sv=L_sv; aux.x_pc=x_pc; aux.lh=lh;
    catch ME
        warning('engine_forward: %s', ME.message);
    end
end

%% --- 生成虚拟观测 ---
function data = generate_virtual_data(theta_true, theta_f, cond, nR, nC)
    [yt,~] = engine_forward(theta_true, theta_f, cond);
    if ~all(isfinite(yt)), error('generate_virtual_data: 前向模型失效'); end
    sR=nR*abs(yt(1)); sC=nC*abs(yt(2));
    data.R_true=yt(1); data.C_true=yt(2);
    data.R_obs=yt(1)+sR*randn(); data.C_obs=yt(2)+sC*randn();
    data.sigma_R=sR; data.sigma_C=sC;
end

%% --- 对数似然（无常数项） ---
function ll = log_likelihood(theta, theta_f, data, cond)
    ll=-Inf;
    [yp,~] = engine_forward(theta, theta_f, cond);
    if ~all(isfinite(yp)), return; end
    rR=(data.R_obs-yp(1))/data.sigma_R;
    rC=(data.C_obs-yp(2))/data.sigma_C;
    ll=-0.5*(rR^2+rC^2);
    if ~isfinite(ll), ll=-Inf; end
end

%% --- 求下一个 beta：二分法使权重 COV = target ---
function beta_next = find_next_beta(log_likes, beta_prev, cov_tgt)
    valid = isfinite(log_likes);
    if sum(valid)<2, beta_next=1.0; return; end
    ll = log_likes(valid);

    % 检查 beta=1 时 COV 是否已满足
    if compute_cov(ll, 1.0-beta_prev) <= cov_tgt
        beta_next=1.0; return;
    end

    lo=beta_prev; hi=1.0;
    for iter_bisect = 1:60
        mid=0.5*(lo+hi);
        if compute_cov(ll, mid-beta_prev) <= cov_tgt
            lo=mid;
        else
            hi=mid;
        end
        if hi-lo<1e-7, break; end
    end
    beta_next = min(max(lo, beta_prev+1e-6), 1.0);
end

function cv = compute_cov(ll_vec, db)
    if db<=0, cv=0; return; end
    lw  = db*ll_vec;
    lw  = lw - max(lw);
    w   = exp(lw);
    mu  = mean(w);
    if mu<=0, cv=Inf; return; end
    cv  = std(w)/mu;
end

%% --- TMCMC 主算法 ---
function results = run_tmcmc(data, theta_f, cond, lb, ub, opts)
    N       = opts.N;
    n_p     = length(lb);
    cov_tgt = opts.COV_target;
    nMCMC   = opts.N_MCMC;
    sc      = opts.scale;

    % ---- 阶段 0：从均匀先验采样 ----
    particles = bsxfun(@plus, lb, bsxfun(@times, ub-lb, rand(N, n_p)));
    log_L     = zeros(N,1);
    fprintf('  初始化粒子似然...');
    for i=1:N
        log_L(i) = log_likelihood(particles(i,:), theta_f, data, cond);
    end
    fprintf(' 完成（有效粒子 %d/%d）\n', sum(isfinite(log_L)), N);

    beta_cur   = 0;
    betas      = 0;
    log_evid   = 0;
    stage      = 0;
    acc_hist   = [];

    while beta_cur < 1.0
        stage = stage+1;
        if stage > opts.max_stages
            warning('达到最大阶段数上限 %d', opts.max_stages); break;
        end

        % ---- 求 beta_next ----
        beta_new   = find_next_beta(log_L, beta_cur, cov_tgt);
        db         = beta_new - beta_cur;

        % ---- 重要性权重 ----
        log_w  = db * log_L;
        lw_max = max(log_w(isfinite(log_w)));
        w      = exp(log_w - lw_max);
        w(~isfinite(w)) = 0;
        w_sum  = sum(w);
        w_n    = w / w_sum;

        log_evid = log_evid + log(w_sum/N) + lw_max;

        % ---- 多项式重采样 ----
        idx       = datasample(1:N, N, 'Weights', w_n);
        pts_r     = particles(idx,:);
        lL_r      = log_L(idx);

        % ---- 提议协方差（加权样本协方差 × scale^2） ----
        mu_w = (w_n' * particles);                       % 1×n_p
        dif  = bsxfun(@minus, particles, mu_w);
        Sig  = (bsxfun(@times, dif, w_n))' * dif;       % n_p×n_p（向量化）
        Sig_prop = sc^2 * Sig + 1e-10*eye(n_p);
        [L_ch, flg] = chol(Sig_prop,'lower');
        if flg~=0
            L_ch = diag(sqrt(max(diag(Sig_prop),1e-12)));
        end

        % ---- MCMC 步（目标：tempered 后验 ∝ prior × L^beta_new） ----
        pts_new = pts_r;
        lL_new  = lL_r;
        n_acc   = 0;
        for i=1:N
            th_c  = pts_r(i,:);
            llc   = beta_new * lL_r(i);          % tempered log-posterior
            for iter_mcmc=1:nMCMC
                z      = randn(n_p,1);
                th_p   = th_c + (L_ch*z)';
                if any(th_p<lb) || any(th_p>ub), continue; end
                ll_raw = log_likelihood(th_p, theta_f, data, cond);
                llp    = beta_new * ll_raw;
                if log(rand()) < llp - llc
                    th_c = th_p;  llc = llp;  n_acc = n_acc+1;
                end
            end
            pts_new(i,:) = th_c;
            lL_new(i)    = log_likelihood(th_c, theta_f, data, cond);
        end

        ar = n_acc/(N*nMCMC);
        acc_hist(end+1) = ar; %#ok<AGROW>
        fprintf('  Stage %2d: beta %.4f→%.4f  |  MCMC 接受率=%.3f\n', ...
            stage, beta_cur, beta_new, ar);

        particles = pts_new;
        log_L     = lL_new;
        beta_cur  = beta_new;
        betas(end+1) = beta_new; %#ok<AGROW>

        if beta_cur >= 1.0, break; end
    end

    % ---- 后验统计 ----
    th_mean = mean(particles,1);
    th_std  = std(particles,0,1);
    [~,bi]  = max(log_L);
    th_map  = particles(bi,:);
    ci95    = zeros(n_p,2);
    for i=1:n_p
        ci95(i,:) = quantile(particles(:,i),[0.025,0.975]);
    end

    results.particles    = particles;
    results.log_L        = log_L;
    results.theta_mean   = th_mean;
    results.theta_std    = th_std;
    results.theta_map    = th_map;
    results.theta_ci95   = ci95;
    results.stage_betas  = betas;
    results.acc_hist     = acc_hist;
    results.log_evidence = log_evid;
    results.n_stages     = stage;
end

%% --- 绘图 ---
function plot_results(results, theta_true, lb, ub, param_labels)
    pts   = results.particles;
    n_p   = size(pts,2);
    nc    = ceil(n_p/2);

    %% 图1: 边缘后验分布
    fig1 = figure('Name','TMCMC 边缘后验','Position',[50,50,1100,750]);
    for i=1:n_p
        subplot(2,nc,i);
        histogram(pts(:,i),60,'Normalization','pdf', ...
            'FaceColor',[0.25,0.55,0.90],'EdgeColor','none','FaceAlpha',0.75);
        hold on;
        xl = xline(theta_true(i),'g-','LineWidth',2.5,'DisplayName','真值');
        xm = xline(results.theta_mean(i),'r--','LineWidth',2.0,'DisplayName','后验均值');
        xp = xline(results.theta_map(i), 'm:','LineWidth',2.0,'DisplayName','MAP');
        xlim([lb(i),ub(i)]);
        xlabel(param_labels{i},'FontSize',11,'Interpreter','latex');
        ylabel('pdf','FontSize',10);
        title(sprintf('%s\n真值=%.4f  均值=%.4f  MAP=%.4f', ...
            param_labels{i},theta_true(i), ...
            results.theta_mean(i),results.theta_map(i)), ...
            'FontSize',9,'Interpreter','latex');
        grid on; grid minor;
        if i==1
            legend([xl,xm,xp],'Location','best','FontSize',8);
        end
    end
    sgtitle('TMCMC 边缘后验分布','FontSize',13);

    %% 图2: 两两联合后验（Corner plot）
    fig2 = figure('Name','TMCMC 联合后验','Position',[100,30,950,880]);
    for row=1:n_p
        for col=1:n_p
            ax = subplot(n_p, n_p, (row-1)*n_p+col);
            if row==col
                % 对角：边缘直方图
                histogram(pts(:,col),50,'Normalization','pdf', ...
                    'FaceColor',[0.25,0.55,0.90],'EdgeColor','none','FaceAlpha',0.8);
                hold on;
                xline(theta_true(col),'g-','LineWidth',2.0);
                xlim([lb(col),ub(col)]);
            elseif row>col
                % 下三角：散点图（密度着色）
                scatter(pts(:,col),pts(:,row),3,'filled', ...
                    'MarkerFaceColor',[0.3,0.55,0.88],'MarkerFaceAlpha',0.12);
                hold on;
                plot(theta_true(col),theta_true(row),'r+','MarkerSize',9,'LineWidth',2);
                xlim([lb(col),ub(col)]);
                ylim([lb(row),ub(row)]);
            else
                % 上三角：Pearson 相关系数
                rho = corr(pts(:,col),pts(:,row));
                text(0.5,0.5,sprintf('\\rho=%.3f',rho), ...
                    'Units','normalized','HorizontalAlignment','center', ...
                    'FontSize',10,'FontWeight','bold', ...
                    'Color', interp1([-1,0,1],[1 0 0;0.8 0.8 0.8;0 0 1],rho));
                set(ax,'XTick',[],'YTick',[]);
            end
            if row==n_p
                xlabel(param_labels{col},'FontSize',9,'Interpreter','latex');
            end
            if col==1
                ylabel(param_labels{row},'FontSize',9,'Interpreter','latex');
            end
            set(ax,'FontSize',8);
        end
    end
    sgtitle('TMCMC 两两联合后验分布（下三角=散点，对角=边缘，上三角=相关系数）', ...
        'FontSize',11);

    %% 图3: beta 演化 + MCMC 接受率
    fig3 = figure('Name','TMCMC 过程诊断','Position',[200,50,950,380]);
    subplot(1,2,1);
    ns = results.n_stages;
    plot(0:ns, results.stage_betas,'bo-','LineWidth',1.8,'MarkerSize',7);
    xlabel('阶段序号','FontSize',11); ylabel('\beta','FontSize',11);
    title('温度参数 \beta 演化','FontSize',12);
    ylim([0,1.05]); grid on; grid minor;

    subplot(1,2,2);
    bar(1:ns, results.acc_hist, 0.6,'FaceColor',[0.3,0.7,0.4]);
    hold on;
    yline(0.234,'r--','LineWidth',1.5,'DisplayName','最优接受率 0.234');
    xlabel('阶段序号','FontSize',11); ylabel('MCMC 接受率','FontSize',11);
    title('各阶段 MCMC 接受率','FontSize',12);
    ylim([0,1]); grid on; grid minor;
    legend('FontSize',8,'Location','best');

    sgtitle('TMCMC 过程诊断','FontSize',12);

    %% 保存图片
    try
        saveas(fig1,'TMCMC_marginal_posterior.png');
        saveas(fig2,'TMCMC_joint_posterior.png');
        saveas(fig3,'TMCMC_diagnostics.png');
        fprintf('图片已保存：TMCMC_marginal_posterior.png | TMCMC_joint_posterior.png | TMCMC_diagnostics.png\n');
    catch
        fprintf('（图片保存失败，可能为无显示环境）\n');
    end
end
