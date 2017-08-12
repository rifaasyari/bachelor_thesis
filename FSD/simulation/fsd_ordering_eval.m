% Effect of the FSD ordering on the Matrix U

addpath('../algorithm');

clear
clc

N_r = 2;
N_t = 2;

P = 4;  % 4-QAM
l_p = ceil(sqrt(N_t)-1);  % Valid iff N_r == N_t
l_1 = N_t - l_p;
n_S = [ones(l_1,1); P*ones(l_p,1)];

trials = 1e4;
diag_mean = zeros(N_t,1);
diag_mean_o = zeros(N_t,1);

for i = 1:trials
    H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
    H_o = fsd_ordering(H,P,n_S);
    [~,U] = qr(H);
    [~,U_o] = qr(H_o);
    diag_mean = diag_mean + diag(U).^2;
    diag_mean_o = diag_mean_o + diag(U_o).^2;
end

diag_mean = diag_mean / trials;  % E[u_{ii}^2] = N_r - i + 1
diag_mean_o = diag_mean_o / trials;
expected_mean = N_r*ones(N_t,1) - (1:N_t)' + ones(N_t,1);

fprintf('Mean of u_ii^2 without ordering: \n\t Expected: [ %s] \n\t Actual: [ %s]\n', ... 
    sprintf('%d ',expected_mean), sprintf('%.3f ', diag_mean));

if N_t == 2
    expected_mean_o = [2.75 0.625];
    fprintf('Mean of u_ii^2 with FSD ordering for 2 x 2 system: '); 
    fprintf('\n\t Expected: [ %s] \n\t Actual: [ %s]\n', ... 
    sprintf('%.3f ',expected_mean_o), sprintf('%.3f ', diag_mean_o));
else
    fprintf('Mean of u_ii^2 with FSD ordering for %d x %d system: ', ...
        N_r, N_t);
    fprintf('\n\t Expected: Increase in levels 1:%d, Decrease in levels %d:%d', ...
        l_1, l_1+1,N_t);
    fprintf('\n\t Actual: [ %s]\n', sprintf('%.3f ', diag_mean_o));
end

% Effect of the FSD ordering on the Outage Probability 

level_step = 2;
trials = 1e4;

P = 4;  % 4-QAM
l_p = ceil(sqrt(N_t)-1);  % Valid iff N_r == N_t
l_1 = N_t - l_p;
n_S = [ones(l_1,1); P*ones(l_p,1)];

signal_level = -40:2:10;

if N_t == 2

    pos = 1;

    outage_1st_anal = zeros(1,length(signal_level));
    outage_2nd_anal = zeros(1,length(signal_level));
    outage_1st_sim = zeros(1,length(signal_level));
    outage_2nd_sim = zeros(1,length(signal_level));

    for i = signal_level
        x = 10^(i/10);
        outage_1st_anal(pos) = 1 - (1+x/2)*exp(-2*x);
        outage_2nd_anal(pos) = 1 - 2*(1+x)*exp(-x) + (1+x)^2*exp(-2*x);

        for j = 1:trials
            H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
            h1 = H(:,1)/norm(H(:,1));
            h2 = H(:,2)/norm(H(:,2));
            phi = subspace(h1,h2);
            eta_2 = sin(phi)^2*min(norm(H(:,1))^2,norm(H(:,2))^2);
            eta_1 = max(norm(H(:,1))^2,norm(H(:,2))^2);
            if eta_2 < x
                outage_1st_sim(pos) = outage_1st_sim(pos) + 1;
            end
            if eta_1 < x
                outage_2nd_sim(pos) = outage_2nd_sim(pos) + 1;
            end
        end

        pos = pos + 1;
    end

    outage_1st_sim = outage_1st_sim / trials;
    outage_2nd_sim = outage_2nd_sim / trials;

    figure
    semilogy(signal_level,outage_1st_anal, 'k-', ...
        'DisplayName', 'FSD 1st step anal.');
    hold on
    semilogy(signal_level,outage_1st_sim, 'kx', ...
        'DisplayName', 'FSD 1st step sim.');
    hold on
    semilogy(signal_level,outage_2nd_anal, 'k--', ...
        'DisplayName', 'FSD 2nd step anal.');
    hold on
    semilogy(signal_level,outage_2nd_sim, 'ko', ...
        'DisplayName', 'FSD 2nd step sim.');
    xticks(-40:10:10);
    ylim([1e-4 1]);
    title(sprintf('$N_r=N_t=%d$', N_r));
    legend('show', 'location', 'northwest');
    grid on
    
    savefig(sprintf('outage_%s.fig', datestr(datetime('now'), 'dd-mm-yyyy_HH-MM-SS')));
    close
end
