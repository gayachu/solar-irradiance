clc;
clear;
close all;

%% ----------------------------
% 1. REAL MONTHLY DATA (kWh/m^2)
% -----------------------------
monthly_real = [ ...
    126 141 163 171 170 116 86 124 142 152 138 89 ...
];  % 

%% ----------------------------
% 2. CONSTANTS
% -----------------------------
omega_season = 2*pi/365;
omega_diurnal = 2*pi/24;

phi_season = -7*pi/9;
phi_diurnal = -pi/2;

I0 = 1000;

days_month = [31 28 31 30 31 30 31 31 30 31 30 31];
hours_month = days_month * 24;

t = 0:1:(365*24 - 1);   % hourly time vector

%% ----------------------------
% 3. OBJECTIVE FUNCTION 
% -----------------------------
objective = @(params) ...
    sum( ...
    ( monthly_real - ...
    compute_monthly(params, t, omega_season, omega_diurnal, ...
                    phi_season, phi_diurnal, I0, hours_month) ...
    ).^2 );

%% ----------------------------
% 4. INITIAL GUESS
% -----------------------------
% [Cavg, Aseason, beta, Adiurnal]
p0 = [0.5, 0.4, 1, 0.9];

%% ----------------------------
% 5. OPTIMISE
% -----------------------------
params_opt = fminsearch(objective, p0);

disp("Optimised Parameters:")
disp(params_opt)

%% ----------------------------
% 6. COMPUTE FINAL MODEL
% -----------------------------
monthly_model = compute_monthly(params_opt, t, omega_season, ...
                                omega_diurnal, phi_season, ...
                                phi_diurnal, I0, hours_month);

error = monthly_real - monthly_model;
SSE = sum(error.^2);
RMSE = sqrt(mean(error.^2));

disp("Final SSE:")
disp(SSE)

disp("Final RMSE:")
disp(RMSE)

%% ----------------------------
% 7. PLOT
% -----------------------------
months = 1:12;
month_labels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

figure
plot(months, monthly_real, 'o-', 'LineWidth',2)
hold on
plot(months, monthly_model, 's-', 'LineWidth',2)

xticks(months)
xticklabels(month_labels)

xlabel('Month')
ylabel('Monthly Insolation (kWh/m^2)')
title('Least Squares Fit')
legend('Observed','Model')
grid on

%% ----------------------------
% Helper block 
% -----------------------------
function monthly_model = compute_monthly(params, t, ...
    omega_season, omega_diurnal, ...
    phi_season, phi_diurnal, I0, hours_month)

    I_model = I0 .* (1 - params(3) .* ...
        (params(1) + params(2).*sin(omega_season.*(t/24) + phi_season))) ...
        .* (params(4).*sin(omega_diurnal.*t + phi_diurnal));

    I_model(I_model < 0) = 0;

    Energy_hourly = I_model / 1000;

    monthly_model = zeros(1,12);
    start_idx = 1;

    for m = 1:12
        end_idx = start_idx + hours_month(m) - 1;
        monthly_model(m) = sum(Energy_hourly(start_idx:end_idx));
        start_idx = end_idx + 1;
    end
end
