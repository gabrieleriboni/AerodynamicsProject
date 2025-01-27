clear
close all
clc
%%
% Noto l'angolo di Theodorsen esatto del NACA 2408 e la sua linea media
% analitica, sono stati validati i metodi di derivazione e di integrazione
% per la stima di tale angolo

Chord = 1;


Corpo = importXfoilProfile('NACA_2408_fit.dat', 1, inf);
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

[~, LE_index] = min(x);
NPannelli = length(x);

% Calcolo linea media
Camber = flipud(load("NACA_2408_camber_fit.dat"));

% Calcolo con interpolante
Camber_fit = spline(Camber(:, 1), Camber(:, 2));
x_eval = linspace(0, Chord, 1001);

m = 0.02;
p = 0.4;

% Funzione analitica della linea media
y_vera = @(x) m / p^2 * (2 * p * x - x.^2) .* (x <= p) + m / (1 - p)^2 * (1 - 2 * p + 2 * p * x - x.^2) .* (x > p);
dy_vera = @(x) 2 * m / p^2 * (p - x) .* (x <= p) + 2 * m / (1 - p)^2 * (p - x) .* (x > p);
phi = acos(1 - 2*p);
dy_vera_eta = @(eta) 2 * m / p^2 * (p - 1/2 * (1 - cos(eta))) .* (eta <= phi) + 2 * m / (1 - p)^2 * (p - 1/2 * (1 - cos(eta))) .* (eta > phi);

dy = fnder(Camber_fit);
dy_plot = ppval(dy, x_eval);

% plot camber
figure
plot(x_eval, ppval(Camber_fit, x_eval));
hold on
plot(x_eval, y_vera(x_eval))
grid on
legend("spline", "vera")
title("Plot y_{lm}(x) con spline semplice")

% plot derivate
figure
plot(x_eval, dy_vera(x_eval))
hold on
plot(x_eval, dy_plot)
grid on
legend("Vera", "Approx")
title("Plot y'_{lm}(x) con spline semplice")


alpha_th = [];

for i = 10:5:1000
    eta = linspace(0, pi, i);

    x_eta = Chord/2 * (1 -cos(eta));

    dy_eta = ppval(dy, x_eta);

    alpha_th = [alpha_th; 1/pi * trapcomp(eta, dy_vera_eta(eta))];
end

% figure
% plot(10:5:1000, alpha_th)
% grid on

alpha_th(end)


%% Check alpha_0 con XFOIL

f = dy_vera_eta(eta) .* cos(eta);

alpha_0 = rad2deg(alpha_th(end) - 1 / pi * trapcomp(eta, f))

%% Ricalcolo linea media
x_eval_plot = linspace(0, Chord, 1001);
err = [];
PanelRatio = [];
y(1) = 0;
y(end) = 0;
for CamberPanels = 50000
    Body = [];
    x_eval = flipud(linspace(0, Chord, CamberPanels/2)')';
    x_eval = [x_eval, flipud(x_eval')'];
    x_eval(CamberPanels/2) = [];

    LEPanels = 1/20 * CamberPanels;
    Body(:, 1) = x_eval';
    Body(:, 2) = [flipud(interp1(x(1:LE_index), y(1:LE_index),  linspace(0, Chord, CamberPanels/2), 'pchip')')',  interp1(x(LE_index+1:end), y(LE_index+1:end), linspace(2 * Chord / CamberPanels, Chord, CamberPanels/2 - 1), 'pchip')];

    clear Camber
    Camber(:, 1) = Body(1:CamberPanels/2, 1);
    for  i = 1:CamberPanels/2
        Camber(i, 2) = (Body(i, 2) + Body(end - i + 1, 2)) / 2;
    end

    figure
    plot(Body(:, 1), Body(:, 2), '.-')
    grid on
    axis equal
    hold on
    plot(Camber(:, 1), Camber(:, 2), '.-')
    title("Corpo con camber line calcolata")


    Camber_fit = csaps(Camber(:, 1), Camber(:, 2), 0.9999998);
    dy = fnder(Camber_fit);
    
    eta = linspace(0, pi, 1001);
    x_eta = Chord/2 * (1 -cos(eta));

    dy_eta = ppval(dy, x_eta);

    alpha_th_stimato = 1/pi * trapcomp(eta, dy_eta)
    alpha_th_vero = 1/pi * trapcomp(eta, dy_vera_eta(eta))
    err = [err, (alpha_th_stimato - alpha_th_vero)/ alpha_th_vero];
    PanelRatio = [PanelRatio, CamberPanels / NPannelli];
end

[errmin, index] = min(abs(err));
PanelRatio(index)
index
errmin



figure
% plot(x_eval_plot, ppval(Camber_fit, x_eval_plot));
plot(Camber(:, 1), Camber(:, 2));
hold on
plot(x_eval_plot, y_vera(x_eval_plot))
grid on
legend("spline", "vera")
title("Plot y_{lm}(x) con interp1")

% plot derivate
figure
plot(x_eval_plot, dy_vera(x_eval_plot))
hold on
plot(x_eval_plot, ppval(dy, x_eval_plot))
grid on
legend("Vera", "Approx")
title("Plot y'_{lm}(x) con interp1 e csaps")