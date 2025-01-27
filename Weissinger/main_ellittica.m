clear
close all
clc

%%

config.RootChord = 1.625;
config.DihedralAngle = 0; % [°]
config.SweepAngle = 0; % [°]
config.TaperRatio = 0.672;
config.AspectRatio = 7.32;
config.Span = 10.92;
config.LEPosition_X = 0;
config.LEPosition_Y = 0;
config.LEPosition_Z = 0;

config.RotationAngle_X = 0;
config.RotationAngle_Y = 3.5;
config.RotationAngle_Z = 0;

% Discretization options
config.SemiSpanwiseDiscr = 10;
config.ChordwiseDiscr = 10;

config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%%

alpha_vect = linspace(-5, 10, 20);
alpha_0 = 0;

S = config.Surface;
b = config.Span;
c_0 = 4*S / (pi * b);

lambda = b^2/S;

N = 100;
theta = pi / N * (1:N) - pi / (2*N);

c = @(theta) c_0 * sin(theta);

Cl_alpha_elliptic = 2*pi;

CL_vect = [];
CD_vect = [];

for alpha = alpha_vect
    A = zeros(N, N);
    d = zeros(N, 1);
    for i = 1:N
        for j = 1:N
            A(i, j) = -4 * b * sin(j * theta(i)) / c(theta(i)) - Cl_alpha_elliptic * j * sin(j * theta(i)) / sin(theta(i));
        end
        d(i) = Cl_alpha_elliptic * deg2rad(alpha);
    end

    B = A\d;
    B1 = B(1);
    alpha_i = B1;

    CL_elliptic = - pi * lambda * B1;
    CD_elliptic = - CL_elliptic * B1;

    CL_vect = [CL_vect; CL_elliptic];
    CD_vect = [CD_vect; CD_elliptic];
end



%%

figure
subplot(1, 2, 1)
plot(alpha_vect, CL_vect, '.-', 'MarkerSize', 8)
grid on
xlabel("\alpha [°]")
ylabel("C_L [-]")

subplot(1, 2, 2)
plot(CD_vect, CL_vect, '.-', 'MarkerSize', 8)
grid on
xlabel("C_D [-]")
ylabel("C_L [-]")

