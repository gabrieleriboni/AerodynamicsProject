%% traccia Hess Smith (2024)

clc
close all
clear
addpath 'mat_functions'

%% Input

% imposto i paramentri di input ed effettuo le necessarie normalizzazioni

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 2;  % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

TestCase = 1;

CodiceProfilo = '0012';
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% importo i punti estremità di ogni pannello necessari per costruire il
% vettore delle x e delle y, carico la linea media dal file testo ed eseguo
% il plot del profilo attraverso i pannelli e della linea media

Corpo = importXfoilProfile('NACA_0012.dat');
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

[~, LE_index] = min(Corpo{:, 1});

figure;
plot(x, y, '.-', 'MarkerSize', 10, 'Color', '#4b808c')
axis equal
grid on

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alphai, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);

%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

    for j = 1:NPannelli
        index_j = indexStart_colonna + j;  % Colonna

        Estremo_1_qui = Estremo_1(j, :)';
        Estremo_2_qui = Estremo_2(j, :)';

        L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
        G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

        matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

        matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);
    end
end


%% Creazione delle componenti dei vettori a_v, c_s e c_v

Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';

b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;


%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));


%% Risoluzione sistema lineare

Soluzione = matriceA\TermineNoto;


%% Calcolo del cp e della velocità sui pannelli

% inizializzo le matrici che conterranno le velocità indotte da ogni
% pannello a causa della presenza delle sorgenti e dei vortici

U_sorg_pesata = zeros(2, NPannelli);
U_vort = zeros(2, NPannelli);
for i = 1:NPannelli
    Centro_qui = Centro(i, :)';
    for j = 1:NPannelli

        % calcolo gli estremi di ogni pannello nel sistema di riferimento
        % locale

        Estremo_1_qui = Estremo_1(j, :)';
        Estremo_2_qui = Estremo_2(j, :)';
        L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
        G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));
        U_sorg_pesata(:, i) = U_sorg_pesata(:, i) + Soluzione(j) * ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
        U_vort(:, i) = U_vort(:, i) + ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
    end
    U_vort_pesata(:, i) = U_vort(:, i) * Soluzione(end);

    % dopo aver calcolato le velocità pesate per le rispettive intensità di
    % sorgenti e vortici, sommo per ogni pannello tutti i contributi
    % considerando anche la corrente asintotica

    U(:, i) = U_inf + U_sorg_pesata(:, i) + U_vort_pesata(:, i);
    U_tang = dot(U(:, i), Tangente(i, :)');

    % calcolo il coefficiente di pressione sul singolo pannello, dopo aver
    % calcolato la componente tangenziale della velocità su ognuno di essi

    Cp(i) = 1 - U_tang^2 / (norm(U_inf))^2;
end

% suddivido il plot tra dorso e ventre, spezzando la pannellizzazione nel
% momento in cui si passa dal bordo di attacco

figure
plot(Centro(LE_index:end, 1), -Cp(LE_index:end))
hold on
plot(Centro(1:LE_index, 1), -Cp(1:LE_index))
grid on
legend("Upper", "Lower")

%% Calcolo Cl

% per calcolare il Cl, integro il Cp sul singolo pannello. Ciò si
% semplifica ad una semplice somma calcolata mediante un ciclo for

Cl = 0;
for i = 1:NPannelli
    Cl = Cl - 1 / Chord * Cp(i) * lunghezza(i) * dot(Normale(i, :), U_inf_normal);
end
Cl

%% Calcolo Cm,AC e confronto con XFOIL

% per calcolare il coefficiente di momento rispetto al centro aerodinamico
% calcolo dapprima la matrice R dei vettori posizione del centro
% aerodinamico rispetto al sistema di riferimento iniziale. Successivamente
% calcolo il Cm,AC come semplice somma dei vari contributi dei pannelli

R = [];
for i = 1:NPannelli
    R(i, :) = Centro(i, :)' - [Chord/4; 0];
end

R = [R, zeros(NPannelli, 1)];
Normale_Cm = [Normale, zeros(NPannelli, 1)];

Cm_AC = 0;
for i = 1:NPannelli
    Cm_AC = Cm_AC + 1 / Chord^2 * Cp(i) * lunghezza(i) * dot(cross(R(i, :), Normale_Cm(i, :)), [0; 0; 1]);
end
Cm_AC

%% Calcolo linea media

Camber = linspace(0, 1, 100)';
Camber(:, 2) = zeros(100, 1);
CamberPanels = length(Camber(:, 1));

% Linea media in funzione di eta:

for i = 1:CamberPanels-1
    Camber_eta(i, :) = acos(-2 / Chord * (Camber(i, :) - Chord/2));
end

dy = zeros(CamberPanels-1, 1);
h = (Camber_eta(end, 1) - Camber_eta(1, 1)) / (CamberPanels-1);
dy(1) = (-3 * Camber_eta(1, 2) + 4 * Camber_eta(2, 2) - Camber_eta(3, 2)) / (2 * h);
dy(end) = (Camber_eta(CamberPanels-3, 2) - 4 * Camber_eta(CamberPanels-2, 2) + 3 * Camber_eta(LE_index-1, 2)) / (2 * h);

for i = 2:CamberPanels-2
    dy(i) = (Camber_eta(i+1, 2) - Camber_eta(i-1, 2)) / (2 * h);
end


%% Calcolo angolo di progetto
alpha_th = 1/pi * trapcomp(Camber_eta(:, 1), dy)


%% Calcolo angolo di incidenza per portanza nulla
for i = 1:LE_index-1
    cos_eta(i) = -2 / Chord * (Camber(i, 1) - Chord/2);
    f(i) = dy(i) * cos_eta(i);
end
alpha_0 = alpha_th - 1/pi * trapcomp(Camber_eta(:, 1), f)

%%

figure
subplot(2, 1, 1)
plot(Centro(LE_index:end, 1), -Cp(LE_index:end))
hold on
plot(Centro(1:LE_index, 1), -Cp(1:LE_index))
grid on
legend("Upper", "Lower")
xlabel("x/c [-]")
ylabel("Cp [-]")

subplot(2, 1, 2)
plot(x, y, '.-', 'MarkerSize', 10, 'Color', '#4b808c')
hold on
plot(Camber(:, 1), Camber(:, 2))
axis equal
grid on
xlabel("x/c [-]")