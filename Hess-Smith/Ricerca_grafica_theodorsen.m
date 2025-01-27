clear
close all
clc

alpha_vect = linspace(-1, 1, 30);


Chord = 1;



LE_X_Position = 0;
LE_Y_Position = 0;

Corpo = importXfoilProfile('FX_63_143_fit.dat', 1, inf);
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

[~, LE_index] = min(x);

NPannelli = length(x) - 1;

x(1) = 1;
x(end) = 1;
y(1) = 0;
y(end) = 0;

CamberPanels = 50000;
x_eval = flipud(linspace(0, Chord, CamberPanels/2)')';
x_eval = [x_eval, flipud(x_eval')'];
x_eval(CamberPanels/2) = [];

LEPanels = 1/20 * CamberPanels;
Body(:, 1) = x_eval';
Body(:, 2) = [flipud(interp1(x(1:LE_index), y(1:LE_index),  linspace(0, Chord, CamberPanels/2), 'pchip')')',  interp1(x(LE_index+1:end), y(LE_index+1:end), linspace(2 * Chord / CamberPanels, Chord, CamberPanels/2 - 1), 'pchip')];

Camber(:, 1) = Body(1:CamberPanels/2, 1);
for  i = 1:CamberPanels/2
    Camber(i, 2) = (Body(i, 2) + Body(end - i + 1, 2)) / 2;
end

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alphai, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);


%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);

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




figure
hold on
grid on

for alpha = alpha_vect
    alpha

    TermineNoto = zeros(NRows, 1);

    U_inf = 1;

    U_inf_x = U_inf * cos(deg2rad(alpha));
    U_inf_y = U_inf * sin(deg2rad(alpha));

    U_inf = [U_inf_x; U_inf_y];
    U_inf_normal = [-U_inf(2); U_inf(1)];
    U_inf_normal = U_inf_normal ./ norm(U_inf_normal);



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

    % faccio scorrere con due cicli for innestati tutti i pannelli per valutare
    % l'influenza che hanno su tutti gli altri pannelli

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

    plot(Centro(LE_index:end, 1), -Cp(LE_index:end))
    plot(Centro(1:LE_index, 1), -Cp(1:LE_index))
end

