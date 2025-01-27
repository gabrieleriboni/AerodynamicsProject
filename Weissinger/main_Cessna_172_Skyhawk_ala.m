clear
close all
clc

%% Inizializzazione variabili

U_Inf_Mag = 1;
alpha_vect = linspace(-5, 10, 20);

rho = 1.225;

config.NCorpi = 1;

% config.[...](1) deve essere l'ala (per calcolo corretto di CL e CD del
% velivolo completo, scalati solo sull'area dell'ala)

config.RootChord = [1.625];
config.DihedralAngle = [1.73]; % [°]
config.SweepAngle = [0]; % [°]
config.TaperRatio = [0.672];
config.AspectRatio = [7.32];
config.Span = [10.92];
config.LEPosition_X = [0];
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];

config.RotationAngle_X = [0];
config.RotationAngle_Y = [3.5];
config.RotationAngle_Z = [0];

% Opzioni di discretizzazione
config.SemiSpanwiseDiscr = [10];
config.ChordwiseDiscr = [10];

%% Calcoli preliminari

% Definisco la semiapertura, la superficie e le posizioni di radice e estremità
config.SemiSpan = config.Span./2;
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
config.TipChord = config.RootChord .* config.TaperRatio;
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Inizializzazione parametri prestazionali
CL_vect = cell(config.NCorpi, 1);
CL_3D = cell(config.NCorpi, 1);
CD_3D = cell(config.NCorpi, 1);
CD_vect = cell(config.NCorpi, 1);
CD_3D_tot = [];
CL_3D_tot = [];
visual = 1;     % Per ottenere la visualizzazione solo per un angolo di incidenza

%% Creazione della geometria

% Costruisco le celle contenenti le informazioni dei pannelli
ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);

% Creo la struttura della pannellizzazione
for iCorpo = 1:config.NCorpi
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end

%% Creazione matrici

% Inizializzo la matrice A
NPanelsTot = 2 * config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
U_j = cell(config.NCorpi, config.NCorpi);

% Determino le componenti di A per ogni pannello
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    % Ciclo per tutti i pannelli lungo la corda
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Ciclo per tutti i pannelli lungo l'apertura
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            rowIndex = rowIndex + 1;
            columnIndex = 0;
            % Calcolo la normale e il punto di collocazione
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;

            % Determino le componenti di velocità indotta per ogni pannello
            for jCorpo = 1:config.NCorpi
                U_ind{iCorpo}(ChordPanel_i, SpanPanel_i, :) = zeros(1, 3);
                % Ciclo per tutti i pannelli lungo la corda
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Ciclo per tutti i pannelli lungo l'apertura
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        columnIndex = columnIndex + 1;

                        % Induzione dovuta alla presenza dei vortici

                        % Vortici semi-infiniti
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Vortici finiti
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Influenza del secondo vortice semi-infinito
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        U_j{iCorpo, jCorpo}(ChordPanel_i, SpanPanel_i, ChordPanel_j, SpanPanel_j, :) = U;

                        % Calcolo la matrice A
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                    end
                end
            end
        end
    end
end

%% Costruzione del termine noto

% Determino le componenti della velocità U_inf al variare dei valori di
% alpha

for alpha = alpha_vect
    U_Inf = [cosd(alpha) 0 sind(alpha)] .* U_Inf_Mag;
    TermineNoto = zeros(NPanelsTot, 1);
    rowIndex = 0;
    % Determino le componenti del termine noto per ogni pannello
    for iCorpo = 1:config.NCorpi
        % Ciclo per tutti i pannelli lungo la corda
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            % Ciclo per tutti i pannelli lungo l'apertura
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                rowIndex = rowIndex + 1;
                NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            end
        end
    end

    %% Risoluzione del sistema lineare

    % Risolvo il sistema lineare Ax = b
    Solution = linsolve(matriceA, TermineNoto);

    % La soluzione del sistema lineare (Gamma) deve essere determinata per ogni
    % pannello
    Gamma = cell(config.NCorpi, 1);
    rowIndex = 0;
    for iCorpo = 1:config.NCorpi
        Gamma{iCorpo} = zeros(config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2);
        % Ciclo per tutti i pannelli lungo la corda
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            % Ciclo per tutti i pannelli lungo l'apertura
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                rowIndex = rowIndex + 1;
                Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
            end
        end
    end

    %% Visualizzazione
    if alpha > alpha_vect(floor(length(alpha_vect) / 2)) && visual == 1
        visual = 0;
        for iCorpo = 1:config.NCorpi
            X = zeros(config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo));
            Y = zeros(config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo));
            Z = zeros(config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo));
            for i = 1:config.ChordwiseDiscr(iCorpo)
                for j = config.SemiSpanwiseDiscr(iCorpo)+1:2*config.SemiSpanwiseDiscr(iCorpo)
                    X(i, j - config.SemiSpanwiseDiscr(iCorpo)) = internalMesh{iCorpo}{i, j}.LERoot(1);
                    Y(i, j - config.SemiSpanwiseDiscr(iCorpo)) = -internalMesh{iCorpo}{i, j}.LERoot(2);
                    Z(i, j - config.SemiSpanwiseDiscr(iCorpo)) = internalMesh{iCorpo}{i, j}.LERoot(3);
                    x(j-config.ChordwiseDiscr(iCorpo)) = internalMesh{iCorpo}{end, j}.TERoot(1);
                    y_te(j - config.SemiSpanwiseDiscr(iCorpo)) = -internalMesh{iCorpo}{end, j}.TERoot(2);
                end
                x_tip(i) = internalMesh{iCorpo}{i, end}.LEtip(1);
                x_tip(i+1) = internalMesh{iCorpo}{i, end}.TEtip(1);
                y(i) = -internalMesh{iCorpo}{i, end}.LEtip(2);
                z(i) = internalMesh{iCorpo}{i, end}.LEtip(3);
            end
            y_te =[y_te, -internalMesh{iCorpo}{end, end}.TEtip(2)];
            X = [X; x];
            X = [X, x_tip'];
            X_corpo{iCorpo}=X;
            Y = [Y, y'];
            Y = [Y; Y(end, :)];
            Y_corpo{iCorpo}=Y;
            Z = [Z, z'];
            Z = [Z; Z(end, :)];
            Z_corpo{iCorpo}=Z;
 
            X_add = flipud(X_corpo{iCorpo}')';
            X_corpo{iCorpo} = [X_add(:,1:config.SemiSpanwiseDiscr(iCorpo)), X_corpo{iCorpo}];
            Y_add = flipud(Y_corpo{iCorpo}')';
            Y_corpo{iCorpo} = [- Y_add(:,1:config.SemiSpanwiseDiscr(iCorpo)), Y_corpo{iCorpo}];
            Z_add = flipud(Z_corpo{iCorpo}')';
            Z_corpo{iCorpo} = [Z_add(:,1:config.SemiSpanwiseDiscr(iCorpo)), Z_corpo{iCorpo}];
        end

        %% Distribuzione di circolazione
        
        figure(40)
        Gamma_plot = Gamma;
        for iCorpo = 1:config.NCorpi
            Gamma_plot{iCorpo} =  [Gamma_plot{iCorpo}; Gamma_plot{iCorpo}(end,:)];
            Gamma_plot{iCorpo} =  [Gamma_plot{iCorpo}, Gamma_plot{iCorpo}(:,end)];

            mesh(X_corpo{iCorpo}, Y_corpo{iCorpo}, Z_corpo{iCorpo}, Gamma_plot{iCorpo},'FaceColor', 'interp','EdgeColor', 'k','FaceAlpha', 0.8)

        end
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colormap('jet')
        colorbar
        hold on
    end

    %% Calcolo di Portanza e Coefficiente di Portanza 2D e 3D

    % inizializzo le celle contenenti Lift, Cl e i segmenti dei pannelli
    L_2D = cell(config.NCorpi, 1);
    Cl_2D = cell(config.NCorpi, 1);
    c = cell(config.NCorpi, 1);

    % Effettuo il calcolo per ogni pannello
    for iCorpo = 1:config.NCorpi

        % Determino il segmento longitudinale costituente il singolo pannello
        for j = 0:config.SemiSpanwiseDiscr(iCorpo)-1
            c{iCorpo}(j+1) = config.RootChord(iCorpo) * (config.TaperRatio(iCorpo) - 1) * j / config.SemiSpanwiseDiscr(iCorpo) + config.RootChord(iCorpo);
        end
        c{iCorpo} = flipud(c{iCorpo}');
        c{iCorpo}(config.SemiSpanwiseDiscr(iCorpo)+1:2*config.SemiSpanwiseDiscr(iCorpo)) = flipud(c{iCorpo});

        % Calcolo Lift e C_l 2D
        for j = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            L_2D{iCorpo}(j) = rho * U_Inf_Mag * sum(Gamma{iCorpo}(:, j)) * cosd(config.DihedralAngle(iCorpo));
            Cl_2D{iCorpo}(j) = L_2D{iCorpo}(j) / (0.5 * rho * U_Inf_Mag^2 * c{iCorpo}(j));
        end
    end

    % Calcolo Lift e C_L 3D
    L_3D = cell(config.NCorpi, 1);
    for iCorpo = 1:config.NCorpi
        L_3D{iCorpo} = sum(L_2D{iCorpo}) * config.Span(iCorpo)/(2*config.SemiSpanwiseDiscr(iCorpo));
    end


    %% Calcolo della Resistenza Indotta e del Coefficiente di Resistenza Indotta 2D e 3D

    % inizializzo le celle contenenti Drag e alpha indotti
    D_2D = cell(config.NCorpi, 1);
    L_2D_panel = cell(config.NCorpi, 1);
    D_2D_panel = cell(config.NCorpi, 1);
    alpha_ind = cell(config.NCorpi, 1);
    U_ind = cell(config.NCorpi, 1);

    % Effettuo il calcolo per ogni pannello
    for iCorpo = 1:config.NCorpi
        U_ind{iCorpo} = zeros(config.ChordwiseDiscr(iCorpo), 2*config.SemiSpanwiseDiscr(iCorpo), 3);
        alpha_ind{iCorpo} = zeros(config.ChordwiseDiscr(iCorpo), 2*config.SemiSpanwiseDiscr(iCorpo));
        % Calcolo gli alpha indotti per i pannelli lungo l'apertura e
        % determino il drag indotto 2D
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            % Ciclo per tutti i pannelli lungo l'apertura
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                % Determino le componenti di velocità indotta per ogni pannello
                for jCorpo = 1:config.NCorpi
                    % Ciclo per tutti i pannelli lungo la corda
                    for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                        % Ciclo per tutti i pannelli lungo l'apertura
                        for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                            U_ind{iCorpo}(ChordPanel_i, SpanPanel_i, 3) = U_ind{iCorpo}(ChordPanel_i, SpanPanel_i, 3) + U_j{iCorpo, jCorpo}(ChordPanel_i, SpanPanel_i, ChordPanel_j, SpanPanel_j, 3) * Gamma{jCorpo}(ChordPanel_j, SpanPanel_j);
                        end
                    end
                end
                alpha_ind{iCorpo}(ChordPanel_i, SpanPanel_i) = alpha_ind{iCorpo}(ChordPanel_i, SpanPanel_i) + atan2(U_ind{iCorpo}(ChordPanel_i, SpanPanel_i, 3), U_Inf_Mag);
                L_2D_panel{iCorpo}(ChordPanel_i, SpanPanel_i) = rho * U_Inf_Mag * Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) * cosd(config.DihedralAngle(iCorpo));
                D_2D_panel{iCorpo}(ChordPanel_i, SpanPanel_i) = L_2D_panel{iCorpo}(ChordPanel_i, SpanPanel_i) * sin(alpha_ind{iCorpo}(ChordPanel_i, SpanPanel_i));
            end
        end
        for j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
            D_2D{iCorpo}(j) = sum(D_2D_panel{iCorpo}(:, j));
        end
    end


    % Calcolo il Drag indotto 3D
    D_3D = cell(config.NCorpi, 1);
    D_3D_tot = 0;
    L_3D_tot = 0;

    for iCorpo = 1:config.NCorpi
        D_3D{iCorpo} = abs(sum(D_2D{iCorpo})) * config.Span(iCorpo)/(2*config.SemiSpanwiseDiscr(iCorpo));

        % Calcolo CL e CD 3D per ogni corpo
        CL_3D{iCorpo} = L_3D{iCorpo} / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
        CD_3D{iCorpo} = D_3D{iCorpo} / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));

        CL_vect{iCorpo} = [CL_vect{iCorpo}; CL_3D{iCorpo}];
        CD_vect{iCorpo} = [CD_vect{iCorpo}; CD_3D{iCorpo}];

        % Calcolo Lift e Drag del velivolo completo
        D_3D_tot = D_3D_tot + D_3D{iCorpo};
        L_3D_tot = L_3D_tot + L_3D{iCorpo};
    end
    % Calcolo CL e CD del velivolo completo
    CD_3D_tot = [CD_3D_tot; D_3D_tot / (0.5 * rho * U_Inf_Mag^2 * config.Surface(1))];
    CL_3D_tot = [CL_3D_tot; L_3D_tot / (0.5 * rho * U_Inf_Mag^2 * config.Surface(1))];
end

%% Plot Cl/alpha e Cl/Cd

% % Plot per singolo corpo
% figure
% plot(alpha_vect, CL_vect{1}, '.-', 'MarkerSize', 8)
% grid on
% xlabel("\alpha [°]")
% ylabel("C_L [-]")
% title('Ala')
% 
% figure
% plot(CD_vect{1}, CL_vect{1}, '.-', 'MarkerSize', 8)
% grid on
% xlabel("C_D_prova [-]")
% ylabel("C_L [-]")
% title('Ala')
% 
% figure
% plot(alpha_vect, CL_vect{2}, '.-', 'MarkerSize', 8)
% grid on
% xlabel("\alpha [°]")
% ylabel("C_L [-]")
% title('Coda')
% 
% figure
% plot(CD_vect{2}, CL_vect{2}, '.-', 'MarkerSize', 8)
% grid on
% xlabel("C_D_prova [-]")
% ylabel("C_L [-]")
% title('Coda')

CL_tot = 0;
CL_vect_tot = zeros(length(alpha_vect), 1);
CD_vect_tot = zeros(length(alpha_vect), 1);

% CL/alpha per singolo corpo
for iCorpo = 1: config.NCorpi
    CL_alpha{iCorpo} = CL_3D{iCorpo} / alpha;
    fprintf("\nCL_alpha_corpo_" + iCorpo + " = " + CL_alpha{iCorpo} + "\n")
    CL_tot = CL_tot + CL_3D{iCorpo};
    CL_vect_tot = CL_vect_tot + CL_vect{iCorpo};
    CD_vect_tot = CD_vect_tot + CD_vect{iCorpo};
end

CL_alpha_velivolo = CL_tot / alpha

figure
subplot(1, 2, 1)
plot(alpha_vect, CL_3D_tot, '.-', 'MarkerSize', 8)
grid on
xlabel("\alpha [°]")
ylabel("C_L [-]")

subplot(1, 2, 2)
plot(CD_3D_tot, CL_3D_tot, '.-', 'MarkerSize', 8)
grid on
xlabel("C_D [-]")
ylabel("C_L [-]")
