% Charger les données
file_path = fullfile('..', 'data', 'signal.dat.txt'); % Chemin relatif
signal = load(file_path); % Chargement du signal

% Création du dossier de sauvegarde s'il n'existe pas
output_dir = fullfile('..', 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Calculer la moyenne temporelle
signal_moy = mean(signal);

% Centrer le signal autour de zéro
signal_centre = signal - signal_moy;

% Affichage des signaux initiaux
f = figure('Visible', 'off'); % Assigner la figure à une variable
subplot(3,1,1);
plot(signal, 'b');
title('Signal Original');
xlabel('Temps');
ylabel('Débit (L/min)');
grid on;

subplot(3,1,2);
plot(signal_centre, 'r');
title('Signal Centré (Moyenne supprimée)');
xlabel('Temps');
ylabel('Débit (L/min)');
grid on;

signal_centre_positif = max(signal_centre, 0);
subplot(3,1,3);
plot(signal_centre_positif, 'g');
title('Signal Centré et Positif');
xlabel('Temps');
ylabel('Débit (L/min)');
grid on;

% Sauvegarde uniquement si le fichier n'existe pas
filepath = fullfile(output_dir, 'signaux.png');
%save_if_not_exists(f, filepath);
saveas(f, filepath);

% Identification des intervalles [t_a, t_b] associés aux valeurs positives
t = (0:length(signal_centre_positif)-1) * 0.01; % Intervalle d'échantillonnage 10ms
indices_positifs = find(signal_centre_positif > 0);
diff_indices = diff(indices_positifs);
separations = find(diff_indices > 1);
ta = [indices_positifs(1); indices_positifs(separations+1)];
tb = [indices_positifs(separations); indices_positifs(end)];
ta = t(ta);
tb = t(tb);

% Détection des maxima initiaux
tmax = [];
for i = 1:length(ta)
    indices_intervalle = find(t >= ta(i) & t <= tb(i));
    signal_intervalle = signal_centre_positif(indices_intervalle);
    [valeur_max, indice_max_local] = max(signal_intervalle);
    tmax = [tmax; t(indices_intervalle(indice_max_local))];
end

% Interpolation quadratique autour des maxima
tmax_interp = [];
for i = 1:length(ta)
    indices_intervalle = find(t >= ta(i) & t <= tb(i));
    signal_intervalle = signal_centre_positif(indices_intervalle);
    [valeur_max, idx_max] = max(signal_intervalle);
    
    if idx_max > 1 && idx_max < length(signal_intervalle)
        x = t(indices_intervalle(idx_max-1:idx_max+1));
        y = signal_intervalle(idx_max-1:idx_max+1);
        p = polyfit(x, y, 2); % Ajustement par un polynôme de degré 2
        tmax_interp = [tmax_interp; -p(2)/(2*p(1))]; % Vertex de la parabole
    else
        tmax_interp = [tmax_interp; t(indices_intervalle(idx_max))];
    end
end

% Comparaison graphique des tmax avec et sans interpolation
f_tmax_comparaison = figure('Visible', 'off');
plot(t, signal_centre_positif, 'k', 'LineWidth', 1.5); % Signal en noir
hold on;
plot(tmax, interp1(t, signal_centre_positif, tmax, 'linear'), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Sans interpolation
plot(tmax_interp, interp1(t, signal_centre_positif, tmax_interp, 'linear'), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5); % Avec interpolation
title('Comparaison des tmax détectés avec et sans interpolation');
xlabel('Temps (s)');
ylabel('Amplitude du signal');
legend('Signal centré positif', 'tmax sans interpolation', 'tmax avec interpolation', 'Location', 'best');
grid on;

% Sauvegarde du graphique
filepath = fullfile(output_dir, 'Comparaison_tmax.png');
saveas(f_tmax_comparaison, filepath);


% Comparaison avec et sans interpolation
nb_classes_list = [5, 10, 20, 50, 100, 200];
colors = ['r', 'g', 'b', 'm', 'c', 'k'];

% Initialisation des figures
f_moyenne_phase = figure('Visible', 'off');
f_ecart_type_phase = figure('Visible', 'off');
f_moyenne_phase_interp = figure('Visible', 'off');
f_ecart_type_phase_interp = figure('Visible', 'off');

% Initialisation des matrices pour stocker les résultats
moyenne_phase_all = cell(1, length(nb_classes_list));
ecart_type_phase_all = cell(1, length(nb_classes_list));
moyenne_phase_interp_all = cell(1, length(nb_classes_list));
ecart_type_phase_interp_all = cell(1, length(nb_classes_list));
centres_classes_all = cell(1, length(nb_classes_list));


for k = 1:length(nb_classes_list)
    nb_classes = nb_classes_list(k);
    classes = linspace(0, 1, nb_classes+1);
    
    % Initialisation
    moyenne_phase = zeros(1, nb_classes);
    moyenne_phase_interp = zeros(1, nb_classes);
    ecart_type_phase = zeros(1, nb_classes);
    ecart_type_phase_interp = zeros(1, nb_classes);
    valeurs_par_classe = cell(1, nb_classes);
    valeurs_par_classe_interp = cell(1, nb_classes);

    % Boucle sur les périodes
    for i = 1:length(tmax)-1
        indices_periode = find(t >= tmax(i) & t < tmax(i+1));
        t_periode = t(indices_periode);
        signal_periode = signal(indices_periode);
        t_normalise = (t_periode - tmax(i)) / (tmax(i+1) - tmax(i));

        indices_periode_interp = find(t >= tmax_interp(i) & t < tmax_interp(i+1));
        t_periode_interp = t(indices_periode_interp);
        signal_periode_interp = signal(indices_periode_interp);
        t_normalise_interp = (t_periode_interp - tmax_interp(i)) / (tmax_interp(i+1) - tmax_interp(i));

        for j = 1:nb_classes
            indices_classe = find(t_normalise >= classes(j) & t_normalise < classes(j+1));
            valeurs_par_classe{j} = [valeurs_par_classe{j}; signal_periode(indices_classe)];

            indices_classe_interp = find(t_normalise_interp >= classes(j) & t_normalise_interp < classes(j+1));
            valeurs_par_classe_interp{j} = [valeurs_par_classe_interp{j}; signal_periode_interp(indices_classe_interp)];
        end
    end

    % Calcul des statistiques pour chaque classe
    for j = 1:nb_classes
        moyenne_phase(j) = mean(valeurs_par_classe{j});
        ecart_type_phase(j) = std(valeurs_par_classe{j});
        moyenne_phase_interp(j) = mean(valeurs_par_classe_interp{j});
        ecart_type_phase_interp(j) = std(valeurs_par_classe_interp{j});
    end

    % Calcul des centres des classes
    centres_classes = (classes(1:end-1) + classes(2:end)) / 2;

    % Stocker les résultats pour la deuxième boucle
    moyenne_phase_all{k} = moyenne_phase;
    ecart_type_phase_all{k} = ecart_type_phase;
    moyenne_phase_interp_all{k} = moyenne_phase_interp;
    ecart_type_phase_interp_all{k} = ecart_type_phase_interp;
    centres_classes_all{k} = centres_classes;

    % Sauvegarde des figures individuelles (inchangé)
end


% Boucle pour créer les graphiques comparatifs
for k = 1:length(nb_classes_list)
    nb_classes = nb_classes_list(k);
    
    % Tracé des résultats sans interpolation
    figure(f_moyenne_phase);
    plot(centres_classes, moyenne_phase, '-', 'Color', colors(k), 'LineWidth', 1.5);
    hold on;

    figure(f_ecart_type_phase);
    plot(centres_classes, ecart_type_phase, '-', 'Color', colors(k), 'LineWidth', 1.5);
    hold on;

    % Tracé des résultats avec interpolation
    figure(f_moyenne_phase_interp);
    plot(centres_classes, moyenne_phase_interp, '-', 'Color', colors(k), 'LineWidth', 1.5);
    hold on;

    figure(f_ecart_type_phase_interp);
    plot(centres_classes, ecart_type_phase_interp, '-', 'Color', colors(k), 'LineWidth', 1.5);
    hold on;
end

% Ouvrir les figures AVANT la boucle pour éviter qu'elles ne se ferment
figure(f_moyenne_phase);
hold on;
figure(f_ecart_type_phase);
hold on;
figure(f_moyenne_phase_interp);
hold on;
figure(f_ecart_type_phase_interp);
hold on;

for k = 1:length(nb_classes_list)
    figure(f_moyenne_phase);
    plot(centres_classes_all{k}, moyenne_phase_all{k}, '-', 'Color', colors(k), 'LineWidth', 1.5);

    figure(f_ecart_type_phase);
    plot(centres_classes_all{k}, ecart_type_phase_all{k}, '-', 'Color', colors(k), 'LineWidth', 1.5);

    figure(f_moyenne_phase_interp);
    plot(centres_classes_all{k}, moyenne_phase_interp_all{k}, '-', 'Color', colors(k), 'LineWidth', 1.5);

    figure(f_ecart_type_phase_interp);
    plot(centres_classes_all{k}, ecart_type_phase_interp_all{k}, '-', 'Color', colors(k), 'LineWidth', 1.5);
end

% Finalisation des figures (ajout des titres, labels, légendes, etc.)
figure(f_moyenne_phase);
title('Moyenne de Phase sans Interpolation');
xlabel('Temps Normalisé (t^*)'); ylabel('Débit Moyen (L/min)');
grid on;
legend(arrayfun(@num2str, nb_classes_list, 'UniformOutput', false), 'Location', 'best');
saveas(f_moyenne_phase, fullfile(output_dir, 'Moyenne_Phase_sans_Interp.png'));

figure(f_ecart_type_phase);
title('Écart-Type sans Interpolation');
xlabel('Temps Normalisé (t^*)'); ylabel('Écart-Type');
grid on;
legend(arrayfun(@num2str, nb_classes_list, 'UniformOutput', false), 'Location', 'best');
saveas(f_ecart_type_phase, fullfile(output_dir, 'Ecart_Type_Phase_sans_Interp.png'));

figure(f_moyenne_phase_interp);
title('Moyenne de Phase avec Interpolation');
xlabel('Temps Normalisé (t^*)'); ylabel('Débit Moyen (L/min)');
grid on;
legend(arrayfun(@num2str, nb_classes_list, 'UniformOutput', false), 'Location', 'best');
saveas(f_moyenne_phase_interp, fullfile(output_dir, 'Moyenne_Phase_avec_Interp.png'));

figure(f_ecart_type_phase_interp);
title('Écart-Type avec Interpolation');
xlabel('Temps Normalisé (t^*)'); ylabel('Écart-Type');
grid on;
legend(arrayfun(@num2str, nb_classes_list, 'UniformOutput', false), 'Location', 'best');
saveas(f_ecart_type_phase_interp, fullfile(output_dir, 'Ecart_Type_Phase_avec_Interp.png'));


% Calcul des périodes successives (différences entre tmax_interp successifs)
T = diff(tmax_interp); 

% Test de Lilliefors pour la normalité (alternative à Shapiro-Wilk)
[h_lillie, p_lillie] = lillietest(T);

% Affichage des résultats
fprintf('Test de Lilliefors (normalité) : h = %d, p = %.5f\n', h_lillie, p_lillie);

% ---- Ajustement à d'autres lois ----
% Loi Normale
mu_T = mean(T);
sigma_T = std(T);
pdf_norm = @(x) normpdf(x, mu_T, sigma_T);

% Loi Log-normale
param_log = fitdist(T, 'Lognormal');
pdf_lognorm = @(x) lognpdf(x, param_log.mu, param_log.sigma);

% Loi Gamma
param_gamma = fitdist(T, 'Gamma');
pdf_gamma = @(x) gampdf(x, param_gamma.a, param_gamma.b);

% Loi de Weibull
param_weibull = fitdist(T, 'Weibull');
pdf_weibull = @(x) wblpdf(x, param_weibull.a, param_weibull.b);

% ---- Tracé des distributions ----
x_vals = linspace(min(T), max(T), 100);

figure;
histogram(T, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
plot(x_vals, pdf_norm(x_vals), 'r', 'LineWidth', 2, 'DisplayName', 'Normale');
plot(x_vals, pdf_lognorm(x_vals), 'g', 'LineWidth', 2, 'DisplayName', 'Log-normale');
plot(x_vals, pdf_gamma(x_vals), 'm', 'LineWidth', 2, 'DisplayName', 'Gamma');
plot(x_vals, pdf_weibull(x_vals), 'c', 'LineWidth', 2, 'DisplayName', 'Weibull');

title('Comparaison des distributions ajustées');
xlabel('Période T (s)');
ylabel('Densité de probabilité');
legend;
grid on;

% Sauvegarde
filepath = fullfile(output_dir, 'Comparaison_Distributions.png');
saveas(gcf, filepath);


% Fonction pour enregistrer une figure seulement si elle n'existe pas
function save_if_not_exists(fig, filepath)
    if ~exist(filepath, 'file')
        saveas(fig, filepath);
    else
        fprintf('Le fichier %s existe déjà, non sauvegardé.\n', filepath);
    end
    close(fig);
end

% Histogramme de densité de probabilité
figure;
histogram(T, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
x = linspace(min(T), max(T), 100);
y = normpdf(x, mu_T, sigma_T);
plot(x, y, 'r', 'LineWidth', 2);
title('Densité de probabilité des périodes successives');
xlabel('Période T (s)');
ylabel('Densité de probabilité');
legend({'Histogramme des périodes', 'Distribution Gaussienne ajustée'}, 'Location', 'northeast');
grid on;

% Chemin du fichier à enregistrer
filepath = fullfile(output_dir, 'Densite_Probabilite_Periods.png');
saveas(gcf, filepath);

% Comparaison des moyennes de phase et des écarts-types
figure;

% Moyenne de phase : sans interpolation
subplot(2, 1, 1);
plot(centres_classes, moyenne_phase, '-r', 'LineWidth', 1.5); % Sans interpolation
hold on;
plot(centres_classes, moyenne_phase_interp, '-b', 'LineWidth', 1.5); % Avec interpolation
title('Comparaison des Moyennes de Phase');
xlabel('Temps Normalisé (t^*)');
ylabel('Débit Moyen (L/min)');
legend('Sans Interpolation', 'Avec Interpolation');
grid on;

% Écart-type : sans interpolation
subplot(2, 1, 2);
plot(centres_classes, ecart_type_phase, '-r', 'LineWidth', 1.5); % Sans interpolation
hold on;
plot(centres_classes, ecart_type_phase_interp, '-b', 'LineWidth', 1.5); % Avec interpolation
title('Comparaison des Écarts-Types');
xlabel('Temps Normalisé (t^*)');
ylabel('Écart-Type');
legend('Sans Interpolation', 'Avec Interpolation');
grid on;

% Sauvegarde de la figure
filepath = fullfile(output_dir, 'Impact_Interpolation.png');
saveas(gcf, filepath);

% Calcul des différences entre tmax avec et sans interpolation
differences_tmax = abs(tmax - tmax_interp);

% Visualisation graphique des différences avec échelle logarithmique
figure;
semilogy(1:length(differences_tmax), differences_tmax, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
title('Différences entre tmax avec et sans interpolation');
xlabel('Index des tmax');
ylabel('|tmax_{sans interp} - tmax_{avec interp}| (s)');
grid on;

% Personnalisation des ticks pour refléter les petites valeurs
yticks = [1e-4, 1e-3, 1e-2, 1e-1]; % Ticks adaptés à la plage des petites valeurs
set(gca, 'YTick', yticks);
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.1e', x), yticks, 'UniformOutput', false));

% Sauvegarde du graphique
filepath = fullfile(output_dir, 'Differences_tmax_log.png');
saveas(gcf, filepath);