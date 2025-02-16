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
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(signal_centre, 'r');
title('Signal Centré (Moyenne supprimée)');
xlabel('Temps');
ylabel('Amplitude');
grid on;

signal_centre_positif = max(signal_centre, 0);
subplot(3,1,3);
plot(signal_centre_positif, 'g');
title('Signal Centré et Positif');
xlabel('Temps');
ylabel('Amplitude');
grid on;

% Sauvegarde uniquement si le fichier n'existe pas
filepath = fullfile(output_dir, 'signaux.png');
save_if_not_exists(f, filepath);

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



% Comparaison avec et sans interpolation
nb_classes_list = [5, 10];
colors = ['r', 'g', 'b', 'm', 'c'];

for k = 1:length(nb_classes_list)
    nb_classes = nb_classes_list(k);
    classes = linspace(0, 1, nb_classes+1);
    moyenne_phase = zeros(1, nb_classes);
    moyenne_phase_interp = zeros(1, nb_classes);
    ecart_type_phase = zeros(1, nb_classes);
    ecart_type_phase_interp = zeros(1, nb_classes);
    valeurs_par_classe = cell(1, nb_classes);
    valeurs_par_classe_interp = cell(1, nb_classes);
    
    for i = 1:length(tmax)-1
        indices_periode = find(t >= tmax(i) & t < tmax(i+1));
        t_periode = t(indices_periode);
        signal_periode = signal_centre(indices_periode);
        t_normalise = (t_periode - tmax(i)) / (tmax(i+1) - tmax(i));
        
        indices_periode_interp = find(t >= tmax_interp(i) & t < tmax_interp(i+1));
        t_periode_interp = t(indices_periode_interp);
        signal_periode_interp = signal_centre(indices_periode_interp);
        t_normalise_interp = (t_periode_interp - tmax_interp(i)) / (tmax_interp(i+1) - tmax_interp(i));
        
        for j = 1:nb_classes
            indices_classe = find(t_normalise >= classes(j) & t_normalise < classes(j+1));
            valeurs_par_classe{j} = [valeurs_par_classe{j}; signal_periode(indices_classe)];
            
            indices_classe_interp = find(t_normalise_interp >= classes(j) & t_normalise_interp < classes(j+1));
            valeurs_par_classe_interp{j} = [valeurs_par_classe_interp{j}; signal_periode_interp(indices_classe_interp)];
        end
    end
    
    for j = 1:nb_classes
        moyenne_phase(j) = mean(valeurs_par_classe{j});
        ecart_type_phase(j) = std(valeurs_par_classe{j});
        moyenne_phase_interp(j) = mean(valeurs_par_classe_interp{j});
        ecart_type_phase_interp(j) = std(valeurs_par_classe_interp{j});
    end
    
    % --- Figure pour Moyenne de Phase ---
    f = figure('Visible', 'off');
    plot(linspace(0,1,nb_classes), moyenne_phase, '-o', 'Color', colors(k), 'LineWidth', 1.5);
    title(sprintf('Moyenne de Phase - %d Classes', nb_classes));
    xlabel('Temps Normalisé (t^*)'); ylabel('Amplitude Moyenne');
    grid on;
    filepath = fullfile(output_dir, sprintf('Moyenne_Phase_%d.png', nb_classes));
    save_if_not_exists(f, filepath);
    close(gcf);

    % --- Figure pour Écart-Type de la Phase ---
    f = figure('Visible', 'off');
    plot(linspace(0,1,nb_classes), ecart_type_phase, '-o', 'Color', colors(k), 'LineWidth', 1.5);
    title(sprintf('Écart-Type de la Phase - %d Classes', nb_classes));
    xlabel('Temps Normalisé (t^*)'); ylabel('Écart-Type');
    grid on;
    filepath = fullfile(output_dir, sprintf('Ecart_Type_Phase_%d.png', nb_classes));
    save_if_not_exists(f, filepath);
    close(gcf);

    % --- Figure pour Moyenne de Phase avec Interpolation ---
    f = figure('Visible', 'off');
    plot(linspace(0,1,nb_classes), moyenne_phase_interp, '-o', 'Color', colors(k), 'LineWidth', 1.5);
    title(sprintf('Moyenne de Phase (Interp) - %d Classes', nb_classes));
    xlabel('Temps Normalisé (t^*)'); ylabel('Amplitude Moyenne');
    grid on;
    filepath = fullfile(output_dir, sprintf('Moyenne_Phase_Interp_%d.png', nb_classes));
    save_if_not_exists(f, filepath);
    close(gcf);

    % --- Figure pour Écart-Type de la Phase avec Interpolation ---
    f = figure('Visible', 'off');
    plot(linspace(0,1,nb_classes), ecart_type_phase_interp, '-o', 'Color', colors(k), 'LineWidth', 1.5);
    title(sprintf('Écart-Type de la Phase (Interp) - %d Classes', nb_classes));
    xlabel('Temps Normalisé (t^*)'); ylabel('Écart-Type');
    grid on;
    filepath = fullfile(output_dir, sprintf('Ecart_Type_Phase_Interp_%d.png', nb_classes));
    save_if_not_exists(f, filepath);
    close(gcf);
end

% Calcul des périodes successives (différences entre tmax_interp successifs)
T = diff(tmax_interp); 

% Estimation des paramètres de la distribution normale
mu_T = mean(T);
sigma_T = std(T);


% Fonction pour enregistrer une figure seulement si elle n'existe pas
function save_if_not_exists(fig, filepath)
    if ~exist(filepath, 'file')
        saveas(fig, filepath);
    else
        fprintf('Le fichier %s existe déjà, non sauvegardé.\n', filepath);
    end
    close(fig);
end

% Exemple pour l'histogramme de densité de probabilité
figure('Visible', 'off');
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
save_if_not_exists(gcf, filepath);