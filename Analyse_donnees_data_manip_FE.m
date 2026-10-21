clc;
close all;
clear; %% supprime les data dans le workspace
chemin=cd; %% détermine le chemin où il y a le programme et les dossiers contenant les data
scrz = get(0,'ScreenSize');

corPR = [0.999961920164244 0.999847648748724 0.999657089962139 0.999390083930246 0.999046406361756 0.998625768078623 0.998127814408720 0.997552124438326 0.996898210121220 0.996165515240539 0.995353414218895 0.994461210771550 0.993488136396673 0.992433348695962 0.991295929518023 0.990074882916012 0.988769132910076 0.987377521044045 0.985898803724714 0.984331649330771 0.982674635077077 0.980926243618512 0.979084859375943 0.977148764565057 0.975116134906829 0.972985034996125 0.970753413302550 0.968419096774858 0.965979785017225 0.963433044002285 0.960776299281993 0.958006828653133 0.955121754229454 0.952118033867032 0.948992451883283 0.945741609003176 0.942361911458281 0.938849559155375 0.935200532821074 0.931410580017323 0.927475199909165 0.923389626650814 0.919148811238316 0.914747401656544 0.910179721124455 0.905439744214833 0.900521070592372 0.895416896076030 0.890119980687013 0.884622613291115 0.878916572381763 0.872993082475940 0.866842765506441 0.860455586487473 0.853820792602204 0.846926844705154 0.839761340042596 0.832310924761613 0.824561194491769 0.816496580927726 0.808100221897082 0.799353811839174 0.790237428913031 0.780729334048886 0.770805736093520 0.760440515686412 0.749604898516988 0.738267065977503 0.726391687688982 0.713939355571815 0.700865892514451 0.687121499445025 0.672649691475946 0.657385954804083 0.641256028056763 0.624173669593653 0.606037707139425 0.586728062719372 0.566100276468179 0.543977765132136 0.520140541313840 0.494308170534155 0.466112868580578 0.435054659798161 0.400421252587635 0.361131020719170 0.315382705803031 0.259702118227043 0.185219029975141 2.10734242554470e-08]; %% tableau de correspondance entre PR (de 1 à 90°) et excentricité (de presque 1 à 0)
[bb,aa] = butter(2,2*3/200) ;

cd([chemin '\Boucles non traitees']); %% change le chemin pour aller dans le dossier des data

dossiers=dir('*.txt'); %% détermine tous les fichiers où il y a le mot boucle dedans

for i =8: size(dossiers,1)
    cd([chemin '\Boucles non traitees']); %% retourne dans le chemin pour aller dans le dossier des data (aucun intérêt pour la première itération mais c'est important pour les autres)
    filename=[dossiers(i,1).folder '\' dossiers(i,1).name]; %% détermine le nom des fichiers
    delimiter = '\t';
    startRow = 12;
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
    fclose(fileID);
    %% Attention, il y a des levers qui font planter l'analyse / Faire une analyse par parties tracées pour la prochaine réunion
    
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    
    for col=[1,2,3,4,5,6,7,8,9,10,11]
        % Converts text in the input cell array to numbers. Replaced non-numeric
        % text with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1)
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\.]*)+[\,]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\.]*)*[\,]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData(row), regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if numbers.contains('.')
                    thousandsRegExp = '^\d+?(\.\d{3})*\,{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'))
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric text to numbers.
                if ~invalidThousandsSeparator
                    numbers = strrep(numbers, '.', '');
                    numbers = strrep(numbers, ',', '.');
                    numbers = textscan(char(numbers), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch
                raw{row, col} = rawData{row};
            end
        end
    end
    
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    data_total = cell2mat(raw);
    
    data_tot = [data_total(:,1:4) data_total(:,6:end)];
    for p=2:size(data_tot,1)
        data_tot(p,6)=data_tot(p-1,6)+data_tot(p,5);
    end
    
    %%% correction d'échelle
    data_tot(:,3)=data_tot(:,3).*0.628;
    data_tot(:,4)=data_tot(:,4).*1.608;
    %%%
    
    
    data=data_tot(data_tot(:,7)>=0,:); %% On ne garde que les parties tracées ! Attention pour le calcul de frequ qui risque d'être anormal
    data=data(11:end,:);
    
    clearvars  filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R startRow formatSpec fileID dataArray ans;%% les 9 lignes précédentes sont là pour télécharger le premier fichier segmenté
    cd(chemin); %% retourne sur le chemin où il y a les programmes d'analyse matlab
    %             figure;subplot(1,2,1);
    %             plot(data(:,6),data(:,3),'r');
    %             hold on
    %             plot(data(:,6),data(:,4));
    % %     coeffdroiteY = polyfit(data(:,6),data(:,3),2);
    % %     data(:,11)=(data(:,6).^2).*coeffdroiteY(1)+data(:,6).*coeffdroiteY(2)+coeffdroiteY(3);
    
    coeffdroiteY = polyfit(data(:,6),data(:,3),3); %% calcule le coefficient directeur de la régression linéaire de tes oscillations sur l'axe x (c'est le 'a' de y = ax +b), et l'ordonnée à l'origine (le 'b' de y = ax + b)
    %     data(:,11)=data(:,6).^3.*coeffdroiteY(1)+data(:,6).^2.*coeffdroiteY(2)+data(:,6).*coeffdroiteY(3)+ coeffdroiteY(4); %% construit la droite de régression dans ton fichier segmenté (colonne 11)
    
    
    
    data(:,12)=data(:,3);%% -data(:,11); %% on soustrait la droite de régression linéaire des oscillations sur l'axe x
    %         subplot(1,2,2)
    %     plot(data(:,6),data(:,3),'r');
    %     hold on
    %     plot(data(:,6),data(:,11),'g')
    %     plot(data(:,6),data(:,12),'g')
    %     plot(data(:,6),data(:,4));
    %     figure
    %     plot(data(:,12),data(:,4))
    %     set(gca,'PlotBoxAspectRatio',[1 1 1]),set(gca,'DataAspectRatio',[1 1 1]) ;
    
    datafilt(:,1)=filtfilt(bb,aa,data(:,4)); %% filtre ("lisse") ta courbe pour déterminer des extrema locaux corrects
    
    if str2num(dossiers(i,1).name(end-10))== 1
        [pks, loc]  = findpeaks(-datafilt); %% calcule les maxima locaux des oscillations sur x pour faire l'analyse cycle par cycle et pour calculer la fréquence des cycles (1/période)
    elseif str2num(dossiers(i,1).name(end-10))== 2
        [pks, loc]  = findpeaks(-datafilt);
    elseif str2num(dossiers(i,1).name(end-10))== 3
        [pks, loc]  = findpeaks(datafilt);
    elseif str2num(dossiers(i,1).name(end-10))== 4
        [pks, loc]  = findpeaks(datafilt);
    end
    
    
    %      plot(data(loc,6), data(loc,12),'or');
    
    figure('Color',[1 1 1],...
        'MenuBar','none',...
        'Toolbar','none', ...
        'Position',[scrz],...
        'Visible','on');
    figue = subplot(2,2,1)
    plot(data(:,3),data(:,4))
    set(gca,'PlotBoxAspectRatio',[1 1 1]),set(gca,'DataAspectRatio',[1 1 1]) ;
    
    subplot(2,2,2)
    plot(data(:,6),data(:,3),'r');
    hold on
%     plot(data(:,6),data(:,11),'g')
%     plot(data(:,6),data(:,12),'g')
    plot(data(:,6),data(:,4));
    
    vitx=diff(data(:,3))./0.005;
    vity=diff(data(:,4))./0.005;
    vit=sqrt(vitx.^2+vity.^2);
    vitmoy = mean(vit);
    vitsd = std(vit);
    dd=1;
    for j= 1: size(loc,1)-1
        subplot(2,size(loc,1)-1,(j+size(loc,1)-1))
        plot(data(loc(j):loc(j+1),12),data(loc(j):loc(j+1),4),'or'); %% affiche la boucle (de chaque cycle)
        set(gca,'PlotBoxAspectRatio',[1 1 1]),set(gca,'DataAspectRatio',[1 1 1]) ;
        
        
        vitc_x=diff(data(loc(j):loc(j+1),3))./0.005;
        vitc_y=diff(data(loc(j):loc(j+1),4))./0.005;
        vit_cycle=sqrt(vitc_x.^2+vitc_y.^2);
        vit_moy_cycle(j) = mean(vit_cycle);
        vit_sd_cycle(j) = std(vit_cycle);
        ratio_cycle(j) = vit_moy_cycle(j)./vit_sd_cycle(j);
        
        %changer le ration à 0,5 pour les boucles qui beugent
        if ratio_cycle(j)> 0.5 && dd<9 %%% variable seuil pour considérer un cycle comme correct
            period2(j) = (data(loc(j+1),6)-data(loc(j),6))./1000; %% calcule la période de ton cycle en s.
            frequ2(j) = 1./ period2(j);
            
            nbre_points_cycle=loc(j+1)-loc(j);
            
            %             portion_ok=round(nbre_points_cycle./5);
            %             x = [data(loc(j)+portion_ok:loc(j+1)-portion_ok,12)';data(loc(j)+portion_ok:loc(j+1)-portion_ok,4)']; %% partitionne les données pour l'anay=lyse cycle par cycle
            %
            
            x = [data(loc(j):loc(j+1),12)';data(loc(j):loc(j+1),4)']; %% partitionne les données pour l'anay=lyse cycle par cycle
            
            if str2num(dossiers(i,1).name(end-10))== 1
                maxx = find(x(2,:)==max(x(2,:)));
            elseif str2num(dossiers(i,1).name(end-10))== 2
                maxx = find(x(2,:)==max(x(2,:)));
            elseif str2num(dossiers(i,1).name(end-10))== 3
                maxx = find(x(2,:)==min(x(2,:)));
            elseif str2num(dossiers(i,1).name(end-10))== 4
                maxx = find(x(2,:)==min(x(2,:)));
            end
                  
            rep_max = maxx(1);
            
            for k = 1:rep_max-5
                compax(k,:)=x(1,k)-x(1,rep_max+1:end);
                compay(k,:)=x(2,k)-x(2,rep_max+1:end);
                compa(k,:)=sqrt(compax(k,:).^2+compay(k,:).^2);
                RepPointCompa{k}=min(compa(k,:));
                RepCompa(k)=RepPointCompa{k}(1);
            end
            intersection1=find(RepCompa==min(RepCompa));
            %             hold on
            %             plot(x(1,intersection1),x(2,intersection1),'ob');
            kk=1;
            for k = rep_max+1:size(x,2)
                compax2(k,:)=x(1,k)-x(1,1:rep_max-5);
                compay2(k,:)=x(2,k)-x(2,1:rep_max-5);
                compa2(k,:)=sqrt(compax2(k,:).^2+compay2(k,:).^2);
                RepPointCompa2{k}=min(compa2(k,:));
                RepCompa2(kk)=RepPointCompa2{k}(1);
                kk=kk+1;
            end
            intersection2=find(RepCompa2==min(RepCompa2));
            hold on
            if intersection2+rep_max+1 <= size(x,2)
                plot(x(1,intersection1:intersection2+rep_max+1),x(2,intersection1:intersection2+rep_max+1),'og');
                x2=x(:,intersection1:intersection2+rep_max+1);
            else
                plot(x(1,intersection1:end),x(2,intersection1:end),'og');
                x2=x(:,intersection1:end);
            end
            
            if size(x2,2)>40
                [zg, ag, bg, alphag] = fitellipse(x2,'linear'); %% calcule les valeurs mathématiques de l'ellipse la plus proche de ta boucle(sans la translation)
                
                
                
                
                plotellipse(zg, ag, bg, alphag, 'k') %% affiche l'ellipse sur ta figure
                orientation2(j)=alphag./pi.*180; %% convertit l'orientation de l'ellipe (calculée en radian) en degré
                
                if ag>bg %% cette boucle if et la suivante normalisent l'excentricité et l'orientation de l'ellipse
                    c(j)=sqrt((ag.^2)-(bg.^2));
                    e(j)=c(j)./ag;
                else
                    c(j)=sqrt((bg.^2)-(ag.^2));
                    e(j)=c(j)./bg;
                    orientation2(j)=orientation2(j)+90;
                end
                
                if orientation2(j) <0
                    orientation2(j) = orientation2(j) + 180;
                end
                if orientation2(j) >180
                    orientation2(j) = orientation2(j) - 180;
                end
                %         e(j) = sqrt(1-((ag.^2)./(bg.^2)));
                hold off
                % plot(zg,'or')
                
                rep = find(corPR<=e(j));
                PR2(j) = rep(1);
                dd=dd+1;
            else
                PR2(j) = NaN;
                orientation2(j) = NaN;
                period2(j) = NaN;
                frequ2(j) = NaN;
            end
        else
            PR2(j) = NaN;
            orientation2(j) = NaN;
            period2(j) = NaN;
            frequ2(j) = NaN;
        end
        
        clear RepCompa2 RepCompa x2 RepPointCompa2 RepPointCompa compax compay compax2 compay2 compa compa2 intersection 1 intersection2 x zg ag bg alphag rep vitc_x vitc_y vit_cycle %% supprime les données qui n'ont plus d'intérêt
    end
    
    %%% calcul des moyennes et ET / ATTENTION : on doit enlever les valeurs
    %%% abérantes liées aux levers...
    orientation=orientation2(orientation2(:)>=0);
    PR=PR2(PR2(:)>=0);
    period=period2(period2(:)>=0);
    frequ=frequ2(frequ2(:)>=0);
    
    essais_reussis=find(orientation2(:)>=0);
    essais_rates=essais_reussis(8)-8;
    
    orientation3=2.*orientation(1:8).*pi./180;
    orient=circMean(orientation3);
    
    orient_moy=orient(1)./2;%% calcule l'orientation moyenne de toutes les boucles de l'essai /// ATTENTION PRENDRE LA MOYENNE CIRCULAIRE !!!!!!!!!!!
    if orient_moy<0
        orient_moy=orient_moy+180;
    end
    orient_std=orient(2); %% calcule la déviation standard (variabilité) des orientations de toutes les boucles de l'essai
    nombre_boucles = size(loc,1)-1;
    period_moy = mean(period(1:8)); %% calcule la période moyenne de toutes tes boucles
    period_std = std(period(1:8)); %% calcule la déviation standard (variabilité) de ta période
    frequ_moy = mean(frequ(1:8)); %% calcule ta freq moyenne
    frequ_std = std(frequ(1:8)); %% calcule la déviation standard de ta fréquence
    PR_moy = mean(PR(1:8)); %% calcule la PR moyenne de ta série de boucles
    PR_std = std(PR(1:8)); %% calcule la DS de la PR
    %             exc_moy=mean(e);
    %             exc_std=std(e);
    
    %% Construction de ton tableau final pour les analyses stat
    donnees{1,1} = 'nom_fichier';
    donnees{1,2} = 'orientation moy';
    donnees{1,3} = 'orientation sd';
    donnees{1,4} = 'nombre_boucles';
    donnees{1,5} = 'frequ_moy';
    donnees{1,6} = 'frequ sd';
    donnees{1,7} = 'PR_moy';
    donnees{1,8} = 'PR sd';
    donnees{1,9} = 'Cte de translation';
    donnees{1,10} = 'Boucles ratées';
    
    donnees{i+1,1} = dossiers(i,1).name;
    donnees{i+1,2} = orient_moy;
    donnees{i+1,3} = orient_std;
    donnees{i+1,4} = nombre_boucles;
    donnees{i+1,5} = frequ_moy;
    donnees{i+1,6} = frequ_std;
    donnees{i+1,7} = PR_moy;
    donnees{i+1,8} = PR_std;
    donnees{i+1,9} = coeffdroiteY(1);
    donnees{i+1,10} = essais_rates;
    i
    pause(2);
    %% KbWait;
    close all
    %%
    clear orientation3 essais_rates coeffdroiteY filename data orient_moy orient_std nombre_boucles frequ_moy frequ_std PR_moy PR_std datafilt
    
end
