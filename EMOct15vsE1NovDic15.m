%% TRAERSE TODOS LOS DATOS
clear
load 'C:\Users\GENETIC SYSTEMS 3\Documents\Abraham\L6Shared_Abraham\phd\Experimentos\E1NovDic15\AnalysisE1NovDic15_13Ene16.mat'
bgdataAllNov = bgdataClean;
timeODNov = timeOD;
gwrateNov = gwrate;
survivalNov = survival;
decayRateNov = decayRate;
load 'C:\Users\GENETIC SYSTEMS 3\Documents\Abraham\L6Shared_Abraham\phd\Experimentos\EMOct15 longevidad de KO en AN\Analisis EMOct15 compu L6\Abraham\Abraham\EMOCT15\9Dic'


%% PLOTS DE CURVAS DE DECAIMIENTO SEPARADAS POR COLORES PARA EMOCT15 (basado en Analysis EMOCT15 2.0 10Nov)
survival=getxvec(bgdataClean,survival,pls,od,odTh,hash,keys,strains)
%%
strains={'S288c(Mosaic/Lab)','NaN','S288c(Mosaic/Lab)','YIIc17_E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM981(Europe/Clinical)','Y12(Asia/Sake)','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L_1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L_1528(Europe/Ferment)','YJM978(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
i=1;
llaves=fieldnames(hashDos);
for plato=pls 
    figure(200+plato)
    %subplot(4,4,i)
    for j=2:length(llaves)
        subplot(3,4,j-1)
    plot(survival(plato).t(:,:), survival(plato).s(:,:),'y-' )
    hold on
    
    vector=hashDos(plato).RefOnly; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'g.-' )
    
    vector=hashDos(plato).Mutonly; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'c.-' )
    
    vector=hashDos(plato).CompRef; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'r.-' )
    
    vector=hashDos(plato).(cell2mat(llaves(j))); %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'k.-' )
    title(llaves(j))
    %pause
    end
    
    title(strains(plato))
    i=i+1;
    
end

%%%Sigue compararestas gráficas con las que se usan para calcular las s
%%
%%
%% post analysis: RoC and s3,s4 plate by plate (96 wells per figure.)
pl=1
longevityPostNoSWAP_EMOCT15(bgdataCleanS,8,12,pl,4,4);

%%  CDF PLOTS!!!
%CDFPlots y barras con error de las Little_S de cada Mutante por plato
%Genera matrixS, que tiene todas las little_S
for plt=[1, 2 3:18]
    y=[]; %media de las réplicas
    err=[]; %error estandar de las replicas
    N=[]; % Número de réplicas tomadas en cuenta  %mediana=[];
    vec=[];
    keys = fieldnames(hashDos);
    for i=1:length(keys) %Para cada mutante
        llave=keys(i);
        figure(300+plt)
        subplot(3,5,i)
        hold off
        y = [y nanmean(bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave))))];
        %mediana = [mediana nanmedian(bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave))))];
        N = [N length(find(~isnan(bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave))))))];
        err = [ nanstd(bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave))))/sqrt(N(i)), err];
        if N(i)
            cdfplotAA2(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)),'r')
            hold on
            cdfplotAA(bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave)))-nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef))  )
            ylim([0 1])
            plot([-.06 -.06], [0 1] )
            plot([.06 .06], [0 1] )
            xlim([-.16 .16])
            title(strcat(mat2str(N(i)),llave))
        end
    end
    ENES(plt,:)=N;
    matrixS(plt,:)=y; %Aquí guardamos todas las little_s
    normalizedS(plt,:)=matrixS(plt,:)-matrixS(plt,4);%Generala la normalizedS a partir de matrixS Ss
    matrixErr(plt,:)=err; %Errores std para todas   
subplot(3,5,i+1)
title(strains(plt))
    %%%%%% ESTA PARTE SON LAS BARRAS VERDES CON ERROR ROJO%%%%%
% figure(plt+100)
% hold off
% h=bar(normalizedS(plt,:), 'g');
% set(gca,'xticklabel',keys)
% hold on
% vec(1:length(keys))=-0.05;
% errorbar(1:length(keys),normalizedS(plt,:), err,'.r')
% text(1:length(keys),vec, mat2cell(N))
% title('media')
% ylim([-.2 .2])
end


%% NECESITO QUE EN UNA SOLA PANTALLA PUEDA VER BIEN
%cdfPLOT DE S COMPREF Y MUTANTE
%curva decaimiento Compref y mutante
%Plots de las RoC y S de Compref y mutante
%plots de RoC
%Plots de RFP y YFP

strains={'S288c(Mosaic/Lab)','NaN','S288c(Mosaic/Lab)','YIIc17_E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM981(Europe/Clinical)','Y12(Asia/Sake)','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L_1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L_1528(Europe/Ferment)','YJM981(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
llaves=fieldnames(hashDos);
for mutante=8%MSN2%1:length(llaves);%BMH1
for plato=1:18%pls 
    plt=plato;
    figure()
    subplot(1,2,1)
%estas son todas las curvas de decaimiento  
    for j=mutante%2:length(llaves)
    
    plot(survival(plato).t(:,:), survival(plato).s(:,:),'y-' )
    hold on
    
    vector=hashDos(plato).RefOnly; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'g.-' )
    
    vector=hashDos(plato).Mutonly; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'c.-' )

    vector=hashDos(plato).AllRef; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'m*-' )
    vector=hashDos(plato).CompRef; %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'r*-' )


    vector=hashDos(plato).(cell2mat(llaves(j))); %[ 31:12:72 36:12:72];%25:12:72
    plot(survival(plato).t(:,vector), survival(plato).s(:,vector),'k*-' )
    title(strains(plato))
    ylabel('% de sobrevivencia')
    %pause
    end
% esta parte hace las RoC y s Rojo CompRef, Negro Mutante
    subplot(1,2,2)

for w = hashDos(plt).AllRef
    plot(bgdataCleanS(plt).tDays,bgdataCleanS(plt).RoC(:,w),'mo')
    hold on
    y2=bgdataCleanS(plt).s(1,w)+bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:), y2,'m-');
    text(bgdataCleanS(plt).tDays(end),y2(end),mat2str(w) )
end
for w = hashDos(plt).CompRef
    plot(bgdataCleanS(plt).tDays,bgdataCleanS(plt).RoC(:,w),'ro')
    hold on
    y2=bgdataCleanS(plt).s(1,w)+bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:), y2,'r*-');
    text(bgdataCleanS(plt).tDays(end),y2(end),mat2str(w) )
end
for w = hashDos(plt).RefOnly
    plot(bgdataCleanS(plt).tDays,bgdataCleanS(plt).RoC(:,w),'go')
    hold on
    y2=bgdataCleanS(plt).s(1,w)+bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:), y2,'g-');
    text(bgdataCleanS(plt).tDays(end),y2(end),mat2str(w) )
end
for w = hashDos(plt).Mutonly
    plot(bgdataCleanS(plt).tDays,bgdataCleanS(plt).RoC(:,w),'co')
    hold on
    y2=bgdataCleanS(plt).s(1,w)+bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:), y2,'c-');
    text(bgdataCleanS(plt).tDays(end),y2(end),mat2str(w) )
end
for w = hashDos(plt).(cell2mat(llaves(mutante)))
    plot(bgdataCleanS(plt).tDays,bgdataCleanS(plt).RoC(:,w),'ko')
    hold on
    y2=bgdataCleanS(plt).s(1,w)+bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:), y2,'k-');
    text(bgdataCleanS(plt).tDays(end),y2(end),mat2str(w) )
end
    title(llaves(j))

end
end

%%

%longevityPostNoSWAP_EMOCT15(bgdata,nr,nc,plates,ylimTh, ylimTh2)
platos=16;
for i=1:length(platos)
    plt=platos(i);
    counter=1;
    figure(500+plt)
    for w=1:96
        mysubplot(8,12,counter,-.020,-.020); %como vas a saber cuantos pozos son para graficar????
        plot(log10(bgdataCleanS(plt).rfp(:,w)), 'r-')
        axis fill
        counter=counter+1;
        ylim([2.5 4.5])
        grid on
    end;  
    counter=1;
    figure(600+plt)
    for w=1:96
        mysubplot(8,12,counter,-.020,-.020); %como vas a saber cuantos pozos son para graficar????
        plot(log10(bgdataCleanS(plt).yfp(:,w)), 'y-')
        hold on
        plot(log10(bgdataCleanS(plt).yfp(:,w)), 'k.')
        ylim([2.5 4.5])
        axis fill
        counter=counter+1;
        grid on
    end;  
	counter=1;
    figure(700+plt)
    for w=1:96
        mysubplot(8,12,counter,-.020,-.020); %como vas a saber cuantos pozos son para graficar????
        plot(bgdataCleanS(plt).od(:,w), 'k-')
        ylim([0 1])
        axis fill
        counter=counter+1;
    end;  
end;

%%
%% S2Ene16 figura con el cambio en la curva de crecimiento conforme pasan los días

colours = ['k';'r';'c';'m';'g';'b'];
syms = ['o-';'*-';'o-';'o-';'o-';'o-';'+-';'+-';'+-';'+-';'+-'];
lc = length(colours);
ls = length(syms);
for pl=1:2;
for well=14:24%16:12:83%22%14:24;
intervalo=intervalsODAA(bgdataClean,-.3,pl);
figure(well);
subplot(1,2,1)
for i=1:length(intervalo)-1
    hold on
    colsym = [colours(mod(i-1,lc)+1), syms(mod(i-1,ls)+1)];
    plot( bgdataClean(pl).t(intervalo(i)+1:intervalo(i+1))-bgdataClean(pl).t(intervalo(i)+1), bgdataClean(pl).OD(intervalo(i)+1:intervalo(i+1),well), colsym )
   %pause
end
legend('Día1','Día2','Día3','Día5','Día7','Día9','Día11','Día13','Día15')
xlim([0 18])
for i=1:length(intervalo)-1
    colsym = [colours(mod(i-1,lc)+1)];
    plot( bgdataClean(pl).t(intervalo(i)+1:intervalo(i+1))-bgdataClean(pl).t(intervalo(i)+1), bgdataClean(pl).OD(intervalo(i)+1:intervalo(i+1),well), colsym )
end
subplot(1,2,2)

end
end

%% ES PARA PONER LAS s  Y LOS RoC  DE UNA SOLA MUTANTE 


strains={'S288c(Mosaic/Lab)','NaN','S288c(Mosaic/Lab)','YIIc17_E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM981(Europe/Clinical)','Y12(Asia/Sake)','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L_1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L_1528(Europe/Ferment)','YJM978(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
lasBuenas=[3 4 6 8 12 16]
llaves=fieldnames(hashDos);
contador=1;
for mutante=5:length(llaves)%MSN2%1:length(llaves);%BMH1
    contador=1
for plato=lasBuenas%1:18%pls 
    plt=plato;
    figure(mutante)
    subplot(2,3,contador)
    contador=contador+1;
for w = hashDos(plt).CompRef
    hold on
    y2=bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:)./24, y2,'k-');
end
for w = hashDos(plt).(cell2mat(llaves(mutante)))
    hold on
    y2=bgdataCleanS(plt).s(2,w)*bgdataCleanS(plt).tDays(:)./24;
    plot(bgdataCleanS(plt).tDays(:)./24, y2,'r-');
end
    %title(llaves(j))
    ylim([-1 1])
    title(strains(plato))
end
%subplot(2,3,6)
%title(llaves(mutante))
end

%%
%% bar lits por cepa
figure(533); con=0; clf;
%lits=bigS;
strains={'S288c(Mosaic/Lab)','NaN','S288c(Mosaic/Lab)','YIIc17-E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM975(Europe/Clinical)','NaN','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L-1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L-1528(Europe/Ferment)','YJM978(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
%lits=lits*-1;
bigS=normalizedS(:,4:end);
errs=matrixErr(:,4:end);

PLBs=lasBuenas;%[3 4 5 6 7 8 10 11 12 13 15 16];
for i=PLBs%pls%[1 3:8 10:18]
con=con+1;
%figure(con)
subplot(2,3,con)
bar(bigS(i,:),'g')
hold on
errorbar( 1:length(bigS(i,:)),bigS(i,:), errs(i,:),'.r')
vec2(1:length(keys(4:end)))=-.1;
%text(vec2,1:length(keys(4:end)),mat2Cell(ENES(i,4:end)) )
%set(gca,'yticklabel',keys(4:end) )
set(gca,'yticklabel', {''})
%xlim([-.12 .12])
title(strains(i))
ylim([-.07 .07])
xlim([0 11])
%pause
end
%% bar lits por gen KO
figure(433); con=0; clf;
%lits=bigS;
strains={'S288c(Mosaic/Lab)',  'NaN',     'S288c(Mosaic/Lab)','YIIc17-E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM975(Europe/Clinical)','NaN','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L-1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L-1528(Europe/Ferment)','YJM978(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
%lits=lits*-1;
bigS=normalizedS(:,4:end);
errs=matrixErr(:,4:end);

PLBs=lasBuenas;%[3 4 5 6 7 8 10 11 12 13 15 16];
for i=5:length(llaves)%pls%[1 3:8 10:18]
con=con+1;
%figure(con)
subplot(3,3,con)
bar(bigS(PLBs,i-3),'g')
hold on
errorbar( 1:length(bigS(PLBs,i-3)),bigS(PLBs,i-3), errs(PLBs,i-3),'.r')
%vec2(1:length(keys(4:end)))=-.1;
%text(vec2,1:length(keys(4:end)),mat2Cell(ENES(i,4:end)) )
%set(gca,'yticklabel',keys(4:end) )
%set(gca,'yticklabel', {''})
%xlim([-.12 .12])
title(llaves(i))
ylim([-.07 .07])
xlim([0 7])
%pause
end
%%
for i=5:length(llaves)%pls%[1 3:8 10:18]
con=con+1;
%figure(con)
subplot(3,3,con)
hold on
for j=1:10
    plot(j, bigS(i,:), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    hold on
    plot(j, nanmedian(bigS(PLBs(j),i)), 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
end
plot([0, 7], [errs(i,1) errs(i,1)])
plot([0, 7], [-errs(i,1) -errs(i,1)])

set(gca, 'xtick',1:6, 'xticklabel', 1:6)
%errorbar( 1:length(bigS(PLBs,i-3)),bigS(PLBs,i-3), errs(PLBs,i-3),'.r')
%vec2(1:length(keys(4:end)))=-.1;
%text(vec2,1:length(keys(4:end)),mat2Cell(ENES(i,4:end)) )
%set(gca,'yticklabel',keys(4:end) )
%set(gca,'yticklabel', {''})
%xlim([-.12 .12])
title(llaves(i))
ylim([-.07 .07])
xlim([0 7])
%pause
end

  
end

%%
%% bar lits por gen KO
figure(433); con=0; clf;
%lits=bigS;
strains={'S288c(Mosaic/Lab)',  'NaN',     'S288c(Mosaic/Lab)','YIIc17-E5(Mosaic/Ferment)','BC187(Europe/Ferment)','Y55(African/Lab)','273614N(Mosaic/Clinical)','YJM975(Europe/Clinical)','NaN','YPS128(American/Wild)','DBVPG6044(African/Ferment)','YPS606(American/Wild)','L-1374(Europe/Ferment)','SK1(African/Lab)','UOWPS87-2421(Mosaic/Wild)','L-1528(Europe/Ferment)','YJM978(Europe/Clinical)','UWOPS05-227.2(Malasia/Wild)'};
%lits=lits*-1;
bigS=normalizedS(:,4:end);
errs=matrixErr(:,4:end);

PLBs=lasBuenas;%[3 4 5 6 7 8 10 11 12 13 15 16];
for i=5:length(llaves)%pls%[1 3:8 10:18]
con=con+1;
%figure(con)
subplot(3,3,con)
 plot(j, bigS(PLBs(j),i), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
hold on
%errorbar( 1:length(bigS(PLBs,i-3)),bigS(PLBs,i-3), errs(PLBs,i-3),'.r')
%vec2(1:length(keys(4:end)))=-.1;
%text(vec2,1:length(keys(4:end)),mat2Cell(ENES(i,4:end)) )
%set(gca,'yticklabel',keys(4:end) )
%set(gca,'yticklabel', {''})
%xlim([-.12 .12])
title(llaves(i))
ylim([-.07 .07])
xlim([0 7])
%pause
end

%%

con=1;

for plt=PLBs
    y=[]; %media de las réplicas
    err=[]; %error estandar de las replicas
    N=[]; % Número de réplicas tomadas en cuenta  %mediana=[];
    vec=[];
    keys = fieldnames(hashDos);
    
    subplot(2,3,con)
    con=con+1;
    for i=6:length(keys) %Para cada mutante

        llave=keys(i);
        y=bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave)));
        if(length(y)>1)
            plot(i-4,y-matrixS(plt,4),'ok', 'MarkerFaceColor','k')
            hold on
            plot(i-4,nanmedian(y-matrixS(plt,4)),'og', 'MarkerFaceColor','g')
        end
    
    end
  %  plot([0, 7], [errs(plt,4) errs(plt,4)])
  %  plot([0, 7], [-errs(plt,4) -errs(plt,4)])

    %xlim([1:6])
end

%%
%% Hace barras, dispersion
clear tobar
for mutante=5%:length(keys)
    llave=keys(mutante);
    figure(mutante)
    i=1;
    clf
    for plt=[3 4 6 8 12 16]
        fondo=nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef));
        %plot(i, nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo), 'go','MarkerSize',10,'MarkerFaceColor','g')
        hold on
        plot(i, bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo, 'ko','MarkerSize',8)
        temp=length(bgdataCleanS(plt).s(2,hashDos(plt).CompRef));
        tobar(i)=nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo);
        tobarerr(i)=std(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo);
        i=i+1;
        y=bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave)))-fondo;
        if length(y)>1
        %plot(i, nanmean(y), 'go', 'MarkerFaceColor','g','MarkerSize',12)
        plot(i, y, 'ro', 'MarkerSize',8)
        tobar(i)=nanmean(y);
        tobarerr(i)=std(y);
                [h p]=ttest2(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo,y) 
        else
            tobar(i)=0;
        tobarerr(i)=0;
        end
        if p<.05
            plot(i-1,max(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo)+.01,'b*','MarkerSize',8)
        end
        i=i+1;
    end
    ylabel('Relative survival')
    title(llave)
    ylim([-0.07 0.07])
    xlim([0 i])
    b2 = bar(tobar,'FaceColor',[0,1,0]);    %bar(tobar)
    errorbar(tobar,tobarerr,'b.')
    set(get(b2,'Children'),'FaceAlpha',0.5)
end
%% Subplot con Todas las barras con su error y estad[istica
clear tobar tobarerr
figure()
for mutante=5:length(keys)
    llave=keys(mutante);
    subplot(3,3,mutante-4)
    hold on
    i=1;
    for plt=[3 4 6 8 12 16]
        fondo=nanmean(bgdataCleanS(plt).s(2,hashDos(plt).CompRef));
        temp=length(bgdataCleanS(plt).s(2,hashDos(plt).CompRef));
        y=bgdataCleanS(plt).s(2,hashDos(plt).(cell2mat(llave)))-fondo;
        if length(y)>1
        tobar(i)=nanmean(y);
        tobarerr(i)=std(y)/sqrt(length(y));
        [h p]=ttest2(bgdataCleanS(plt).s(2,hashDos(plt).CompRef)-fondo,y) 
           if p<.05
            plot(i,.06,'k*','MarkerSize',8)
           end
           if p<.01
            plot(i,.05,'k*','MarkerSize',8)
           end
        else
            tobar(i)=0;
            tobarerr(i)=0;
        end
        i=i+1;
    end
%    ylabel('Relative survival')
    title(llave)
    xlim([0 6.5])
    b2 = bar(tobar,'FaceColor',[0,1,0]);    %bar(tobar)
    hold on
    errorbar(tobar,tobarerr,'b.')
    set(gca,'xtick', 1:6)
    set(get(b2,'Children'),'FaceAlpha',0.5)
    ylim([-0.07 0.07])
end

%% Hace las curvas de supervivencia para LGAC
close all
i=1;
for plt=[3 4 6 8 12 16 14 15 17 18]
    subplot(2,5,i)
    hold on
    i=i+1;
   % figure()
    plot(survival(plt).t(:,hashDos(plt).AllRef),survival(plt).s(:,hashDos(plt).AllRef),'b' )
    hold on
    plot(survival(plt).t(:,hashDos(plt).AllRef),survival(plt).s(:,hashDos(plt).AllRef),'k.' )
    
    ylim([0 100])
    title(plt)
end
%%
save 20Ene16
% \\10.10.20.7\lab6_shared\PERSONAL FOLDERS\Abraham\phd\EMOct15_vs_E1NovDic15
%% Comparar la variación que observo con la que encontro Garay en las 2bles
%% mutantes

cd '\\10.10.20.7\lab6_shared\PERSONAL FOLDERS\Abraham\phd\EMOct15_vs_E1NovDic15\'
load 20Ene16
cd '\\Labseis\lab6_shared\PERSONAL FOLDERS\Abraham\phd\Epistasis Longevity\GENERAL DATA\GENERAL DATA\cluster_analysis\'
edit AnalysisEpistasis25June.m
load 'little_s (2).mat'
cd '\\10.10.20.7\lab6_shared\PERSONAL FOLDERS\Abraham\phd\EMOct15_vs_E1NovDic15\'

%%
a=[];b=[];c=[];d=[];
for i = 5:length(keys)
    for j=1:length(namesx)
        if strcmp(keys(i),namesx(j))
            a=[a i]
            keys(i)
            b=[b j]
            namesx(j)
        end
    end
	for k=1:length(namesy)
        if strcmp(keys(i),namesy(k))
            c=[c i]
            keys(i)
            d=[d k]
            namesy(k)
        end
    end
end


xs=x(:,1); %deja de ser matrix x y se vuelve vector x
ys=y(1,:); %deja de ser matrix y y se vuelve vector y
Xerr
Yerr
%%
con=1;
figure()
for i=d(1)%solo el primero porque los dem[as est[anen la otra direccion, o sea que en la otra direccion hay mas
    subplot(2,3,con)
    %figure(i)
    con=con+1;
    [vec I]=sort(xy(:,i)-xs)
    bar([ vec ],'m')
    hold on
%          plot([ xs(I) ],'c')%La otra mutante sencilla
%          plot([ xy(I,i) ],'k')%Doble mutante
    plot([1 length(xs)+2],[ys(i) ys(i)],'b','linewidth', 1 )
    title(namesy(i))
    xlim([0 75])
    ylim([-1.5 1.5])
    ylabel('s. doble mutante (k) - "fondo genético" (c)')
end
con=1;
for i=[ 1:3:71 ]%b%[12 14 23 24 65 76]%
    %figure(i)
    subplot(4,6,con)
    con=con+1;
    [vec I]=sort(xy(i,:)-ys )
    bar(vec,'m')
     hold on
%          plot([ys(I)],'c')
%          plot([xy(i,I)],'k')
    plot([1 length(ys)+2],[xs(i) xs(i)],'b','linewidth', 1)
    title(namesx(i))
    xlim([0 145])
        ylim([-1.5 1.5])

end






