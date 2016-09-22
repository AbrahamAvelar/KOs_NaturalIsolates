%This will find the lines where there is a change in days of experiment

function intervals = intervalsODAA(bgdata,odTh,pl)

diffOD = diff(bgdata(pl).od,1,1); %create matrix where large negative numbers indicate a change in the day of measurement
meandiffOD = mean(diffOD');%%No estaba en el código original
%Sirve para que un outlier no sea el culpable de la missidentificación de
%días, y en cambio que sea el promedio. así un outlier solo afecta 1/96
minOD = (min(meandiffOD',[],2)); %sort data of minimums (2nd dim measure points not wells!) in order to get all the lines with lower od than threshold od
ind = find(minOD<odTh);
intervals = [0;ind;(size(bgdata(pl).od,1)-1)]; %intervals in each day of measurement
%figure();%%estas sirven para que una vez que lo veas, identifiques el odTh
%plot(mean(diffOD'))%%
end

