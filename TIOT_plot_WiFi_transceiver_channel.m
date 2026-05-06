clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

% SNR = 0:10:40
%% Mecha
Mecha_low_threeTags_Tag1 = [0.391398913043478;0.285838405797101;0.254851811594203;0.249388405797101;0.248703623188406];
Mecha_low_threeTags_Tag2 = [0.392121739130435;0.286634782608696;0.255654347826087;0.250392028985507;0.249648550724638];
Mecha_low_threeTags_Tag3 = [0.391971739130435;0.286892391304348;0.255930072463768;0.250605434782609;0.249917391304348];

Mecha_medium_threeTags_Tag1 = [0.0523869565217391;0.0240460144927536;0.0197050724637681;0.0190550724637681;0.0189630434782609];
Mecha_medium_threeTags_Tag2 = [0.0530572463768116;0.0238699275362319;0.0194192028985507;0.0188250000000000;0.0187394927536232];
Mecha_medium_threeTags_Tag3 = [0.0534735507246377;0.0242166666666667;0.0199264492753623;0.0192949275362319;0.0192003623188406];

Mecha_high_threeTags_Tag1 = [0.0200231884057971;0.00180760869565217;0.000215942028985507;4.20289855072464e-05;3.69565217391304e-05];
Mecha_high_threeTags_Tag2 = [0.0208271739130435;0.00196050724637681;0.000176811594202899;4.34782608695652e-05;3.22463768115942e-05];
Mecha_high_threeTags_Tag3 = [0.0209876811594203;0.00165181159420290;0.000216666666666667;9.78260869565217e-05;7.60869565217391e-05];

%% MultiRider
MultiRider_low_threeTags_Tag1 = [0.550000724637681;0.547101811594203;0.546622101449275;0.546559057971015;0.546541666666667];
MultiRider_low_threeTags_Tag2 = [0.553911956521739;0.551211594202899;0.550804347826087;0.550711231884058;0.550737318840580];
MultiRider_low_threeTags_Tag3 = [0.556535869565217;0.553845289855072;0.553720289855073;0.553701811594203;0.553721376811594];

MultiRider_medium_threeTags_Tag1 = [0.0455554347826087;0.0424250000000000;0.0421692028985507;0.0421735507246377;0.0421800724637681];
MultiRider_medium_threeTags_Tag2 = [0.0665307971014493;0.0490478260869565;0.0470858695652174;0.0469057971014493;0.0468634057971015];
MultiRider_medium_threeTags_Tag3 = [0.0853586956521739;0.0546583333333333;0.0508043478260870;0.0504431159420290;0.0503608695652174];

MultiRider_high_threeTags_Tag1 = [0.000878985507246377;3.62318840579710e-06;0;0;0];
MultiRider_high_threeTags_Tag2 = [0.0184989130434783;0.00187934782608696;0.000138405797101449;1.63043478260870e-05;7.24637681159420e-06];
MultiRider_high_threeTags_Tag3 = [0.0364971014492754;0.00430543478260870;0.000495652173913043;6.73913043478261e-05;1.37681159420290e-05];


%% FreeCollision
FreeCollision_low_threeTags_Tag1 = [0.574456313131313;0.580109595959596;0.591351010101010;0.595955555555556;0.596337121212121];
FreeCollision_low_threeTags_Tag2 = [0.587772474747475;0.593767929292929;0.604078787878788;0.607892424242424;0.608497474747475];
FreeCollision_low_threeTags_Tag3 = [0.611811363636364;0.618719949494950;0.628788636363636;0.631305050505051;0.631094949494950];

FreeCollision_medium_threeTags_Tag1 = [0.151832323232323;0.0841684343434343;0.0753813131313131;0.0763861111111111;0.0767404040404040];
FreeCollision_medium_threeTags_Tag2 = [0.164201767676768;0.0900651515151515;0.0809636363636364;0.0815146464646465;0.0820025252525253];
FreeCollision_medium_threeTags_Tag3 = [0.210568434343434;0.109180808080808;0.0948909090909091;0.0948669191919192;0.0951068181818182];

FreeCollision_high_threeTags_Tag1 = [0.0963633838383838;0.0157636363636364;0.00245631313131313;0.00101994949494950;0.000748484848484849];
FreeCollision_high_threeTags_Tag2 = [0.106756060606061;0.0174494949494950;0.00273055555555556;0.00108005050505051;0.000877777777777778];
FreeCollision_high_threeTags_Tag3 = [0.149704797979798;0.0246075757575758;0.00387222222222222;0.00160202020202020;0.00146616161616162];


%% ParaFi
ParaFi_low_threeTags_Tag1 = [0.216625000000000;0.0418619565217391;0.00329094202898551;0.000336231884057971;9.27536231884058e-05];
ParaFi_low_threeTags_Tag2 = [0.213634782608696;0.0418605072463768;0.00323913043478261;0.000314855072463768;0.000113768115942029];
ParaFi_low_threeTags_Tag3 = [0.216109782608696;0.0421231884057971;0.00359239130434783;0.000396376811594203;0.000111594202898551];

ParaFi_medium_threeTags_Tag1 = [0.216625000000000;0.0418619565217391;0.00329094202898551;0.000336231884057971;9.27536231884058e-05];
ParaFi_medium_threeTags_Tag2 = [0.213634782608696;0.0418605072463768;0.00323913043478261;0.000314855072463768;0.000113768115942029];
ParaFi_medium_threeTags_Tag3 = [0.216109782608696;0.0421231884057971;0.00359239130434783;0.000396376811594203;0.000111594202898551];

ParaFi_high_threeTags_Tag1 = [0.216625000000000;0.0418619565217391;0.00329094202898551;0.000336231884057971;9.27536231884058e-05];
ParaFi_high_threeTags_Tag2 = [0.213634782608696;0.0418605072463768;0.00323913043478261;0.000314855072463768;0.000113768115942029];
ParaFi_high_threeTags_Tag3 = [0.216109782608696;0.0421231884057971;0.00359239130434783;0.000396376811594203;0.000111594202898551];


%% ConcurScatter
ConcurScatter_low_threeTags_Tag1 = [0.566626535087719;0.566404385964912;0.566494298245614;0.566420614035088;0.566389035087719];
ConcurScatter_low_threeTags_Tag2 = [0.566723026315790;0.566480921052632;0.566540350877193;0.566566008771930;0.566602631578947];
ConcurScatter_low_threeTags_Tag3 = [0.567165789473684;0.567071710526316;0.567176754385965;0.567206359649123;0.567146929824561];

ConcurScatter_medium_threeTags_Tag1 = [0.499958552631579;0.500224780701754;0.500332017543860;0.500260307017544;0.500267982456140];
ConcurScatter_medium_threeTags_Tag2 = [0.499564254385965;0.499484429824561;0.499318859649123;0.499377192982456;0.499270833333333];
ConcurScatter_medium_threeTags_Tag3 = [0.500127192982456;0.500050877192983;0.499849780701754;0.499919298245614;0.499888377192982];

ConcurScatter_high_threeTags_Tag1 = [0.499963596491228;0.499937280701754;0.500076535087719;0.500024122807018;0.500062719298246];
ConcurScatter_high_threeTags_Tag2 = [0.499396491228070;0.499273026315789;0.499099561403509;0.499155701754386;0.499025438596491];
ConcurScatter_high_threeTags_Tag3 = [0.500268859649123;0.500263596491228;0.500121491228070;0.500249561403509;0.500203508771930];

%% BER
%%% low
BER_Mecha_low = (Mecha_low_threeTags_Tag1 + Mecha_low_threeTags_Tag2 + Mecha_low_threeTags_Tag3)./3;
BER_FreeCollision_low = (FreeCollision_low_threeTags_Tag1 + FreeCollision_low_threeTags_Tag2 + FreeCollision_low_threeTags_Tag3)./3;
BER_ParaFi_low = (ParaFi_low_threeTags_Tag1 + ParaFi_low_threeTags_Tag2 + ParaFi_low_threeTags_Tag3)./3;
BER_MultiRider_low = (MultiRider_low_threeTags_Tag1 + MultiRider_low_threeTags_Tag2 + MultiRider_low_threeTags_Tag3)./3;
BER_ConcurScatter_low = (ConcurScatter_low_threeTags_Tag1 + ConcurScatter_low_threeTags_Tag2 + ConcurScatter_low_threeTags_Tag3)./3;

%%% medium 
BER_Mecha_medium = (Mecha_medium_threeTags_Tag1 + Mecha_medium_threeTags_Tag2 + Mecha_medium_threeTags_Tag3)./3;
BER_FreeCollision_medium = (FreeCollision_medium_threeTags_Tag1 + FreeCollision_medium_threeTags_Tag2 + FreeCollision_medium_threeTags_Tag3)./3;
BER_ParaFi_medium = (ParaFi_medium_threeTags_Tag1 + ParaFi_medium_threeTags_Tag2 + ParaFi_medium_threeTags_Tag3)./3;
BER_MultiRider_medium = (MultiRider_medium_threeTags_Tag1 + MultiRider_medium_threeTags_Tag2 + MultiRider_medium_threeTags_Tag3)./3;
BER_ConcurScatter_medium = (ConcurScatter_medium_threeTags_Tag1 + ConcurScatter_medium_threeTags_Tag2 + ConcurScatter_medium_threeTags_Tag3)./3;

%%% high
BER_Mecha_high = (Mecha_high_threeTags_Tag1 + Mecha_high_threeTags_Tag2 + Mecha_high_threeTags_Tag3)./3;
BER_FreeCollision_high = (FreeCollision_high_threeTags_Tag1 + FreeCollision_high_threeTags_Tag2 + FreeCollision_high_threeTags_Tag3)./3;
BER_ParaFi_high = (ParaFi_high_threeTags_Tag1 + ParaFi_high_threeTags_Tag2 + ParaFi_high_threeTags_Tag3)./3;
BER_MultiRider_high = (MultiRider_high_threeTags_Tag1 + MultiRider_high_threeTags_Tag2 + MultiRider_high_threeTags_Tag3)./3;
BER_ConcurScatter_high = (ConcurScatter_high_threeTags_Tag1 + ConcurScatter_high_threeTags_Tag2 + ConcurScatter_high_threeTags_Tag3)./3;

%% Throughput
k = 250; % Kbps
%%% low
TP_Mecha_low = (1-Mecha_low_threeTags_Tag1).*k + (1-Mecha_low_threeTags_Tag2).*k + (1-Mecha_low_threeTags_Tag3).*k;
TP_MultiRider_low = (1-MultiRider_low_threeTags_Tag1).*k + (1-MultiRider_low_threeTags_Tag2).*k + (1-MultiRider_low_threeTags_Tag3).*k;
TP_FreeCollision_low = (1-FreeCollision_low_threeTags_Tag1).*k + (1-FreeCollision_low_threeTags_Tag2).*k + (1-FreeCollision_low_threeTags_Tag3).*k;
TP_ParaFi_low = (1-ParaFi_low_threeTags_Tag1).*k + (1-ParaFi_low_threeTags_Tag2).*k + (1-ParaFi_low_threeTags_Tag3).*k;
TP_ConcurScatter_low = (1-ConcurScatter_low_threeTags_Tag1).*k + (1-ConcurScatter_low_threeTags_Tag2).*k + (1-ConcurScatter_low_threeTags_Tag3).*k;

%%% medium
TP_Mecha_medium = (1-Mecha_medium_threeTags_Tag1).*k + (1-Mecha_medium_threeTags_Tag2).*k + (1-Mecha_medium_threeTags_Tag3).*k;
TP_MultiRider_medium = (1-MultiRider_medium_threeTags_Tag1).*k + (1-MultiRider_medium_threeTags_Tag2).*k + (1-MultiRider_medium_threeTags_Tag3).*k;
TP_FreeCollision_medium = (1-FreeCollision_medium_threeTags_Tag1).*k + (1-FreeCollision_medium_threeTags_Tag2).*k + (1-FreeCollision_medium_threeTags_Tag3).*k;
TP_ParaFi_medium = (1-ParaFi_medium_threeTags_Tag1).*k + (1-ParaFi_medium_threeTags_Tag2).*k + (1-ParaFi_medium_threeTags_Tag3).*k;
TP_ConcurScatter_medium = (1-ConcurScatter_medium_threeTags_Tag1).*k + (1-ConcurScatter_medium_threeTags_Tag2).*k + (1-ConcurScatter_medium_threeTags_Tag3).*k;

%%% high
TP_Mecha_high = (1-Mecha_high_threeTags_Tag1).*k + (1-Mecha_high_threeTags_Tag2).*k + (1-Mecha_high_threeTags_Tag3).*k;
TP_MultiRider_high = (1-MultiRider_high_threeTags_Tag1).*k + (1-MultiRider_high_threeTags_Tag2).*k + (1-MultiRider_high_threeTags_Tag3).*k;
TP_FreeCollision_high = (1-FreeCollision_high_threeTags_Tag1).*k + (1-FreeCollision_high_threeTags_Tag2).*k + (1-FreeCollision_high_threeTags_Tag3).*k;
TP_ParaFi_high = (1-ParaFi_high_threeTags_Tag1).*k + (1-ParaFi_high_threeTags_Tag2).*k + (1-ParaFi_high_threeTags_Tag3).*k;
TP_ConcurScatter_high = (1-ConcurScatter_high_threeTags_Tag1).*k + (1-ConcurScatter_high_threeTags_Tag2).*k + (1-ConcurScatter_high_threeTags_Tag3).*k;

%% plot BER
vals_BER_low = [BER_FreeCollision_low,BER_Mecha_low,BER_MultiRider_low,BER_ParaFi_low,BER_ConcurScatter_low];
vals_BER_medium = [BER_FreeCollision_medium,BER_Mecha_medium,BER_MultiRider_medium,BER_ParaFi_medium,BER_ConcurScatter_medium];
vals_BER_high = [BER_FreeCollision_high,BER_Mecha_high,BER_MultiRider_high,BER_ParaFi_high,BER_ConcurScatter_high];


%% plot Throughput
vals_TP_low = [TP_FreeCollision_low,TP_Mecha_low,TP_MultiRider_low,TP_ParaFi_low,TP_ConcurScatter_low];
vals_TP_medium = [TP_FreeCollision_medium,TP_Mecha_medium,TP_MultiRider_medium,TP_ParaFi_medium,TP_ConcurScatter_medium];
vals_TP_high = [TP_FreeCollision_high,TP_Mecha_high,TP_MultiRider_high,TP_ParaFi_high,TP_ConcurScatter_high];

x = 0:10:40;

%%%** plot BER **%%%
figure(1)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_BER_low);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
ax1.YLim = [1e-5 1e0];
ax1.YTick = [1e-6,1e-4,1e-2,1e0];
ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

figure(2)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_BER_medium);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
ax1.YLim = [1e-5 1e0];
ax1.YTick = [1e-6,1e-4,1e-2,1e0];
ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

figure(3)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_BER_high);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
ax1.YLim = [1e-5 1e0];
ax1.YTick = [1e-6,1e-4,1e-2,1e0];
ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 



figure(4)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_TP_low);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
% ax1.YLim = [1e-5 1e0];
% ax1.YTick = [1e-6,1e-4,1e-2,1e0];
% ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

figure(5)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_TP_medium);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
% ax1.YLim = [1e-5 1e0];
% ax1.YTick = [1e-6,1e-4,1e-2,1e0];
% ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

figure(6)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_TP_high);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255;
b1(5).FaceColor = [132/255 165/255 157/255];
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.XLim = [5 35];
% ax1.YLim = [1e-5 1e0];
% ax1.YTick = [1e-6,1e-4,1e-2,1e0];
% ax1.YScale = 'log';
% ax1.XTickLabel = {'low quality','medium quality','high quality'};
grid on
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

