clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

% SNR = 30; packet number: 10000;

%% Mecha
Mecha_twoTags_0 = [1.60714285714286e-05,2.23214285714286e-05];
Mecha_twoTags_1 = [0.0972187500000000,0.0978517857142857,0.410667857142857];
Mecha_twoTags_2 = [0.335913392857143,0.334403571428571,0.321233928571429,0.379992857142857];

%% MultiRider
MultiRider_twoTags_0 = [0,8.92857142857143e-07];
MultiRider_twoTags_1 = [0.0891169642857143,0.0893303571428571,0.410795535714286];
MultiRider_twoTags_2 = [0.178772321428571,0.178495535714286,0.322210714285714,0.321552678571429];

%% ParaFi
ParaFi_twoTags_0 = [0.000150000000000000,0.000183035714285714];
ParaFi_twoTags_1 = [0.106267857142857,0.106650892857143,0.423707142857143];
ParaFi_twoTags_2 = [0.388066071428571,0.385753571428571,0.394010714285714,0.394200000000000];

%% FreeCollision
FreeCollision_twoTags_0 = [0.000204545454545455,0.000337121212121212];
FreeCollision_twoTags_1 = [0.000715909090909091,0.00107575757575758,0.00158863636363636];
FreeCollision_twoTags_2 = [0.0106068181818182,0.0109022727272727,0.0141356060606061,0.0210287878787879];

%% ConcurScatter
ConcurScatter_twoTags_0 = [0.499370394736842,0.500127631578947];
ConcurScatter_twoTags_1 = [0.500286842105263,0.499392763157895,0.499557894736842];
ConcurScatter_twoTags_2 = [0.498957236842105,0.499279605263158,0.499615789473684,0.499401315789474];

%% avg BER
BER_Mecha_twoTags_0 = mean(Mecha_twoTags_0);
BER_Mecha_twoTags_1 = mean(Mecha_twoTags_1);
BER_Mecha_twoTags_2 = mean(Mecha_twoTags_2);

BER_MultiRider_twoTags_0 = mean(MultiRider_twoTags_0);
BER_MultiRider_twoTags_1 = mean(MultiRider_twoTags_1);
BER_MultiRider_twoTags_2 = mean(MultiRider_twoTags_2);

BER_ParaFi_twoTags_0 = mean(ParaFi_twoTags_0);
BER_ParaFi_twoTags_1 = mean(ParaFi_twoTags_1);
BER_ParaFi_twoTags_2 = mean(ParaFi_twoTags_2);

BER_FreeCollision_twoTags_0 = mean(FreeCollision_twoTags_0);
BER_FreeCollision_twoTags_1 = mean(FreeCollision_twoTags_1);
BER_FreeCollision_twoTags_2 = mean(FreeCollision_twoTags_2);

BER_ConcurScatter_twoTags_0 = mean(ConcurScatter_twoTags_0);
BER_ConcurScatter_twoTags_1 = mean(ConcurScatter_twoTags_1);
BER_ConcurScatter_twoTags_2 = mean(ConcurScatter_twoTags_2);

BER_Mecha = [BER_Mecha_twoTags_0;BER_Mecha_twoTags_1;BER_Mecha_twoTags_2];
BER_MultiRider = [BER_MultiRider_twoTags_0;BER_MultiRider_twoTags_1;BER_MultiRider_twoTags_2];
BER_ParaFi = [BER_ParaFi_twoTags_0;BER_ParaFi_twoTags_1;BER_ParaFi_twoTags_2];
BER_FreeCollision = [BER_FreeCollision_twoTags_0;BER_FreeCollision_twoTags_1;BER_FreeCollision_twoTags_2];
BER_ConcurScatter = [BER_ConcurScatter_twoTags_0;BER_ConcurScatter_twoTags_1;BER_ConcurScatter_twoTags_2];


%% Throughput
k = 250; % Kbps
TP_Mecha_twoTags_0 = (2-sum(Mecha_twoTags_0)).*k;
TP_Mecha_twoTags_1 = (3-sum(Mecha_twoTags_1)).*k;
TP_Mecha_twoTags_2 = (4-sum(Mecha_twoTags_2)).*k;

TP_MultiRider_twoTags_0 = (2-sum(MultiRider_twoTags_0)).*k;
TP_MultiRider_twoTags_1 = (3-sum(MultiRider_twoTags_1)).*k;
TP_MultiRider_twoTags_2 = (4-sum(MultiRider_twoTags_2)).*k;

TP_ParaFi_twoTags_0 = (2-sum(ParaFi_twoTags_0)).*k;
TP_ParaFi_twoTags_1 = (3-sum(ParaFi_twoTags_1)).*k;
TP_ParaFi_twoTags_2 = (4-sum(ParaFi_twoTags_2)).*k;

TP_FreeCollision_twoTags_0 = (2-sum(FreeCollision_twoTags_0)).*k;
TP_FreeCollision_twoTags_1 = (3-sum(FreeCollision_twoTags_1)).*k;
TP_FreeCollision_twoTags_2 = (4-sum(FreeCollision_twoTags_2)).*k;

TP_ConcurScatter_twoTags_0 = (2-sum(ConcurScatter_twoTags_0)).*k;
TP_ConcurScatter_twoTags_1 = (3-sum(ConcurScatter_twoTags_1)).*k;
TP_ConcurScatter_twoTags_2 = (4-sum(ConcurScatter_twoTags_2)).*k;

TP_Mecha = [TP_Mecha_twoTags_0;TP_Mecha_twoTags_1;TP_Mecha_twoTags_2];
TP_MultiRider = [TP_MultiRider_twoTags_0;TP_MultiRider_twoTags_1;TP_MultiRider_twoTags_2];
TP_ParaFi = [TP_ParaFi_twoTags_0;TP_ParaFi_twoTags_1;TP_ParaFi_twoTags_2];
TP_FreeCollision = [TP_FreeCollision_twoTags_0;TP_FreeCollision_twoTags_1;TP_FreeCollision_twoTags_2];
TP_ConcurScatter = [TP_ConcurScatter_twoTags_0;TP_ConcurScatter_twoTags_1;TP_ConcurScatter_twoTags_2];


%% plot BER
vals_BER = [BER_FreeCollision,BER_Mecha,BER_MultiRider,BER_ParaFi,BER_ConcurScatter];
vals_TP = [TP_FreeCollision,TP_Mecha,TP_MultiRider,TP_ParaFi,TP_ConcurScatter];


x = 0:1:2;
%%%** plot BER **%%%
figure(1)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b1 = bar(x,vals_BER);
b1(1).FaceColor = [69/255 123/255 157/255];
b1(2).FaceColor = [202/255 55/255 83/255];
b1(3).FaceColor =  [128/255 128/255 128/255];
b1(4).FaceColor = [153 153 204]./255; 
b1(5).FaceColor = [132/255 165/255 157/255]; 
ax1 = gca;
ax1.FontSize = 14;
ax1.FontName = 'Times New Roman';
ax1.YLim = [4e-8 1e0];
ax1.YTick = [1e-6,1e-4,1e-2,1e0];
ax1.YScale = 'log';
ax1.XTickLabel = {'0','1','2'};
grid on
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('The number of newly added tags','FontSize',14,'FontName','Times New Roman');

lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

%% Plot throughput
x = 0:1:2;
%%%** plot BER **%%%
figure(2)
set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
b2 = bar(x,vals_TP);
b2(1).FaceColor = [69/255 123/255 157/255];
b2(2).FaceColor = [202/255 55/255 83/255];
b2(3).FaceColor =  [128/255 128/255 128/255];
b2(4).FaceColor = [153 153 204]./255; 
b2(5).FaceColor = [132/255 165/255 157/255]; 
ax2 = gca;
ax2.FontSize = 14;
ax2.FontName = 'Times New Roman';
% ax2.YLim = [4e-6 1e0];
% ax2.YTick = [1e-6,1e-4,1e-2,1e0];
% ax2.YScale = 'log';
ax2.XTickLabel = {'0','1','2'};
grid on
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('The number of newly added tags','FontSize',14,'FontName','Times New Roman');

lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 


