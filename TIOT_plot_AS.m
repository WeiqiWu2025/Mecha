clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

% SNR = 30 dB; tag_count = 2; delay = 0:0.1:1;
%% Mecha
Mecha_twoTags_Tag1 = [1.60714285714286e-05;3.03571428571429e-05;5.80357142857143e-05;0.000372321428571429;0.000551785714285714;0.00106517857142857;0.00104910714285714;0.00136160714285714;0.00149464285714286;0.00185803571428571;0.00208392857142857];
Mecha_twoTags_Tag2 = [2.23214285714286e-05;0;0;1.78571428571429e-06;0;4.64285714285714e-05;3.39285714285714e-05;0;0;0;9.82142857142857e-06];

%% MultiRider
MultiRider_twoTags_Tag1 = [0;0;0;0;0;0;0;0;0;0;0];
MultiRider_twoTags_Tag2 = [8.92857142857143e-07;0;0;2.67857142857143e-06;0;3.57142857142857e-06;4.46428571428572e-06;3.57142857142857e-06;2.67857142857143e-06;8.92857142857143e-06;8.92857142857143e-06];

%% ParaFi
ParaFi_twoTags_Tag1 = [0.000150000000000000;0.000158035714285714;0.000433035714285714;0.0113928571428571;0.0158696428571429;0.0243973214285714;0.0325625000000000;0.0441660714285714;0.0489589285714286;0.0608928571428571;0.0633669642857143];
ParaFi_twoTags_Tag2 = [0.000183035714285714;0.000153571428571429;0.000125892857142857;0.000244642857142857;0.000316071428571429;0.000693750000000000;0.00128839285714286;0.00485714285714286;0.00500982142857143;0.0120500000000000;0.0137669642857143];

%% FreeCollision
FreeCollision_twoTags_Tag1 = [0.000204545454545455;0.00427727272727273;0.0115348484848485;0.0553234848484849;0.0571121212121212;0.0907136363636364;0.103487878787879;0.0859856060606061;0.0904113636363636;0.0963962121212121;0.127466666666667];
FreeCollision_twoTags_Tag2 = [0.000337121212121212;0.00546439393939394;0.0157681818181818;0.0701386363636364;0.0692939393939394;0.115155303030303;0.133810606060606;0.102770454545455;0.111821969696970;0.116209848484848;0.160575757575758];

%% ConcurScatter
ConcurScatter_twoTags_Tag1 = [0.499370394736842;0.499705263157895;0.501597368421053;0.499719736842105;0.498927631578947;0.500357236842105;0.501093421052632;0.499948684210526;0.501324342105263;0.500568421052632;0.500374342105263];
ConcurScatter_twoTags_Tag2 = [0.500127631578947;0.499607236842105;0.500210526315789;0.500200000000000;0.498441447368421;0.499601973684211;0.500615789473684;0.500462500000000;0.499649342105263;0.500084868421053;0.499417763157895];


%% Avg BER
BER_Mecha_twoTags = (Mecha_twoTags_Tag1 + Mecha_twoTags_Tag2)./2;
BER_MultiRider_twoTags = (MultiRider_twoTags_Tag1 + MultiRider_twoTags_Tag2)./2;
BER_ParaFi_twoTags = (ParaFi_twoTags_Tag1 + ParaFi_twoTags_Tag2)./2;
BER_FreeCollision_twoTags = (FreeCollision_twoTags_Tag1 + FreeCollision_twoTags_Tag2)./2;
BER_ConcurScatter_twoTags = (ConcurScatter_twoTags_Tag1 + ConcurScatter_twoTags_Tag2)./2;

%% Throughput
k = 250; % kbps
TP_Mecha_twoTags = (1-Mecha_twoTags_Tag1).*k +(1-Mecha_twoTags_Tag2).*k;
TP_MultiRider_twoTags = (1-MultiRider_twoTags_Tag1).*k +(1-MultiRider_twoTags_Tag2).*k;
TP_ParaFi_twoTags = (1-ParaFi_twoTags_Tag1).*k +(1-ParaFi_twoTags_Tag2).*k;
TP_FreeCollision_twoTags = (1-FreeCollision_twoTags_Tag1).*k +(1-FreeCollision_twoTags_Tag2).*k;
TP_ConcurScatter_twoTags = (1-ConcurScatter_twoTags_Tag1).*k +(1-ConcurScatter_twoTags_Tag2).*k;

%% plot BER
vals_BER = [BER_FreeCollision_twoTags,BER_Mecha_twoTags,BER_MultiRider_twoTags,BER_ParaFi_twoTags,BER_ConcurScatter_twoTags];
vals_TP = [TP_FreeCollision_twoTags,TP_Mecha_twoTags,TP_MultiRider_twoTags,TP_ParaFi_twoTags,TP_ConcurScatter_twoTags];

vals_BER(vals_BER == 0) = 1e-6;

colors = {[69/255 123/255 157/255],[202/255 55/255 83/255],[128/255 128/255 128/255],[153 153 204]./255,[132/255 165/255 157/255],[229,152,155]/255,[247,237,226]/255,[246,189,96]/255};
lineStyles = {'-','--','-.',':'};
lineMarker = {'square','o','^','pentagram','diamond','*','h','o','+','x'};
markersize = 8;
linewidth = 1.5;

x =0:0.1:1;

%%%** plot BER **%%%
figure(1)
plot(x,vals_BER(:,1),'Marker',lineMarker{1},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{1},'MarkerFaceColor',colors{1});
grid on
set(gcf,'unit','centimeters','position',[3 5 14.5 8]);
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.YScale = 'log';
% ax.XTickLabel = {'100','200','300','500','1000','1500','2000','3000'};
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('Time difference ($\mu$s)','FontSize',14,'FontName','Times New Roman','Interpreter','latex');
hold on
plot(x,vals_BER(:,2),'Marker',lineMarker{2},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{2},'MarkerFaceColor',colors{2});
hold on
plot(x,vals_BER(:,3),'Marker',lineMarker{3},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{3},'MarkerFaceColor',colors{3});
hold on
plot(x,vals_BER(:,4),'Marker',lineMarker{4},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{4},'MarkerFaceColor',colors{4});
hold on
plot(x,vals_BER(:,5),'Marker',lineMarker{5},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{5},'MarkerFaceColor',colors{5});
ax.XLim = [0 0.5];
ax.YLim = [1e-7 1e0];
ax.YTick = [1e-6,1e-4,1e-2,1e0];
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 


% figure(2)
% set(gcf,'unit','centimeters','position',[3 5 14.5 8]); 
% b2 = bar(x,vals_TP);
% b2(1).FaceColor = [69/255 123/255 157/255];
% b2(2).FaceColor = [202/255 55/255 83/255];
% b2(3).FaceColor =  [128/255 128/255 128/255];
% b2(4).FaceColor = [153 153 204]./255; 
% b2(5).FaceColor = [132/255 165/255 157/255]; 
% ax2 = gca;
% ax2.FontSize = 14;
% ax2.FontName = 'Times New Roman';
% ax2.XLim = [-0.05 0.55];
% ax2.YLim = [200 515];
% grid on
% ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
% xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');
% 
% lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
%     'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
% lgd.Location = 'northoutside'; 

figure(2)
plot(x,vals_TP(:,1),'Marker',lineMarker{1},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{1},'MarkerFaceColor',colors{1});
grid on
set(gcf,'unit','centimeters','position',[3 5 14.5 8]);
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
% ax.XTickLabel = {'100','200','300','500','1000','1500','2000','3000'};
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('Time difference ($\mu$s)','FontSize',14,'FontName','Times New Roman','Interpreter','latex');
hold on
plot(x,vals_TP(:,2),'Marker',lineMarker{2},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{2},'MarkerFaceColor',colors{2});
hold on
plot(x,vals_TP(:,3),'Marker',lineMarker{3},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{3},'MarkerFaceColor',colors{3});
hold on
plot(x,vals_TP(:,4),'Marker',lineMarker{4},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{4},'MarkerFaceColor',colors{4});
hold on
plot(x,vals_TP(:,5),'Marker',lineMarker{5},'MarkerSize',markersize,'LineWidth',linewidth,'color',colors{5},'MarkerFaceColor',colors{5});
ax.XLim = [0 0.5];
ax.YLim = [240 515];
% ax.YTick = [1e-6,1e-4,1e-2,1e0];
lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 





