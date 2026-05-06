clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

%% BER record
FreeCollision_threeTags_Tag1 = [0.0959303030303030;0.0386909090909091;0.00588939393939394;0.00109393939393939;0.000792424242424242];
FreeCollision_threeTags_Tag2 = [0.107446969696970;0.0428469696969697;0.00566666666666667;0.00174848484848485;0.00104696969696970];
FreeCollision_threeTags_Tag3 = [0.145737878787879;0.0590606060606061;0.00956060606060606;0.00257424242424242;0.00157575757575758];

MultiRider_threeTags_Tag1 = [0.00115217391304348;0.000100000000000000;0;0;0];
MultiRider_threeTags_Tag2 = [0.0184717391304348;0.00608695652173913;0.000643478260869565;8.04347826086957e-05;3.26086956521739e-05];
MultiRider_threeTags_Tag3 = [0.0361304347826087;0.0127347826086957;0.00133478260869565;0.000184782608695652;5.86956521739130e-05];

Mecha_threeTags_Tag1 = [0.0193717391304348;0.00554782608695652;0.000732608695652174;0.000173913043478261;6.52173913043478e-05];
Mecha_threeTags_Tag2 = [0.0198413043478261;0.00590434782608696;0.000693478260869565;0.000141304347826087;6.52173913043478e-05];
Mecha_threeTags_Tag3 = [0.0189152173913043;0.00425869565217391;0.000284782608695652;0.000100000000000000;8.26086956521739e-05];

ParaFi_threeTags_Tag1 = [0.215773913043478;0.116173913043478;0.0114260869565217;0.00116304347826087;0.000241304347826087];
ParaFi_threeTags_Tag2 = [0.214241304347826;0.115023913043478;0.0121565217391304;0.00130000000000000;0.000250000000000000];
ParaFi_threeTags_Tag3 = [0.215723913043478;0.114828260869565;0.0107760869565217;0.000776086956521739;0.000150000000000000];

ConcurScatter_threeTags_Tag1 = [0.498973684210526;0.499111842105263;0.499539473684211;0.499876315789474;0.499828947368421];
ConcurScatter_threeTags_Tag2 = [0.499486842105263;0.498951315789474;0.498798684210526;0.499102631578947;0.499086842105263];
ConcurScatter_threeTags_Tag3 = [0.499747368421053;0.499136842105263;0.499497368421053;0.499125000000000;0.498959210526316];


FreeCollision_threeTags = (FreeCollision_threeTags_Tag1 + FreeCollision_threeTags_Tag2 + FreeCollision_threeTags_Tag3)./3;
MultiRider_threeTags = (MultiRider_threeTags_Tag1 + MultiRider_threeTags_Tag2 + MultiRider_threeTags_Tag3)./3;
Mecha_threeTags = (Mecha_threeTags_Tag1 + Mecha_threeTags_Tag2 + Mecha_threeTags_Tag3)./3;
ParaFi_threeTags = (ParaFi_threeTags_Tag1 + ParaFi_threeTags_Tag2 + ParaFi_threeTags_Tag3)./3;
ConcurScatter_threeTags = (ConcurScatter_threeTags_Tag1 + ConcurScatter_threeTags_Tag2 + ConcurScatter_threeTags_Tag3)./3;

%% Throughput
k = 250; % kbps
TP_FreeCollision_threeTags = (1-FreeCollision_threeTags_Tag1).*k + (1-FreeCollision_threeTags_Tag2).*k + (1-FreeCollision_threeTags_Tag3).*k;
TP_MultiRider_threeTags = (1-MultiRider_threeTags_Tag1).*k + (1-MultiRider_threeTags_Tag2).*k + (1-MultiRider_threeTags_Tag3).*k;
TP_ParaFi_threeTags = (1-ParaFi_threeTags_Tag1).*k + (1-ParaFi_threeTags_Tag2).*k + (1-ParaFi_threeTags_Tag3).*k;
TP_ConcurScatter_threeTags = (1-ConcurScatter_threeTags_Tag1).*k + (1-ConcurScatter_threeTags_Tag2).*k + (1-ConcurScatter_threeTags_Tag3).*k;
TP_Mecha_threeTags = (1-Mecha_threeTags_Tag1).*k + (1-Mecha_threeTags_Tag2).*k + (1-Mecha_threeTags_Tag3).*k;


%% BER plot
vals_BER = [FreeCollision_threeTags,Mecha_threeTags,MultiRider_threeTags,ParaFi_threeTags,ConcurScatter_threeTags];
vals_TP = [TP_FreeCollision_threeTags,TP_Mecha_threeTags,TP_MultiRider_threeTags,TP_ParaFi_threeTags,TP_ConcurScatter_threeTags];

x = 1:1:5;
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
ax1.YLim = [4e-6 1e0];
ax1.YTick = [1e-6,1e-4,1e-2,1e0];
ax1.YScale = 'log';
ax1.XTickLabel = {'0','5','15','25','35'};
grid on
ylabel('BER','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');

lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 

%% Plot throughput
x = 1:1:5;
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
ax2.XTickLabel = {'0','5','15','25','35'};
grid on
ylabel('Throughput (Kbps)','FontSize',14,'FontName','Times New Roman');
xlabel('Backscattered signal SNR (dB)','FontSize',14,'FontName','Times New Roman');

lgd = legend('FreeCollision','Mecha','MultiRider','ParaFi','ConcurScatter',...
    'FontSize',12,'FontName','Times New Roman','NumColumns',3,'Interpreter','latex');
lgd.Location = 'northoutside'; 


