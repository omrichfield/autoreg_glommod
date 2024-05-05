clear all
close all

load("train_data");
load("autoreg_params");
load("Myogenic_tone_fit_data");
load("TGF_tone_fit_data");
load("Myo_TGF_model_curves_20230911.mat");

P_aa = 100:25:150;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
plot(P_aa,train_data.D_cont,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(P_aa,train_data.D_furo,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k');
plot(P_aa,train_data.D_dilt,'kv','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(Try_Pf,D_pred_dilt,'k-','LineWidth',3);
plot(Try_Pf,D_pred_furo,'b-','LineWidth',3);
plot(Try_Pf,D_pred_cont,'m--','LineWidth',3);
plot(Try_Pf,D_pred_cont_int,'r-','LineWidth',3);
hold off;
set(gca, 'box','off','FontSize',15,'YTick',5:1:11);
xlim([90 160]);
ylim([5 11]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("D_{AA} (\mu{m})",'FontSize',20);
axis square

   % saveas(gcf,'figures/P_D_takenaka.fig')
   % saveas(gcf,'figures/P_D_takenaka','epsc')
   % saveas(gcf,'figures/P_D_takenaka','jpeg')
   

subplot(2,3,2)
plot(P_aa,train_data.Q_cont,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(P_aa,train_data.Q_furo,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k');
plot(P_aa,train_data.Q_dilt,'kv','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');

plot(Try_Pf,Q_pred_dilt/0.6,'k-','LineWidth',3);
plot(Try_Pf,Q_pred_furo/0.6,'b-','LineWidth',3);
plot(Try_Pf,Q_pred_cont/0.6,'m--','LineWidth',3);
plot(Try_Pf,Q_pred_cont_int/0.6,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
ylim([100 650]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Q_{AA} (nl/min)",'FontSize',20);
axis square

subplot(2,3,3)
plot(P_aa,train_data.Q_cont,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(P_aa,train_data.Q_furo,'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k');
plot(P_aa,train_data.Q_dilt,'kv','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');

plot(Try_Pf,Q_pred_dilt/0.6,'k-','LineWidth',3);
plot(Try_Pf,Q_pred_furo/0.6,'b-','LineWidth',3);
plot(Try_Pf,Q_pred_cont/0.6,'m--','LineWidth',3);
plot(Try_Pf,Q_pred_cont_int/0.6,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
ylim([1000 100000]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Afferent Arteriole Blood Flow (nl/min)",'FontSize',20);
legend('Control, Takenaka 1994', ...
        'Furosemide, Takenaka 1994', ...
        'Diltiazem, Takenaka 1994', ...
        'Model Passive',...
        'Model Myo', ...
        'Model Myo + TGF', ...
        'Model Myo* + TGF', ...
        'location','west')
axis off

subplot(2,3,4)
plot(Try_Pf,PG_dilt,'k-','LineWidth',3); hold on;
plot(Try_Pf,PG_furo,'b-','LineWidth',3);
plot(Try_Pf,PG_cont,'m--','LineWidth',3);
plot(Try_Pf,PG_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
ylim([40 110]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("P_G (mmHg)",'FontSize',20);
axis square

subplot(2,3,5)
plot(Try_Pf,SNGFR_dilt,'k-','LineWidth',3); hold on;
plot(Try_Pf,SNGFR_furo,'b-','LineWidth',3);
plot(Try_Pf,SNGFR_cont,'m--','LineWidth',3);
plot(Try_Pf,SNGFR_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',0:30:180);
xlim([90 160]);
ylim([10 180]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("SNGFR (nl/min)",'FontSize',20);
axis square

% subplot(2,3,6)
% plot(Try_Pf,C_MD_dilt,'k-','LineWidth',3); hold on;
% plot(Try_Pf,C_MD_furo,'b-','LineWidth',3);
% plot(Try_Pf,C_MD_cont,'m-','LineWidth',3);
% plot(Try_Pf,C_MD_int,'r-','LineWidth',3);hold off;
% set(gca, 'box','off','FontSize',15);
% xlim([90 160]);
% ylim([0 650]);
% xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
% ylabel("C_{MD} (mM)",'FontSize',20);
% axis square

    saveas(gcf,'figures/P_QD_takenaka.fig')
    saveas(gcf,'figures/P_QD_takenaka','epsc')
    saveas(gcf,'figures/P_QD_takenaka','jpeg')


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1);  
patch([Try_Pf' fliplr(Try_Pf')], [shear_mn_furo' fliplr(shear_mn_dilt')],...
    [169/255 255/255 142/255], 'EdgeColor',[169/255 255/255 142/255]); hold on;
patch([Try_Pf' fliplr(Try_Pf')], [shear_mn_cont' fliplr(shear_mn_furo')],...
    [165/255 199/255 255/255], 'EdgeColor',[169/255 255/255 142/255]); 
patch([Try_Pf' fliplr(Try_Pf')], [shear_mn_int' fliplr(shear_mn_cont')],...
    [238/255 180/255 180/255], 'EdgeColor',[238/255 180/255 180/255]); 

plot(Try_Pf,shear_mn_dilt,'k-','LineWidth',3); 
plot(Try_Pf,shear_mn_furo,'b-','LineWidth',3); 
plot(Try_Pf,shear_mn_cont,'m--','LineWidth',3);
plot(Try_Pf,shear_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',10:20:90);
xlim([90 160]);
ylim([10 90]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Shear Stress (dynes/cm^2)",'FontSize',20);
axis square

subplot(2,4,2)    
patch([Try_Pf' fliplr(Try_Pf')], [hoop_mn_furo' fliplr(hoop_mn_dilt')],...
    [169/255 255/255 142/255], 'EdgeColor',[169/255 255/255 142/255]); hold on;
patch([Try_Pf' fliplr(Try_Pf')], [hoop_mn_cont' fliplr(hoop_mn_furo')],...
    [165/255 199/255 255/255], 'EdgeColor',[169/255 255/255 142/255]); 
patch([Try_Pf' fliplr(Try_Pf')], [hoop_mn_int' fliplr(hoop_mn_cont')],...
    [238/255 180/255 180/255], 'EdgeColor',[238/255 180/255 180/255]); 

plot(Try_Pf,hoop_mn_dilt,'k-','LineWidth',3); hold on;
plot(Try_Pf,hoop_mn_furo,'b-','LineWidth',3); 
plot(Try_Pf,hoop_mn_cont,'m--','LineWidth',3);
plot(Try_Pf,hoop_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',50:20:170);
xlim([90 160]);
ylim([50 170]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Hoop Stress (kPa)",'FontSize',20);
axis square

subplot(2,4,3)   
patch([Try_Pf' fliplr(Try_Pf')], [CSGFR_mn_furo' fliplr(CSGFR_mn_dilt')],...
    [169/255 255/255 142/255], 'EdgeColor',[169/255 255/255 142/255]); hold on;
patch([Try_Pf' fliplr(Try_Pf')], [CSGFR_mn_cont' fliplr(CSGFR_mn_furo')],...
    [165/255 199/255 255/255], 'EdgeColor',[169/255 255/255 142/255]); 
patch([Try_Pf' fliplr(Try_Pf')], [CSGFR_mn_int' fliplr(CSGFR_mn_cont')],...
    [238/255 180/255 180/255], 'EdgeColor',[238/255 180/255 180/255]); 

plot(Try_Pf,CSGFR_mn_dilt,'k-','LineWidth',3); hold on;
plot(Try_Pf,CSGFR_mn_furo,'b-','LineWidth',3); 
plot(Try_Pf,CSGFR_mn_cont,'m--','LineWidth',3);
plot(Try_Pf,CSGFR_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);%,'YTick',-1:1:4);
xlim([90 160]);
%ylim([0 0.2]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("CSGFR (nl/min)",'FontSize',20);
axis square

subplot(2,4,5)    
plot(Try_Pf,(shear_mn_dilt-shear_mn_furo),'Color',[169/255 255/255 142/255],'LineWidth',5); hold on;
plot(Try_Pf,(shear_mn_furo-shear_mn_cont),'Color',[165/255 199/255 255/255],'LineWidth',5); 
plot(Try_Pf,(shear_mn_cont-shear_mn_int),'Color',[238/255 180/255 180/255],'LineWidth',5); hold off;
%plot(Try_Pf,shear_mn_cont,'m-','LineWidth',3);
%plot(Try_Pf,shear_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',0:10:100);
xlim([90 160]);
%ylim([0 100]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Shear Stress Contribution (dynes/cm^2)",'FontSize',17);
axis square

subplot(2,4,6)    
plot(Try_Pf,(hoop_mn_dilt-hoop_mn_furo),'Color',[169/255 255/255 142/255],'LineWidth',5); hold on;
plot(Try_Pf,(hoop_mn_furo-hoop_mn_cont),'Color',[165/255 199/255 255/255],'LineWidth',5); 
plot(Try_Pf,(hoop_mn_cont-hoop_mn_int),'Color',[238/255 180/255 180/255],'LineWidth',5); hold off;
%plot(Try_Pf,shear_mn_cont,'m-','LineWidth',3);
%plot(Try_Pf,shear_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',0:20:100);
xlim([90 160]);
%ylim([0 100]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Hoop Stress Contribution (kPa)",'FontSize',17);
axis square

subplot(2,4,7)    
plot(Try_Pf,(CSGFR_mn_dilt-CSGFR_mn_furo),'Color',[169/255 255/255 142/255],'LineWidth',5); hold on;
plot(Try_Pf,(CSGFR_mn_furo-CSGFR_mn_cont),'Color',[165/255 199/255 255/255],'LineWidth',5); 
plot(Try_Pf,(CSGFR_mn_cont-CSGFR_mn_int),'Color',[238/255 180/255 180/255],'LineWidth',5); hold off;
%plot(Try_Pf,shear_mn_cont,'m-','LineWidth',3);
%plot(Try_Pf,shear_mn_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15,'YTick',0:0.1:0.5);
xlim([90 160]);
%ylim([0 100]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("CSGFR Contribution (nl/min)",'FontSize',17);
axis square

subplot(2,4,4)
plot(Try_Pf,Q_pred_dilt/0.6,'k-','LineWidth',3);hold on;
plot(Try_Pf,Q_pred_furo/0.6,'b-','LineWidth',3); 
plot(Try_Pf,Q_pred_cont/0.6,'m--','LineWidth',3);
plot(Try_Pf,Q_pred_cont_int/0.6,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
ylim([1000 100000]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Afferent Arteriole Blood Flow (nl/min)",'FontSize',20);
legend('Model Passive',...
        'Model Myo', ...
        'Model Myo + TGF', ...
        'Model Myo* + TGF', ...
        'location','west')
axis off

subplot(2,4,8)
plot(Try_Pf,Q_pred_dilt/0.6,'Color',[169/255 255/255 142/255],'LineWidth',3);hold on;
plot(Try_Pf,Q_pred_furo/0.6,'Color',[165/255 199/255 255/255],'LineWidth',3); 
plot(Try_Pf,Q_pred_cont/0.6,'Color',[238/255 180/255 180/255],'LineWidth',3);
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
ylim([1000 100000]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Afferent Arteriole Blood Flow (nl/min)",'FontSize',20);
legend('Myo Contribution', ...
        'TGF Contribution', ...
        'Myo* (Intxn) Contribution', ...
        'location','west')
axis off

    saveas(gcf,'figures/P_mech_takenaka.fig')
    saveas(gcf,'figures/P_mech_takenaka','epsc')
    saveas(gcf,'figures/P_mech_takenaka','jpeg')

P_aa = 100:25:150;
figure('units','normalized','outerposition',[0 0 0.8 0.55])
subplot(1,3,2)
plot(Try_Pf,S_Myo_furo,'b-','LineWidth',3); hold on;
plot(Try_Pf,S_Myo_cont,'m-','LineWidth',3);
plot(Try_Pf,S_Myo_cont_int,'r-','LineWidth',3);
plot(Try_Pf,S_TGF_cont,'m--','LineWidth',3);
plot(Try_Pf,S_TGF_cont_int,'r--','LineWidth',3); hold off;
hold off;
set(gca, 'box','off','FontSize',15,'YTick',-1:1:4);
xlim([90 160]);
ylim([-1 4]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Tone",'FontSize',20);
axis square

subplot(1,3,1)
plot(Try_Pf,SNGFR_dilt,'k-','LineWidth',3); hold on;
plot(Try_Pf,SNGFR_furo,'b-','LineWidth',3);
plot(Try_Pf,SNGFR_cont,'m-','LineWidth',3);
plot(Try_Pf,SNGFR_int,'r-','LineWidth',3);hold off;
set(gca, 'box','off','FontSize',15);
xlim([90 160]);
%ylim([100 650]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("SNGFR (nl/min)",'FontSize',20);
axis square

   % saveas(gcf,'figures/P_D_takenaka.fig')
   % saveas(gcf,'figures/P_D_takenaka','epsc')
   % saveas(gcf,'figures/P_D_takenaka','jpeg')
    


   load("Tubule_schem_L_mat.mat");
   load("Tubule_schem_P_mat.mat");

figure()

tt = tiledlayout(1,1);
ax1=axes(tt);
ax2=axes(tt);


dummy_vec = Tubule_schem_L(:,3)*5000;
plot(ax2,Tubule_schem_L(:,1), dummy_vec,'b','LineWidth',3);hold on;
plot(ax1,Tubule_schem_L(:,3),Tubule_schem_L(:,3));
plot(ax2,Tubule_schem_L(:,1), Tubule_schem_L(:,2),'k','LineWidth',3);
xline(0.5,'k--','LineWidth',2);
xline(0.9,'k--','LineWidth',2);
hold off;
legend('Velocity','Osmolality','location','SouthWest')

%scatter(mtfd.Tp(2:4), mtfd.S_tone(2:4),1000,'k.');
ax1.YAxisLocation='right';
set(ax2, 'FontSize',15,'XTick',0:0.2:1.4);
set(ax1, 'FontSize',15,'YTick',0:0.02:0.12,'box','on','XTick',[]);
%xticks(ax1,0:50:200);
text(0.15,550,'PCT','fontsize',15);
text(0.15+0.5,550,'DL','fontsize',15);
text(0.15+1,550,'AL','fontsize',15);
xlim([0 1.4]);
%ylim([-1 4]);
xlabel("Length along Tubule (cm)",'FontSize',20);
ylabel(ax2,"Osmolality (mosmol/kg H_2O)",'FontSize',20);
ylabel(ax1,"Fluid Velocity (cm/s)",'FontSize',20);

    saveas(gcf,'figures/tubule_schem_L.fig')
    saveas(gcf,'figures/tubule_schem_L','epsc')
    saveas(gcf,'figures/tubule_schem_L','jpeg')

figure()

tt = tiledlayout(1,1);
ax1=axes(tt);
ax2=axes(tt);


dummy_vec = Tubule_schem(2,:)*5000;
plot(ax2,Tubule_schem(1,:), dummy_vec,'b','LineWidth',3);hold on;
plot(ax1,Tubule_schem(2,:),Tubule_schem(2,:));
plot(ax2,Tubule_schem(1,:), Tubule_schem(3,:),'k','LineWidth',3);
hold off;
legend('Velocity','Osmolality','location','SouthEast')

%scatter(mtfd.Tp(2:4), mtfd.S_tone(2:4),1000,'k.');
ax1.YAxisLocation='right';
set(ax2, 'FontSize',15,'XTick',90:10:160);
set(ax1, 'FontSize',15,'YTick',0:0.02:0.12,'box','on','XTick',[]);
%xticks(ax1,0:50:200);
%text(0.15,550,'PCT','fontsize',15);
%text(0.15+0.5,550,'DL','fontsize',15);
%ext(0.15+1,550,'AL','fontsize',15);
xlim([100 150]);
ylim(ax1,[0 0.12]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel(ax2,"Osmolality (mosmol/kg H_2O)",'FontSize',20);
ylabel(ax1,"Fluid Velocity (cm/s)",'FontSize',20);

    saveas(gcf,'figures/tubule_schem_P.fig')
    saveas(gcf,'figures/tubule_schem_P','epsc')
    saveas(gcf,'figures/tubule_schem_P','jpeg')

   
figure()
plot(Tubule_schem(1,:), Tubule_schem(4,:),'k','LineWidth',3);hold on;
plot(Tubule_schem(1,:),Tubule_schem(4,:)-Tubule_schem(5,:), 'b','LineWidth',3);
hold off;
legend('SNGFR','Reabsorbed Fluid','location','SouthEast')

%scatter(mtfd.Tp(2:4), mtfd.S_tone(2:4),1000,'k.');
set(gca,'FontSize',15,'XTick',90:10:160,'box','off');
xlim([100 150]);
ylim([0 70]);
xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
ylabel("Fluid Flow (nl/min)",'FontSize',20);
axis square

    saveas(gcf,'figures/tubule_schem_P_2.fig')
    saveas(gcf,'figures/tubule_schem_P_2','epsc')
    saveas(gcf,'figures/tubule_schem_P_2','jpeg')

% subplot(1,3,3)
% plot(Try_Pf,SNGFR_dilt,'k-','LineWidth',3); hold on;
% plot(Try_Pf,SNGFR_furo,'b-','LineWidth',3);
% plot(Try_Pf,SNGFR_cont,'m-','LineWidth',3);
% plot(Try_Pf,SNGFR_int,'r-','LineWidth',3);hold off;
% set(gca, 'box','off','FontSize',15);
% xlim([90 160]);
% ylim([20 70]);
% xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
% ylabel("Afferent Arteriole Blood Flow (nl/min)",'FontSize',20);
% axis square
%     saveas(gcf,'figures/P_Q_takenaka.fig')
%     saveas(gcf,'figures/P_Q_takenaka','epsc')
%     saveas(gcf,'figures/P_Q_takenaka','jpeg')



%     figure(5)
% plot(Try_Pf,S_Myo_furo,'b-','LineWidth',3);hold on;
% plot(Try_Pf,S_Myo_cont,'m-','LineWidth',3);
% plot(Try_Pf,S_Myo_cont_int,'r-','LineWidth',3);hold off;
% legend('Myo', ...
%         'Myo with TGF', ...
%         'Myo* with TGF', ...
%         'location','northwest')
% set(gca, 'box','off','FontSize',15);
% xlim([90 160]);
% ylim([-1 4]);
% xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
% ylabel("Myogenic Tone",'FontSize',20);
%     saveas(gcf,'figures/P_S_Myo.fig')
%     saveas(gcf,'figures/P_S_Myo','epsc')
%     saveas(gcf,'figures/P_S_Myo','jpeg')
% 
%     figure(6)
% plot(Try_Pf,S_TGF_furo*0,'b-','LineWidth',3);hold on;
% plot(Try_Pf,S_TGF_cont,'m-','LineWidth',3);
% plot(Try_Pf,S_TGF_cont_int,'r-','LineWidth',3);hold off;
% legend('No TGF', ...
%         'TGF with Myo', ...
%         'TGF with Myo*', ...
%         'location','northwest')
% set(gca, 'box','off','FontSize',15);
% xlim([90 160]);
% ylim([-1 4]);
% xlabel("Perfusion Pressure (mmHg)",'FontSize',20);
% ylabel("TGF Tone",'FontSize',20);
%     saveas(gcf,'figures/P_S_TGF.fig')
%     saveas(gcf,'figures/P_S_TGF','epsc')
%     saveas(gcf,'figures/P_S_TGF','jpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
C_try = 0:180;
Tp_try = 200:450;
ap = autoreg_params;

S_tone_myo = ap.Splus_myo./(1+exp(-ap.C_myo.*(Tp_try - ap.Tp_myo))) + ap.Smin_myo;
S_tone_TGF = ap.Splus_TGF./(1+exp(-ap.C_TGF.*(C_try - ap.C_MD_TGF))) + ap.Smin_TGF;


figure('units','normalized','outerposition',[0 0 0.6 0.5])
%tt = tiledlayout(1,1);
%ax1=axes(tt);
%ax2=axes(tt);

%dummy_vec = C_MD_show*2.5+100;
subplot(1,2,1)
%plot(ax2,dummy_vec, S_TGF_show,'k','LineWidth',3);hold on;
%plot(ax1,C_MD_show,C_MD_show);
plot(Tp_show, S_Myo_show_furo,'b','LineWidth',3); hold on;
plot(Tp_show, S_Myo_show_int,'r','LineWidth',3); hold off;
legend('Myo (No TGF)','Myo* (TGF)','location','SouthEast')

%scatter(mtfd.Tp(2:4), mtfd.S_tone(2:4),1000,'k.');
%ax1.XAxisLocation='top';
set(gca, 'FontSize',15);
set(gca, 'FontSize',15,'box','off');
%xticks(ax1,0:50:200);

xlim([100 600]);
ylim([-1 4]);
xlabel("Afferent Arteriole Wall Tension (mmHg \mu m)",'FontSize',20);
ylabel("Tone",'FontSize',20);

subplot(1,2,2)
plot(C_MD_show, S_TGF_show,'k','LineWidth',3);hold on;
%plot(ax1,C_MD_show,C_MD_show);
%plot(Tp_show, S_Myo_show_furo,'b','LineWidth',3); hold on;
%plot(Tp_show, S_Myo_show_int,'r','LineWidth',3); hold off;
legend('TGF','location','SouthEast')

%scatter(mtfd.Tp(2:4), mtfd.S_tone(2:4),1000,'k.');
%ax1.XAxisLocation='top';
set(gca, 'FontSize',15);
set(gca, 'FontSize',15,'box','off');
ylim([-1 4]);
xlabel("Macula Densa Osmolality (mOsmol/kg H_2O)",'FontSize',20);

ylabel("Tone",'FontSize',20);

    saveas(gcf,'figures/autoreg_curves.fig')
    saveas(gcf,'figures/autoreg_curves','epsc')
    saveas(gcf,'figures/autoreg_curves','jpeg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    TRANSIENT   MECHANICS   %%%%%%%%%%%%%%%%%%%%
load("transient_mech_20230914.mat")
load("transient_mech_130_20230915.mat")

load("transient_mech_100_20231012.mat")
load("transient_mech_130_20231012.mat")

figure('units','normalized','outerposition',[0 0 1 0.8]);
% 
 set(gca,'box','off','FontSize',20)
% 
label_1 = { 'Passive, 100','Passive, 130',...
                    'Myo, 100', 'Myo, 130',...
                    'TGF, 100', 'TGF, 130',...
                    'Myo + TGF, 100','Myo + TGF, 130',...
                    'Myo* + TGF, 100','Myo* + TGF, 130'};
label = categorical(label_1);
label = reordercats(label,label_1);

strain_100 = strain_100([3 2 5 1 4],:)';
delt_shear_100 = delt_shear_100([3 2 5 1 4],:)';
delt_CSGFR_100 = delt_CSGFR_100([3 2 5 1 4],:)';

strain_100_myo_contrib = (strain_100(:,1)-strain_100(:,2));%./(strain_100(:,1)-strain_100(:,4))*100;
strain_100_TGF_contrib = (strain_100(:,2)-strain_100(:,3));%./(strain_100(:,1)-strain_100(:,4))*100;
strain_100_int_contrib = (strain_100(:,3)-strain_100(:,4));%./(strain_100(:,1)-strain_100(:,4))*100;

strain_100_contrib = [strain_100_myo_contrib strain_100_TGF_contrib strain_100_int_contrib];

delt_shear_100_myo_contrib = (delt_shear_100(:,1)-delt_shear_100(:,2));%./(delt_shear_100(:,1)-delt_shear_100(:,4))*100;
delt_shear_100_TGF_contrib = (delt_shear_100(:,2)-delt_shear_100(:,3));%./(delt_shear_100(:,1)-delt_shear_100(:,4))*100;
delt_shear_100_int_contrib = (delt_shear_100(:,3)-delt_shear_100(:,4));%./(delt_shear_100(:,1)-delt_shear_100(:,4))*100;

delt_shear_100_contrib = [delt_shear_100_myo_contrib delt_shear_100_TGF_contrib delt_shear_100_int_contrib];

delt_CSGFR_100_myo_contrib = (delt_CSGFR_100(:,1)-delt_CSGFR_100(:,2));%./(delt_CSGFR_100(:,1)-delt_CSGFR_100(:,4))*100;
delt_CSGFR_100_TGF_contrib = (delt_CSGFR_100(:,2)-delt_CSGFR_100(:,3));%./(delt_CSGFR_100(:,1)-delt_CSGFR_100(:,4))*100;
delt_CSGFR_100_int_contrib = (delt_CSGFR_100(:,3)-delt_CSGFR_100(:,4));%./(delt_CSGFR_100(:,1)-delt_CSGFR_100(:,4))*100;

delt_CSGFR_100_contrib = [delt_CSGFR_100_myo_contrib delt_CSGFR_100_TGF_contrib delt_CSGFR_100_int_contrib];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strain_130 = strain_130([3 2 5 1 4],:)';
delt_shear_130 = delt_shear_130([3 2 5 1 4],:)';
delt_CSGFR_130 = delt_CSGFR_130([3 2 5 1 4],:)';

strain_130_myo_contrib = (strain_130(:,1)-strain_130(:,2));%./(strain_130(:,1)-strain_130(:,4))*100;
strain_130_TGF_contrib = (strain_130(:,2)-strain_130(:,3));%./(strain_130(:,1)-strain_130(:,4))*100;
strain_130_int_contrib = (strain_130(:,3)-strain_130(:,4));%./(strain_130(:,1)-strain_130(:,4))*100;

strain_130_contrib = [strain_130_myo_contrib strain_130_TGF_contrib strain_130_int_contrib];

delt_shear_130_myo_contrib = (delt_shear_130(:,1)-delt_shear_130(:,2));%./(delt_shear_130(:,1)-delt_shear_130(:,4))*100;
delt_shear_130_TGF_contrib = (delt_shear_130(:,2)-delt_shear_130(:,3));%./(delt_shear_130(:,1)-delt_shear_130(:,4))*100;
delt_shear_130_int_contrib = (delt_shear_130(:,3)-delt_shear_130(:,4));%./(delt_shear_130(:,1)-delt_shear_130(:,4))*100;

delt_shear_130_contrib = [delt_shear_130_myo_contrib delt_shear_130_TGF_contrib delt_shear_130_int_contrib];

delt_CSGFR_130_myo_contrib = (delt_CSGFR_130(:,1)-delt_CSGFR_130(:,2));%./(delt_CSGFR_130(:,1)-delt_CSGFR_130(:,4))*100;
delt_CSGFR_130_TGF_contrib = (delt_CSGFR_130(:,2)-delt_CSGFR_130(:,3));%./(delt_CSGFR_130(:,1)-delt_CSGFR_130(:,4))*100;
delt_CSGFR_130_int_contrib = (delt_CSGFR_130(:,3)-delt_CSGFR_130(:,4));%./(delt_CSGFR_130(:,1)-delt_CSGFR_130(:,4))*100;

delt_CSGFR_130_contrib = [delt_CSGFR_130_myo_contrib delt_CSGFR_130_TGF_contrib delt_CSGFR_130_int_contrib];

strain = 0*[strain_100 strain_130];
strain(:,1:2:end)=strain_100;
strain(:,2:2:end)=strain_130;

delt_shear = 0*[delt_shear_100 delt_shear_130];
delt_shear(:,1:2:end)=delt_shear_100;
delt_shear(:,2:2:end)=delt_shear_130;

delt_CSGFR = 0*[delt_CSGFR_100 delt_CSGFR_130];
delt_CSGFR(:,1:2:end)=delt_CSGFR_100;
delt_CSGFR(:,2:2:end)=delt_CSGFR_130;

delt_shear(delt_shear < 1) = NaN;

color_mat=repmat([0 0 0],length(label_1),1);
color_mat(2:2:end,1)=1;
axis square

subplot(1,3,1)
boxplot(strain, 'Labels',label_1,'Symbol','o','color',color_mat);%,'BoxStyle','filled'); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 3]);
%set(a, 'Color', color_mat);   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("Circumferential Strain, \epsilon_\theta (%)");
axis square

subplot(1,3,2)
boxplot(delt_shear, 'Labels',label_1,'Symbol','o','color',color_mat); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 5]);
set(gca, 'YScale', 'log')
%set(a, 'Color', 'k');   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("\Delta{\tau} (dynes/cm^2)");
axis square


subplot(1,3,3)
boxplot(delt_CSGFR, 'Labels',label_1,'Symbol','o','color',color_mat); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 5]);
%set(a, 'Color', 'k');   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("\Delta{CSGFR} (nl/min)");
axis square

    saveas(gcf,'figures/mech_box.fig')
    saveas(gcf,'figures/mech_box','epsc')
    saveas(gcf,'figures/mech_box','jpeg')

strain_contrib = 0*[strain_100_contrib strain_130_contrib];
strain_contrib(:,1:2:end)=strain_100_contrib;
strain_contrib(:,2:2:end)=strain_130_contrib;

delt_shear_contrib = 0*[delt_shear_100_contrib delt_shear_130_contrib];
delt_shear_contrib(:,1:2:end)=delt_shear_100_contrib;
delt_shear_contrib(:,2:2:end)=delt_shear_130_contrib;

delt_CSGFR_contrib = 0*[delt_CSGFR_100_contrib delt_CSGFR_130_contrib];
delt_CSGFR_contrib(:,1:2:end)=delt_CSGFR_100_contrib;
delt_CSGFR_contrib(:,2:2:end)=delt_CSGFR_130_contrib;

delt_CSGFR_contrib(abs(delt_CSGFR_contrib)>100) = NaN;


figure('units','normalized','outerposition',[0 0 1 0.8]);

label_2 = {'Myogenic, 100','Myogenic, 130', 'TGF, 100', 'TGF,130', 'Myo* (Intxn), 100', 'Myo* (Intxn), 130'};

subplot(1,3,1)
boxplot(abs(strain_contrib), 'Labels',label_2,'Symbol','o','color',color_mat); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 100]);
%set(a, 'Color', 'k');   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("Strain Contribution (%)");
axis square

subplot(1,3,2)
boxplot(abs(delt_shear_contrib), 'Labels',label_2,'Symbol','o','color',color_mat); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 100]);
%set(a, 'Color', 'k');   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("\Delta{\tau} Contribution (dynes/cm^2)");
axis square

%indx = findRows(delt_CSGFR_100_contrib>0 & delt_CSGFR_100_contrib<100) ;

subplot(1,3,3)
boxplot(abs(delt_CSGFR_contrib), 'Labels',label_2,'Symbol','o','color',color_mat); 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
%ylim([0 110]);
%set(a, 'Color', 'k');   % Set the color of the first box to green
set(gca,'box','off','FontSize',20)
ylabel("\Delta{CSGFR} Contribution (nl/min)");
axis square

    saveas(gcf,'figures/mech_contrib.fig')
    saveas(gcf,'figures/mech_contrib','epsc')
    saveas(gcf,'figures/mech_contrib','jpeg')



load("G_base_20230925.mat")

min_strain_130 = min(min(strain_130(strain_130>0)));
strain_130norm = strain_130/min_strain_130;
strain_130norm = [strain_130norm(2,:); strain_130norm(2:end-1,:); strain_130norm(end-1,:)];

sn=G_base.src;
tn=G_base.trg;
[sn,I] = sort(sn);
tn = tn(I);
Ggraph = digraph(sn, tn);

figure('units','normalized','outerposition',[0 0 0.9 1])

for i=1:4
subplot(2,2,i)

if i == 4
    i=5;
end

color_mat = [0*strain_100(:,3) 0*strain_100(:,3) 0*strain_100(:,3)];
color_mat(:,1) = (strain_130norm(:,i))/(max(max(strain_130norm)));
color_mat(:,2)= 0;
color_mat(:,3) = 0.5;%1-(strain_130norm(:,i))/(max(max(strain_130norm)));

size_vec = strain_130norm(:,i)/max(max(strain_130norm))*23;
size_vec(size_vec==0)=min(size_vec(size_vec>0));


p1=plot(Ggraph,'k', 'Layout', 'layered', ...
    'EdgeColor', color_mat, ...
    'LineWidth',size_vec,...
    'ArrowSize', 0); hold on;
p2=plot(Ggraph,'k', 'Layout', 'layered', ...
    'ArrowSize', 12); hold off;
text(p2.XData(195),p2.YData(195)-1.5,'EA','fontsize',30);
text(p2.XData(1),p2.YData(1)+1.5,'AA','fontsize',30);
set(gca,'Visible','off')

color_mat = [0*strain_100(:,3) 0*strain_100(:,3) 0*strain_100(:,3)];
color_mat(:,1) = ((strain_130(:,1))/(max(max(strain_130))));
color_mat(:,2)= 0;
color_mat(:,3) = 0.5;%1-((strain_130(:,1))/(max(max(strain_130))));
colormap(sort(color_mat,'ascend'));
h=colorbar('FontSize',40);
caxis([0 3]);
set(get(h,'label'),'string','Strain, \epsilon_\theta (%)',...
    'FontSize', 40);
end


saveas(gcf,'figures/shea_strain_autoreg.fig')
saveas(gcf,'figures/shea_strain_autoreg','epsc')
saveas(gcf,'figures/shea_strain_autoreg','jpeg')



figure('units','normalized','outerposition',[0 0 0.9 0.5])

subplot(1,2,2)
color_mat = [0*strain_100(:,1) 0*strain_100(:,1) 0*strain_100(:,1)];
color_mat(:,1) = log10(strain_130(:,1))/log10(max(max(strain_130)));
color_mat(:,2)= 0;
color_mat(:,3) = 0.5;

size_vec = strain_130(:,1)/max(max(strain_130))*20;
size_vec(size_vec==0)=min(size_vec(size_vec>0));

p1=plot(Ggraph,'k', 'Layout', 'layered', ...
    'EdgeColor', color_mat, ...
    'LineWidth',size_vec,...
    'ArrowSize', 0); hold on;
p2=plot(Ggraph,'k', 'Layout', 'layered', ...
    'ArrowSize', 12); hold off;
text(p2.XData(195),p2.YData(195)-1.5,'EA','fontsize',30);
text(p2.XData(1),p2.YData(1)+1.5,'AA','fontsize',30);
set(gca,'Visible','off')
colormap(sort(color_mat,'ascend'));
h=colorbar('FontSize',40);
caxis([0 3]);
set(get(h,'label'),'string','Strain, \epsilon_\theta (%)',...
    'FontSize', 40);

saveas(gcf,'figures/shea_strain_autoreg_for_colorbar.fig')
saveas(gcf,'figures/shea_strain_autoreg_for_colorbar','epsc')
saveas(gcf,'figures/shea_strain_autoreg_for_colorbar','jpeg')

   