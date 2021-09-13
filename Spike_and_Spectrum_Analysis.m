%  review_run_git_hub 
% Created by Quinton Skilling and  Bolaji Eniwaye
% Plots Spike Rasters, calculates and plots Power Spectrum and Local Field Potential estimates


clear

%%


pre_str="C:/Users/eniwbola/Desktop/Research/Quinton_Proj/noiseSweep/noiseSweep/new_het_scan_stability_4_4_2021/review_run/cut_off_inh_34/7_21_2021/tau_300_real_mem_scan_34/"
spike_dir=strcat(pre_str,'spikeData/');
pre_dir_str=spike_dir;

weight_dir=strcat(pre_str,'Weights/');
conn_dir=strcat(pre_str,'Connectivity/');
spike_dir=strcat(pre_str,'spikeData/');
run_num=1;%3%0%0%3;
my_xlim=[0 500]+2000;

num_offset=0;



%%


A=dlmread(strcat(pre_dir_str,'1750_T1NL_',num2str(run_num),'.dat'));
A_exc=A(A(:,1)>0,:);
A_inh=A(A(:,1)<0,:);


B=dlmread(strcat(pre_dir_str,'1750_T2L_',num2str(run_num),'.dat'));
B_exc=B(B(:,1)>0,:);
B_inh=B(B(:,1)<0,:);

C=dlmread(strcat(pre_dir_str,'1750_T1AL_',num2str(run_num),'.dat'));
C_exc=C(C(:,1)>0,:);
C_inh=C(C(:,1)<0,:);

rast_total=[A_exc;B_exc;C_exc];
rast_total_inh=[A_inh;B_inh;C_inh];

 [sigFreqr, fullSigr,find_times,lfp1_new]    =calculateSignalFrequency(rast_total, 3300)
 [sigFreq2r, fullSig2r,find_times_2,lfp2_new]=calculateSignalFrequency(rast_total_inh, 3300)
 


%%

A=dlmread(strcat(pre_dir_str,'1750_T1NL_',num2str(run_num),'.dat'));

A_exc=A(A(:,1)>0,:);
A_inh=A(A(:,1)<0,:);


un_exc=unique(A_exc(:,4));
un_inh=unique(A_inh(:,4));



[val_exc rank_exc]=sort(un_exc);
[val_inh rank_inh]=sort(un_inh);
A_rank_exc=A_exc(rank_exc,:);
A_rank_inh=A_inh(rank_inh,:);

A_exc_sort=A_exc;
A_inh_sort=A_inh;
corresponding_rank_id_exc=[];
for ii =1:length(val_exc)
    %so find
    un_asort=unique(A_exc_sort(find(A_exc(:,4)==val_exc(ii)),1));
    corresponding_rank_id_exc=[corresponding_rank_id_exc; un_asort(1) rank_exc(ii)];
   A_exc_sort(find(A_exc(:,4)==val_exc(ii)),1)=rank_exc(ii);
   
end
for ii =1:length(val_inh)
   A_inh_sort(find(A_inh(:,4)==val_inh(ii)),1)=rank_inh(ii);
    
end

    
   
%---------
B=dlmread(strcat(pre_dir_str,'1750_T2L_',num2str(run_num),'.dat'));

B_exc=B(B(:,1)>0,:);
B_inh=B(B(:,1)<0,:);

un_exc=unique(B_exc(:,4));
un_inh=unique(B_inh(:,4));

[val_exc rank_exc]=sort(un_exc);
[val_inh rank_inh]=sort(un_inh);
B_rank_exc=B_exc(rank_exc,:);
B_rank_inh=B_inh(rank_inh,:);

B_exc_sort=B_exc;
B_inh_sort=B_inh;
corresponding_rank_id_exc=[]
for ii =1:length(val_exc)
    %so find
    try
    un_asort=unique(B_exc_sort(find(B_exc(:,4)==val_exc(ii)),1));
    corresponding_rank_id_exc=[corresponding_rank_id_exc; un_asort(1) rank_exc(ii)];
   B_exc_sort(find(B_exc(:,4)==val_exc(ii)),1)=rank_exc(ii);
    catch
       continue 
    end
   
end
for ii =1:length(val_inh)
   B_inh_sort(find(B_inh(:,4)==val_inh(ii)),1)=rank_inh(ii);
    
end


C=dlmread(strcat(pre_dir_str,'1750_T1AL_',num2str(run_num),'.dat'));

C_exc=C(C(:,1)>0,:);
C_inh=C(C(:,1)<0,:);


un_exc=unique(C_exc(:,4));
un_inh=unique(C_inh(:,4));

[val_exc rank_exc]=sort(un_exc);
[val_inh rank_inh]=sort(un_inh);
C_rank_exc=C_exc(rank_exc,:);
C_rank_inh=C_inh(rank_inh,:);

C_exc_sort=C_exc;
C_inh_sort=C_inh;
corresponding_rank_id_exc=[]
for ii =1:length(val_exc)
    %so find
    un_asort=unique(C_exc_sort(find(C_exc(:,4)==val_exc(ii)),1));
    corresponding_rank_id_exc=[corresponding_rank_id_exc; un_asort(1) rank_exc(ii)];
   C_exc_sort(find(C_exc(:,4)==val_exc(ii)),1)=rank_exc(ii);
   
end
for ii =1:length(val_inh)
   C_inh_sort(find(C_inh(:,4)==val_inh(ii)),1)=rank_inh(ii);
    
end


    
%%
rast_total=[A;B;C];
rast_total=[A_exc;B_exc;C_exc];

stab_before= calculateStability(rast_total(:,1:3), 2000, 3300);
stab_after=calculateStability(rast_total(:,1:3), 2000, 3500)

[sigFreq, fullSig]=calculateSignalFrequency(rast_total, 3300)
[sigFreq2, fullSig2]=calculateSignalFrequency(rast_total, 3300)

stab_diff=stab_after-stab_before
%%


%% new

figure(6)   
clf
plot(fullSig(:,1),fullSig(:,2),'b-')
hold on
plot(fullSig2(:,1),fullSig2(:,2),'r-')
hold off

xlabel('Frequency (hz)')
ylabel('Power Spectral density')

 
 


%%

figure(718)
plot(A_exc_sort(:,3), A_exc_sort(:,1)+(800-max(A_exc_sort(:,1))) ,'k.')
hold on
plot(B_exc_sort(:,3), B_exc_sort(:,1)+(800-max(B_exc_sort(:,1))),'k.');plot(C_exc_sort(:,3), C_exc_sort(:,1)+num_offset+(800-max(C_exc_sort(:,1)+num_offset)),'k.')
plot(A_inh_sort(:,3),A_inh_sort(:,1)+800,'b.');plot(B_inh_sort(:,3),B_inh_sort(:,1)+800,'b.'); plot(C_inh_sort(:,3),C_inh_sort(:,1)+800,'b.');
title('Type 2 With Implanted Memory (High DC)','fontsize',15)
ylabel('Cell ID','fontsize',22)
set(gca,'xtick',[]);set(gca,'xticklabel',[]);set(gca,'ytick',[]);set(gca,'yticklabel',[])
ylim([0 1000])
yline(800,'r--','linewidth',5)
hold off
xlim(my_xlim+600)
set(gcf,'color','white')
ae=gcf;

figure(719)

a6=plot([1:length(lfp1_new)]*.05,-lfp1_new*5+300-150,'k-','linewidth',2)
 hold on
 a7=plot([1:length(lfp2_new)]*.05,-lfp2_new*5+300-150,'b-','linewidth',2)
 hold off

my_xlim=[100  600]
xlim(my_xlim+600)

ylabel('LFP','fontsize',22)
set(gcf,'color','white')
 ylim([-150 450]+50)
set(gca,'ytick',[])
set(gca,'yticklabel',[])

x0=360; y0=246; width=560; height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
 xticks([600 850 1100 ])
 xticklabels({'5000','5250','5500'})
af=gcf; ax = gca; ax.FontSize = 16; 
saveas(gcf,'C:\Users\eniwbola\Desktop\Research\Quinton_Proj\Figures\Review_Figs\Final\LFP_with_mem.epsc')


