%  Stability_vs_Noise_git  
% Created by Quinton Skilling and  Bolaji Eniwaye
% loads and creates stability vs noise plot and stability vs memory size


clearvars
close all

save_dir='C:/Users/eniwbola/Desktop/Research/Quinton_Proj/Figures/'


load('stability_freq.mat')
z=het_val;
x=stab_val%het_end;
y=spec_dens%freq_dens_after;
yourValues=z;
figure(185)
clf
numValues = length(z);
cmap = jet(numValues);
v = rescale(z, 1, numValues); 


markerColors = zeros(numValues, 3)
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    %markerColors(k, :) = cmap(row, :);
      markerColors(k, :) = cmap(k, :);
end
% Create the scatter plot.
figure(14)
clf
scatter(y, x-min(x), [], markerColors,'filled');
colorbar

caxis([0 400])
colormap('jet')
ylim([-.2 1]);  xlim([0 .4]);
xlabel('Change in Spectral Power','fontsize',14); ylabel('Change in Stability','fontsize',14);



hcb = colorbar;
set(hcb, 'Ytick', [0 200 400]); %// 4 yticks, each "in the middle" of one color
set(hcb, 'YTickLabel', {'0.0' ,     '0.25' ,    '0.5'})

caxis([0 400])
ax = gca; box(ax,'on') ;add = gca;
extra = axes('Position',get(add,'Position'),'Color','none','XTick',[],'YAxisLocation','right');

xticks([0.02 .2 .38]); xticklabels({'0','0.2','0.4 '});
yticks([0 1]); yticklabels({'',' '});

ed=ylabel('Memory Size   (Network Fraction) ','fontsize',14)

ed.Position=ed.Position+[.14 0 0]


%%
%so i need to analyze the stability before and after
% I need to put the title lines so I can switch them out
% I need to plot the before and after weight

load('stability_nz.mat')
figure(1138)


h_1= errorbar(noise_stab_hz,stab_before_nz_run,stab_bef_err,'ko','MarkerFaceColor','k');
hold on

h_2 = errorbar(noise_stab_hz,stab_after_nz_run,stab_aft_err,'bo','MarkerFaceColor','b');
set(get(h_1,'Parent'), 'XScale', 'log')
set(get(h_2,'Parent'), 'XScale', 'log')


hold off
xlabel('Noise Frequency [Hz]','fontsize',14)
ylabel('Stability','fontsize',14)
legend('Post-Learning','Pre-Learning')

ylim([-.1 .9])
set(gcf,'color','w');




%%


pre_str='C:/Users/eniwbola/Desktop/Research/Quinton_Proj/noiseSweep/noiseSweep/new_het_scan_stability_4_4_2021/noise_scan_2_longer/gks_1p5/'
run_num=4;%3%0%0%3;


weight_dir=strcat(pre_str,'Weights/');
conn_dir=strcat(pre_str,'Connectivity/');
spike_dir=strcat(pre_str,'spikeData/');


%%
pre_dir_str='spikeData/';
pre_dir_str=spike_dir;
A=dlmread(strcat(pre_dir_str,'1850_T1NL_',num2str(run_num),'.dat'));
A_exc=A(A(:,1)>0,:);
A_inh=A(A(:,1)<0,:);


B=dlmread(strcat(pre_dir_str,'1850_T2L_',num2str(run_num),'.dat'));
B_exc=B(B(:,1)>0,:);
B_inh=B(B(:,1)<0,:);

C=dlmread(strcat(pre_dir_str,'1850_T1AL_',num2str(run_num),'.dat'));
C_exc=C(C(:,1)>0,:);
C_inh=C(C(:,1)<0,:);



rast_total=[A;B;C];
stab_before= calculateStability(rast_total(:,1:3), 8000, 10000);
% 
[sigFreq, fullSig,find_times,lfp1]=calculateSignalFrequency(rast_total, 5000)
[sigFreq2, fullSig2,find_times_2,lfp2]=calculateSignalFrequency(rast_total, 10000)
figure(61)
clf
plot(fullSig(:,1),fullSig(:,2),'b-')
hold on
plot(fullSig2(:,1),fullSig2(:,2),'r-')
hold off

%%


A=dlmread(strcat(pre_dir_str,'1850_T1NL_',num2str(run_num),'.dat'));

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

B=dlmread(strcat(pre_dir_str,'1850_T2L_',num2str(run_num),'.dat'));

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

C=dlmread(strcat(pre_dir_str,'1850_T1AL_',num2str(run_num),'.dat'));


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
stab_before= calculateStability(rast_total(:,1:3), 5900, 9900);
stab_after=calculateStability(rast_total(:,1:3), 15000, 19000)
stab_diff=stab_after-stab_before

C_inh_sort

%--------------

[sigFreq, fullSig,find_times,lfp1]=calculateSignalFrequency(rast_total, 5000)
[sigFreq2, fullSig2,find_times,lfp2]=calculateSignalFrequency(rast_total, 20000)




