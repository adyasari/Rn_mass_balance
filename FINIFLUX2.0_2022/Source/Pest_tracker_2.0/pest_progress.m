clear
path=import_string('FINIFLUX.xlsx','Tabelle1','C4:C4');  
str=strjoin(path);
slave='\slave1\Rn_modeled.csv';
file='\Rn_modeled.csv';
file2='\slave1\gw_inflow.dat';
path_l=strcat(str,slave);
path_c=strcat(str,file);
path_b=strcat(str,file2);

obs_conc=import_obs_conc('FINIFLUX.xlsx','Tabelle1',4,1000);
obs=import_obs('FINIFLUX.xlsx','Tabelle1',4,1000);
obs=obs(~isnan(obs));
obs_conc=obs_conc(~isnan(obs_conc));

delimiter = '\t';
formatSpec = '%*s%f%[^\n\r]';
delimiter2 = '';
formatSpec2 = '%f%[^\n\r]';

figure(1)
ax1=subplot(3,1,1);
plot(obs,obs_conc,'O')
ylabel(ax1,'222Rn (Bq/m3)')
xlabel(ax1,'Distance (m)')
Rt='Correlation Coefficient: ';
hold on


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

go = true;
i=1;
   
   
   while go
   

   
   if exist(path_l, 'file')
       s=dir(path_l);
       if s.bytes ~= 0 
           %copyfile(path_l, path_c);
           fileID = fopen(path_l,'r');
           dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
           fclose(fileID);
           simulated = dataArray{:, 1};
           R = corrcoef(obs_conc,simulated);
           RR=num2str(R(1,2));
           CC(i)=R(1,2);
           XX(i)=i;
           t=strcat(Rt,RR);
            
                      
           plot(obs,simulated,'-.');
           title(ax1,t);
           if i==1
           legend('observed','simulated');
           else
           end
            

           ax2=subplot(3,1,2);
           plot(XX,CC,'-X');
           ylabel(ax2,'Correlation Coefficient [-]')
           xlabel(ax2,'Iteration [-]')
           %refreshdata(hh);
           i=i+1;
           clear dataArray simulated R RR
           
           
       else
           
       end
   else
   end
   
    if exist(path_b, 'file')
       ss=dir(path_b);
       if ss.bytes ~= 0 
           %copyfile(path_l, path_c);
          
           fileIID = fopen(path_b,'r');
           dataArray2 = textscan(fileIID, formatSpec2, 'Delimiter', delimiter2,  'ReturnOnError', false);
           fclose(fileIID);
           simulated_gw = dataArray2{:, 1};
           reach=1:numel(simulated_gw);
          
           ax3=subplot(3,1,3);
           bar(reach,simulated_gw);
           ylabel(ax3,'Groundwater Inflow [m²/s]')
           xlabel(ax3,'Reach [-]')
           clear dataArray2 simulated_gw 
           
       else
       end
    else
    end
   
   
 
   
   pause(0.5);
    end



