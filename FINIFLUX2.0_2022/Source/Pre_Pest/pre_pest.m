path=import_string('FINIFLUX.xlsx','Tabelle1','C4:C4');
CPU=import_float('FINIFLUX.xlsx','Tabelle1',6,6);
degasing_mode=import_float('FINIFLUX.xlsx','Tabelle1',8,8);
HZ_mode=import_float('FINIFLUX.xlsx','Tabelle1',13,13);
RT_mode=import_float('FINIFLUX.xlsx','Tabelle1',17,17);
GAMMA_alpha=import_float('FINIFLUX.xlsx','Tabelle1',22,22);
upstream_weigth=import_float('FINIFLUX.xlsx','Tabelle1',24,24);

obs=import_obs('FINIFLUX.xlsx','Tabelle1',4,1000);
obs=obs(~isnan(obs));

obs_conc=import_obs_conc('FINIFLUX.xlsx','Tabelle1',4,1000);
obs_conc=obs_conc(~isnan(obs_conc));

data=importfile_matrix('FINIFLUX.xlsx','Tabelle1','J4:T1000');
data=data(~isnan(data));

n_obs=numel(obs);

river_depth=data(n_obs:2*n_obs-2);
river_width=data(2*n_obs-1:3*n_obs-3);
RN_GW=data(3*n_obs-2:4*n_obs-4);
depth_hz=data(4*n_obs-3:5*n_obs-5);
por_hz=data(5*n_obs-4:6*n_obs-6);
surf_inflow=data(6*n_obs-5:7*n_obs-7);
inflow_length=data(7*n_obs-6:8*n_obs-8);
inflow_conc=data(8*n_obs-7:9*n_obs-9);
discharge=data(9*n_obs-8:10*n_obs-10);
k_values_def=data(10*n_obs-9:11*n_obs-11);


dlmwrite('processor_nodes.txt',CPU,'delimiter','\t','precision', 16);
dlmwrite('obs_concentrations.txt',obs_conc,'delimiter','\t','precision', 16);
dlmwrite('n_observation.dat',n_obs,'delimiter','\t','precision', 16);
dlmwrite('observation_points.txt',obs,'delimiter','\t','precision', 16);
dlmwrite('river_depth.txt',river_depth,'delimiter','\t','precision', 16);
dlmwrite('river_width.txt',river_width,'delimiter','\t','precision', 16);
dlmwrite('gw_concentration.txt',RN_GW,'delimiter','\t','precision', 16);
dlmwrite('discharge.dat',discharge,'delimiter','\t','precision', 16);
dlmwrite('initial_c.dat',obs_conc(1),'delimiter','\t','precision', 16);
dlmwrite('hz_on_off.dat',HZ_mode,'delimiter','\t','precision', 16);
dlmwrite('k_degas.dat',degasing_mode,'delimiter','\t','precision', 16);
dlmwrite('distribution.dat',RT_mode,'delimiter','\t','precision', 16);
dlmwrite('river_inflow.dat',surf_inflow,'delimiter','\t','precision', 16);
dlmwrite('inflow_conc.dat',inflow_conc,'delimiter','\t','precision', 16);
dlmwrite('upstream_weight.dat',upstream_weigth,'delimiter','\t','precision', 16);
dlmwrite('gamma_alpha.dat',GAMMA_alpha,'delimiter','\t','precision', 16);
dlmwrite('k_values.dat',k_values_def,'delimiter','\t','precision', 16);
dlmwrite('depth_hz.dat',depth_hz,'delimiter','\t','precision', 16);
dlmwrite('porosity_hz.dat',por_hz,'delimiter','\t','precision', 16);
dlmwrite('inflow_length.dat',inflow_length,'delimiter','\t','precision', 16);




str=strjoin(path);
fid = fopen('path_definition.txt','wt');
fprintf(fid,'%s',str);
fclose(fid);

system('write_pest.exe &');


