close all;
clearvars;

path_data = '../data/';
path_output = '../processed/';

% model_name = 'test';
% domain_name = 'test_rectangle';
% ic_bacteria_name = 'line_rectangle';
% ic_fields_name = 'constant';
% parameter_set_name = 'test_rectangle';
model_name = 'test';
domain_name = 'test_ring';
ic_bacteria_name = 'line_ring';
ic_fields_name = 'constant';
parameter_set_name = 'test_ring';
run_nr = '0';

nr_particles = 1000;
dim = 2;
    
if strcmp(domain_name, 'test_rectangle')
    dimensions = [200, 100, 1];
    nr_points = [200, 100, 1];
    lengths = dimensions./nr_points;
end
if strcmp(domain_name, 'test_ring')
    dimensions = [100, 100, 1];
    nr_points = [200, 200, 1];
    lengths = dimensions./nr_points;
end

data_void=load(strcat(path_data,"Data_grid_void_",domain_name,"_",...
    parameter_set_name,".dat"));
data_solid=load(strcat(path_data,"Data_grid_solid_",domain_name,...
    "_",parameter_set_name,".dat"));

data_nutrient=load(strcat(path_data,"Data_nutrient_",model_name,"_",...
    domain_name,"_",ic_bacteria_name,"_",ic_fields_name,"_",...
    parameter_set_name,"_",run_nr,".dat"));
data_particle_state=load(strcat(path_data,"Data_particle_state_",...
    model_name,"_",domain_name,"_",ic_bacteria_name,"_",ic_fields_name,...
    "_",parameter_set_name,"_",run_nr,".dat"));

%%
close all;

f=figure(1);
c=parula();
nr_colors = size(c,1);
nr_times = size(data_particle_state,1);
time = data_particle_state(:,1);
positions=zeros(nr_times,nr_particles,dim);
orientations=zeros(nr_times,nr_particles);
states=zeros(nr_times,nr_particles);

for ii = 1:nr_times
    for dd = 1:dim
        positions(ii,:,dd) = data_particle_state(ii,(1+dd):(dim+2):end);
    end
    orientations(ii,:) = data_particle_state(ii,(2+dim):(dim+2):end);
    states(ii,:) = data_particle_state(ii,(3+dim):(dim+2):end);
end

plot(0,0,'visible','off');
hold on;
if strcmp(domain_name, 'test_rectangle')
    xlim([0-lengths(1),dimensions(1)+lengths(1)]);
    ylim([0-lengths(1),dimensions(2)+lengths(1)]);
end
if strcmp(domain_name, 'test_ring')
    xlim([-dimensions(1)/2-lengths(1),dimensions(1)/2+lengths(1)]);
    ylim([-dimensions(2)/2-lengths(2),dimensions(2)/2+lengths(2)]);
end
pbaspect(dimensions);
max_nutrient = max(data_nutrient(1,2:end));
writerObj = VideoWriter(strcat(path_output,"Video_",...
    model_name,"_",domain_name,"_",ic_bacteria_name,"_",ic_fields_name,...
    "_",parameter_set_name,"_",run_nr,".avi"));
writerObj.FrameRate = 10;
open(writerObj);
for i=1:nr_times
    time(i)
    color_idx=floor(data_nutrient(i,2:end)/max_nutrient*(nr_colors-1))+1;
    cla;
    scatter(data_void(:,1),data_void(:,2),20,c(color_idx,:),'filled','sq');
    hold on;
    scatter(data_solid(:,1),data_solid(:,2),20,[0 0 0],'filled','sq');
    color_idx=floor((data_nutrient(i,2:end))/max_nutrient*(nr_colors-1))+1;
    scatter(positions(i,:,1),positions(i,:,2),10,[1 0 0],'filled');
    frame=getframe(gcf);
    writeVideo(writerObj, frame);
end
close(writerObj);
