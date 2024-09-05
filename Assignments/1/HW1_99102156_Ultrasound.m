%% Q1
clc;close all;



%********1D Environment***********
%defining environment
x_len = 10; %total x length [mm]
dx = 0.05;

Nx = length(0:dx:x_len);
kgrid1 = kWaveGrid(Nx, dx);

%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_r = round(0.1/dx);    % [grid points]
disc1 = zeros(Nx,1);
disc1(disc_x-disc_r:disc_x+disc_r) = disc_amp;
source1.p0 = disc1;



%defining reciever sensor, a centered circular sensor
sensor_r = 4;   % [mm]
num_sensor_points = 50;
sensor1.mask = linspace(-sensor_r,sensor_r,num_sensor_points);


%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

bone_speed = 1540;  % [m/s]
bone_density = 1500; %[kg/m^3]

air_speed = 340;  % [m/s]
air_density = 1.2; %[kg/m^3]

medium.sound_speed = water_speed;
medium.density = water_density;

% % run the simulation
sensor_data1 = kspaceFirstOrder1D(kgrid1, medium, source1, sensor1);








% 
% %********2D Environment***********
%defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
dx = 0.05;
dy = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
kgrid2 = kWaveGrid(Nx, dx, Ny, dy);

%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_y = round(1/dx) + 1;    % [grid points]
disc_r = round(0.1/dx);    % [grid points]

disc2 = disc_amp * makeDisc(Nx, Ny, disc_x, disc_y, disc_r);

source2.p0 = disc2;



%defining reciever sensor, a centered circular sensor
sensor_r = 4;   % [mm]
num_sensor_points = 50;
sensor2.mask = makeCartCircle(sensor_r, num_sensor_points); %setting points at which sensor elements are placed.



% run the simulation
sensor_data2 = kspaceFirstOrder2D(kgrid2, medium, source2, sensor2);









% % 
% % %********3D Environment***********
% %defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
z_len = 10; %total z length [mm]

dx = 0.05;
dy = 0.05;
dz = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
Nz = length(0:dz:z_len);

kgrid3 = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% 
% %defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_y = round(1/dx) + 1;    % [grid points]
disc_z = round(5/dx) + 1;   % [grid points]
disc_r = round(0.1/dx);    % [grid points]

disc3 = disc_amp * makeSphere(Nx, Ny, Nz, disc_r);

source3.p0 = disc3;
% 
% 
% 
% %defining reciever sensor, a centered circular sensor
sensor_r = 4;   % [mm]
num_sensor_points = 50;
sensor3.mask = makeCartSphere(sensor_r, num_sensor_points); %setting points at which sensor elements are placed.




%NOTE: PLEASE UNCOMMENT THE FOLLOWING LINE TO PLOT FOR 3D SITUATION (My computer couldn't handle it so I commented it, but it's probably correct!)

% sensor_data3 = kspaceFirstOrder3D(kgrid3, water_medium, source3, sensor3,'PlotSim',false,'PlotPML',false);








% plot the simulated sensor data
figure;

subplot(131);
imagesc(sensor_data1, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('1D simulation');

subplot(132);
imagesc(sensor_data2, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('2D simulation');

if (exist('sensor_data3'))
    subplot(133);
    imagesc(sensor_data3, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    colorbar;
    title('3D simulation');
end



%% Q2
clc; close all;
%defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
dx = 0.05;
dy = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_y = round(1/dx) + 1;    % [grid points]
disc_r = round(0.1/dx);    % [grid points]

disc = disc_amp * makeDisc(Nx, Ny, disc_x, disc_y, disc_r);

source.p0 = disc;


sensor.mask = [0 0; -4 4];




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

bone_speed = 1540;  % [m/s]
bone_density = 1500; %[kg/m^3]

air_speed = 340;  % [m/s]
air_density = 1.2; %[kg/m^3]


medium.sound_speed = bone_speed * ones(Nx, Ny);   % [m/s]
medium.sound_speed(:,1:floor(Ny/2)) = water_speed;       % [m/s]
medium.density = bone_density * ones(Nx, Ny);       % [kg/m^3]
medium.density(:,1:floor(Ny/2)) = water_density;          % [kg/m^3]


sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);


figure;

% imagesc(sensor_data, [-1, 1]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step');
% colorbar;
% title('Water-Bone Simulation');

plot(sensor_data(1,:));
hold on;
plot(sensor_data(2,:)+5);
yline(5,'--');
yline(0,'--');
ylabel('Sensor Position & Amplitude');
ylim([-10,22]);
xlabel('Time Step');
title('Water-Bone Simulation');
legend('sensor1 (5mm,1mm)','sensor2 (5mm,9mm)');




%% Q3
clc; close all;
%defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
dx = 0.05;
dy = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_y = round(1/dx) + 1;    % [grid points]
disc_r = round(0.1/dx);    % [grid points]

disc = disc_amp * makeDisc(Nx, Ny, disc_x, disc_y, disc_r);

source.p0 = disc;


%defining reciever sensor, a centered circular sensor
sensor_r = 4;   % [mm]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_r, num_sensor_points); %setting points at which sensor elements are placed.




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

bone_speed = 1540;  % [m/s]
bone_density = 1500; %[kg/m^3]

air_speed = 340;  % [m/s]
air_density = 1.2; %[kg/m^3]

medium.sound_speed = air_speed * ones(Nx, Ny);   % [m/s]
medium.sound_speed(:,1:floor(0.2*Ny)) = bone_speed;       % [m/s]
medium.sound_speed(:,floor(0.7*Ny):end) = water_speed;       % [m/s]

medium.density = air_density * ones(Nx, Ny);       % [kg/m^3]
medium.density(:,1:floor(0.2*Ny)) = bone_density;          % [kg/m^3]
medium.density(:,floor(0.7*Ny):end) = water_density;          % [kg/m^3]


sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);


figure;

imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('Water-Bone Simulation');

% plot(sensor_data(1,:));
% hold on;
% plot(sensor_data(2,:)+5);
% yline(5,'--');
% yline(0,'--');
% ylabel('Sensor Position & Amplitude');
% ylim([-10,22]);
% xlabel('Time Step');
% title('Water-Bone Simulation');
% legend('sensor1 (5mm,1mm)','sensor2 (5mm,9mm)');



%% Q4
clc; close all;
%defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
dx = 0.05;
dy = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x = round(5/dx) + 1;    % [grid points]
disc_y = round(1/dx) + 1;    % [grid points]
disc_r = round(0.1/dx);    % [grid points]

disc = disc_amp * makeDisc(Nx, Ny, disc_x, disc_y, disc_r);

source.p0 = disc;


%defining reciever sensor, a centered circular sensor
sensor_r = 4;   % [mm]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_r, num_sensor_points); %setting points at which sensor elements are placed.




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

bone_speed = 1540;  % [m/s]
bone_density = 1500; %[kg/m^3]

air_speed = 340;  % [m/s]
air_density = 1.2; %[kg/m^3]

medium.sound_speed = air_speed * ones(Nx, Ny);   % [m/s]
medium.sound_speed(:,1:floor(0.2*Ny)) = bone_speed;       % [m/s]
medium.sound_speed(:,floor(0.8*Ny):end) = bone_speed;       % [m/s]

medium.density = air_density * ones(Nx, Ny);       % [kg/m^3]
medium.density(:,1:floor(0.2*Ny)) = bone_density;          % [kg/m^3]
medium.density(:,floor(0.8*Ny):end) = bone_density;          % [kg/m^3]


sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);


figure;

imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('Bone-Air-Bone Simulation');


%% Q5
clc; close all;
%defining environment
x_len = 10; %total x length [mm]
y_len = 10; %total y length [mm]
dx = 0.05;
dy = 0.05;

Nx = length(0:dx:x_len);
Ny = length(0:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%defining source disc
disc_amp = 20; %disk initial pressure amplitude [pa]
disc_x_mm = 5; %source center [mm]
disc_y_mm = 1; %source center [mm]
disc_x = round(disc_x_mm/dx) + 1;    % source center [grid points]
disc_y = round(disc_y_mm/dx) + 1;    % source center [grid points]
disc_r = round(0.1/dx);    % source center [grid points]

disc = disc_amp * makeDisc(Nx, Ny, disc_x, disc_y, disc_r);

source.p0 = disc;


%defining reciever sensor, a centered circular sensor
sensors_dist = 0.1;   % sensors distance from source & eachother [mm]
num_sensor_points = 50;
sensor.mask = [disc_x_mm * ones(1,num_sensor_points) - x_len/2 ; [disc_y_mm + sensors_dist: sensors_dist: disc_y_mm + num_sensor_points*sensors_dist] - y_len/2];




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

bone_speed = 1540;  % [m/s]
bone_density = 1500; %[kg/m^3]

air_speed = 340;  % [m/s]
air_density = 1.2; %[kg/m^3]

medium.sound_speed = water_speed;   % [m/s]
medium.density = water_density;       % [kg/m^3]


sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

figure;

imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('Water Simulation, linear shape sensors');


sensorsPeak_data = max(abs(sensor_data.')); %finding peaks in sensors data
dist_from_source = sensors_dist: sensors_dist: num_sensor_points*sensors_dist; %sensors' Distance from source
figure;

plot(dist_from_source, sensorsPeak_data,'LineWidth',1.5);
ylabel('Peak value');
xlabel('Sensor distance from source (mm)');
title('Sensors Peak vs their distance from source');
