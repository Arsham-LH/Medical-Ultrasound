%% Q1, Q5, Q6.
clc; close all;



%defining environment (grid)
x_len = 22e-3; %total x length [m]
y_len = 22e-3; %total y length [m]
dx = 0.05e-3;
dy = 0.05e-3;

Nx = length(dx:dx:x_len);
Ny = length(dy:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);



% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny);   % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);         % [m/s]
medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000;     % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 800;      % [m/s]
medium.density(Nx/2:end, :)                 = 1200;     % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1;   % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.5;   % [dB/(MHz^2 cm)]

center_x = Nx/2;
center_y = Ny/2;

% create initial pressure distribution using makeDisc
disc_magnitude = 25; % [Pa]
disc_x_pos_grid = center_x - 50;    % [grid points]
disc_y_pos_grid = center_y - 175;    % [grid points]
% disc_radius = 5;    % [grid points]
source.p0 = zeros(Nx,Ny);
source.p0(disc_x_pos_grid,disc_y_pos_grid) = disc_magnitude;


disc_x_pos = (disc_x_pos_grid-center_x)/(Nx/x_len);
disc_y_pos = (disc_y_pos_grid-center_y)/(Ny/y_len);

alpha = linspace(pi/12,pi/3.2,50);
arrowsNum = length(alpha);

beta = asin(sin(alpha) .* medium.sound_speed_compression(end,end) ./ medium.sound_speed_compression(1,1));
removed_ind = find(imag(beta) ~= 0);
beta (removed_ind) = [];
tmp_alpha = alpha;
tmp_alpha(removed_ind) = [];

gamma = asin(sin(alpha) * medium.sound_speed_shear(end,end) / medium.sound_speed_compression(1,1));

%defining reciever sensor, a centered circular sensor
imping_sensors_dist = 1e-3; %impinging sensors distance from the source

imping_sensors = [disc_x_pos + imping_sensors_dist*cos(alpha); disc_y_pos + imping_sensors_dist*sin(alpha)];  
reflect_sensors = [disc_x_pos * ones(1,arrowsNum); disc_y_pos + 2 * tan(alpha) * abs(disc_x_pos)];
transmit_compression_sensors = [-disc_x_pos * ones(1,length(beta)); disc_y_pos + tan(tmp_alpha) * abs(disc_x_pos) + tan(beta) * abs(disc_x_pos)];
transmit_shear_sensors = [-disc_x_pos * ones(1,arrowsNum); disc_y_pos + tan(alpha) * abs(disc_x_pos) + tan(gamma) * abs(disc_x_pos)];


sensor.mask = [imping_sensors, reflect_sensors, transmit_compression_sensors, transmit_shear_sensors];





% create the time array
% t_end = 7e-6;   % [s]
kgrid.makeTime(medium.sound_speed_compression(:));




% define input arguments
input_args = {'PlotPML', false,...
     'DataCast', 'single'};

% run the simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});


sensors_r = radius(disc_x_pos, disc_y_pos, sensor.mask(1,:), sensor.mask(2,:)); %all sensors radius

figure;

imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('2D simulation');

figure;

slct_reflect_data = sensor_data(arrowsNum+1:2*arrowsNum,:);
reflect_square_ratio_compen = (max(sensor_data(arrowsNum+1:2*arrowsNum,440:end).') ./ max(sensor_data(1:arrowsNum,:).')).^2 .*(sensors_r(arrowsNum+1:2*arrowsNum) ./ sensors_r(1:arrowsNum));
subplot(131);
plot(alpha*180/pi,reflect_square_ratio_compen);
title('Reflection square ratio');
xlabel('Alpha (deg)');
ylabel('Ratio');

tmp_sensors_r = sensors_r(1:arrowsNum);
tmp_sensors_r(removed_ind) = [];
tmp_imping_data = sensor_data(1:arrowsNum,:);
tmp_imping_data(removed_ind,:) = [];
transmit_compression_square_ratio_compen = (max(sensor_data(2*arrowsNum+1:2*arrowsNum + length(beta),:).') ./ max(tmp_imping_data.')).^2 .*(sensors_r(2*arrowsNum+1:2*arrowsNum + length(beta)) ./ tmp_sensors_r);
subplot(132);
plot(tmp_alpha*180/pi, transmit_compression_square_ratio_compen);
title('Transmission (compression) square ratio');
xlabel('Alpha (deg)');
ylabel('Ratio');


transmit_shear_square_ratio_compen = (max(sensor_data(2*arrowsNum + length(beta) + 1 : end,:).') ./ max(sensor_data(1:arrowsNum,:).')).^2 .*(sensors_r(2*arrowsNum + length(beta) + 1 : end) ./ sensors_r(1:arrowsNum));
subplot(133);
plot(alpha*180/pi,transmit_shear_square_ratio_compen);
title('Transmission (shear) square ratio');
xlabel('Alpha (deg)');
ylabel('Ratio');




%% Q2: Scattering wave
%clear;
clc; close all;

%defining environment (grid)
x_len = 10e-3; %total x length [mm]
y_len = 10e-3; %total y length [mm]
dx = 0.05e-3;
dy = 0.05e-3;

Nx = length(dx:dx:x_len);
Ny = length(dy:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);



% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny);   % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);         % [m/s]
medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000;     % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 800;      % [m/s]
medium.density(Nx/2:end, :)                 = 1200;     % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.0;   % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.0;   % [dB/(MHz^2 cm)]


% create initial pressure distribution using makeDisc
disc_magnitude = 25; % [Pa]
disc_x_pos_grid = 50;    % [grid points]
disc_y_pos_grid = 50;    % [grid points]


% disc_radius = 5;    % [grid points]
source.p0 = zeros(Nx,Ny);
source.p0(disc_x_pos_grid,disc_y_pos_grid) = disc_magnitude;

center_x = Nx/2;
center_y = Ny/2;


disc_x_pos = (disc_x_pos_grid-center_x)/(Nx/x_len);
disc_y_pos = (disc_y_pos_grid-center_y)/(Ny/y_len);



%defining reciever sensor, a centered circular sensor
sensor_r = 2e-3;   % [mm]
num_sensor_points = 10;
test_sensor.mask =[disc_x_pos * ones(1,num_sensor_points) ;linspace(-sensor_r,sensor_r,num_sensor_points)];

sensors_r = radius(disc_x_pos, disc_y_pos, test_sensor.mask(1,:), test_sensor.mask(2,:));

kgrid.makeTime(medium.sound_speed_compression(:));

% define input arguments
input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
     'DataCast', 'single'};

% run the simulation
test_sensor_data = pstdElastic2D(kgrid, medium, source, test_sensor, input_args{:});

test_sensor_data_compen = test_sensor_data .* sqrt((sensors_r.'));
figure;

subplot(131);
imagesc(test_sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('2D simulation');

subplot(132);
plot(test_sensor.mask(2,:) - disc_y_pos, max(test_sensor_data.'));

title('pressure reduction vs sensors y position');
xlabel('y position (mm)');
ylabel('pressure (Pa)');

subplot(133);
plot(test_sensor.mask(2,:) - disc_y_pos, max(test_sensor_data_compen.'));

title('pressure reduction (compensated) vs sensors y position');
xlabel('y position (mm)');
ylabel('pressure (Pa)');
ylim([0,0.1]);

figure;
x = linspace(0.5,4.5,num_sensor_points);
plot(x,1.88 ./ sqrt(x));
title('f(x) = 1.88/sqrt(x)');
xlabel('x(mm)');
ylabel('f(x)');


%% Q3: absorption
clc; close all;

%defining environment (grid)
x_len = 10e-3; %total x length [mm]
y_len = 10e-3; %total y length [mm]
dx = 0.05e-3;
dy = 0.05e-3;

Nx = length(dx:dx:x_len);
Ny = length(dy:dy:y_len);
kgrid = kWaveGrid(Nx, dx, Ny, dy);



% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny);   % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);         % [m/s]
medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000;     % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 800;      % [m/s]
medium.density(Nx/2:end, :)                 = 1200;     % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1;   % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.5;   % [dB/(MHz^2 cm)]


% create initial pressure distribution using makeDisc
disc_magnitude = 25; % [Pa]
disc_x_pos_grid = 50;    % [grid points]
disc_y_pos_grid = 50;    % [grid points]


% disc_radius = 5;    % [grid points]
source.p0 = zeros(Nx,Ny);
source.p0(disc_x_pos_grid,disc_y_pos_grid) = disc_magnitude;

center_x = Nx/2;
center_y = Ny/2;


disc_x_pos = (disc_x_pos_grid-center_x)/(Nx/x_len);
disc_y_pos = (disc_y_pos_grid-center_y)/(Ny/y_len);



%defining reciever sensor, a centered circular sensor
sensor_r = 2e-3;   % [mm]
num_sensor_points = 10;
test_sensor.mask =[disc_x_pos * ones(1,num_sensor_points) ;linspace(-sensor_r,sensor_r,num_sensor_points)];

sensors_r = radius(disc_x_pos, disc_y_pos, test_sensor.mask(1,:), test_sensor.mask(2,:));

kgrid.makeTime(medium.sound_speed_compression(:));

% define input arguments
input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
     'DataCast', 'single'};

% run the simulation
test_sensor_data_absorp = pstdElastic2D(kgrid, medium, source, test_sensor, input_args{:});

test_sensor_data_compen_absorp = test_sensor_data_absorp .* sqrt((sensors_r.'));
figure;

subplot(131);
imagesc(test_sensor_data_absorp, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
title('2D simulation');

subplot(132);
plot(test_sensor.mask(2,:) - disc_y_pos, max(test_sensor_data_absorp.'));

title('pressure reduction vs sensors y position');
xlabel('y position (mm)');
ylabel('pressure (Pa)');

subplot(133);
plot(test_sensor.mask(2,:) - disc_y_pos, max(test_sensor_data_compen_absorp.'));

title('pressure reduction (compensated) vs sensors y position');
xlabel('y position (mm)');
ylabel('pressure (Pa)');


figure;

plot(test_sensor.mask(2,:) - disc_y_pos, max(test_sensor_data_compen_absorp.') ./ max(test_sensor_data_compen.'));
title('absorped/non-absorped Pressure ratio (compensated scattering) vs sensors y position');
xlabel('y position (mm)');
ylabel('Ratio');


%% Q4: computing reflect. and transmit. ratio theoritically
clc; close all;

alpha = linspace(0,pi/2,100);

beta = asin(sin(alpha) * medium.sound_speed_compression(end,end) / medium.sound_speed_compression(1,1));
gamma = asin(sin(alpha) * medium.sound_speed_shear(end,end) / medium.sound_speed_compression(1,1));


% computing reflection and transmission coeffs
A = sin(gamma) .* sin(2*gamma) .* (cos(gamma) - medium.sound_speed_shear(end,end) / medium.sound_speed_compression(end,end) .* cos(beta));
B = cos(2*gamma).^2;
G = (medium.density(end, end) * medium.sound_speed_compression(end,end) .* cos(alpha)) ./ (medium.density(1, 1) * medium.sound_speed_compression(1,1) .* cos(beta));
D = (medium.sound_speed_shear(end,end) / medium.sound_speed_compression(1,1))^2 .* sin(2*alpha) .* sin(2*gamma);

reflect_square_ratio = ((1-G.*(1-2*A))./(1+G.*(1-2*A))).^2;
transmit_shear_square_ratio = 4*B.*G./((1+G.*(1-2*A)).^2);
transmit_compression_square_ratio = 4 * medium.density(end, end)/medium.density(1, 1) .* D ./ ((1+G.*(1-2*A)).^2);


figure;

subplot(131);
hold on;
plot(alpha*180/pi,real(reflect_square_ratio),'Linewidth',1.5);
plot(alpha*180/pi,imag(reflect_square_ratio));
plot(alpha*180/pi,abs(reflect_square_ratio));
title('square reflection ratio vs alpha');
xlabel('alpha(deg)');
ylabel('Ratio');
legend('real part', 'imaginary part', 'abs');

subplot(132);
hold on;
plot(alpha*180/pi,real(transmit_compression_square_ratio),'Linewidth',1.5);
plot(alpha*180/pi,imag(transmit_compression_square_ratio));
plot(alpha*180/pi,abs(transmit_compression_square_ratio));
title('square longtudinal ratio vs alpha');
xlabel('alpha(deg)');
ylabel('Ratio');
legend('real part', 'imaginary part', 'abs');

subplot(133);
hold on;
plot(alpha*180/pi,real(transmit_shear_square_ratio),'Linewidth',1.5);
plot(alpha*180/pi,imag(transmit_shear_square_ratio));
plot(alpha*180/pi,abs(transmit_shear_square_ratio));
title('square shear transmit ratio vs alpha');
xlabel('alpha(deg)');
ylabel('Ratio');
legend('real part', 'imaginary part', 'abs');






%% functions
function r = radius(source_x, source_y, x_arr,y_arr)
r = zeros(1,length(x_arr));
r = sqrt((x_arr - source_x).^2 + (y_arr - source_y).^2);
end