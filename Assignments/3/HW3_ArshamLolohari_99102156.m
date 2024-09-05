%% Q1 and Q2
clc;


% create the computational grid
Nx = 250;           % number of grid points in the x (row) direction
Ny = 250;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

medium.sound_speed = water_speed;
medium.density = water_density;


% load the initial pressure distribution from an image and scale the magnitude
p0_magnitude = 3;
source_pic = loadImage('EDITED_example_pr_2D_tr_circular_sensor_01.png');
p0 = p0_magnitude * source_pic;
% resize the image to match the size of the computational grid and assign to the source input structure
resized_source = resize(p0, [Nx, Ny]);
source.p0 = resized_source;


% create a binary sensor mask of an equivalent continuous circle
sensor_radius = 8 * 10^(-3); %[m]
sensor_angle = 4*pi/2;

gridCenter_x = round(kgrid.Nx/2) + 1;
gridCenter_y = round(kgrid.Ny/2) + 1;

sensor_radius_grid_points = round(sensor_radius / kgrid.dx);
binary_sensor_mask = makeCircle(kgrid.Nx, kgrid.Ny, gridCenter_x, gridCenter_y, sensor_radius_grid_points, sensor_angle);

% assign to sensor structure
% num_sensor_points = length(find(binary_sensor_mask));
num_sensor_points = 400; %Change it to 100 for question 1
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points);

sensor.mask = cart_sensor_mask;

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor); %dimension = num_sensor_points*Nt

dt = kgrid.dt;
Nt = kgrid.Nt;



sensors_x = sensor.mask(1,:); %dimension = 1*num_sensor_points
sensors_y = sensor.mask(2,:); %dimension = 1*num_sensor_points

shifted_sensor_data = zeros(num_sensor_points,Nt);
grid_vals = zeros(Nx,Ny); %Grid points final picture value

for i=1:Nx
    disp(i); %Show the progress
    for j=1:Ny
        point_x = kgrid.x_vec(i);
        point_y = kgrid.y_vec(j);
        
        sensors_dist = sqrt((sensors_x - point_x).^2 + (sensors_y - point_y).^2); %dimension = 1*num_sensor_points
        sensors_timeShift = sensors_dist / medium.sound_speed; %[sec]. Dimension = 1*num_sensor_points
        sensors_sampleShift = round(sensors_timeShift / dt); %[samples]. Dimension = 1*num_sensor_points
        
        for s = 1:num_sensor_points
           shifted_sensor_data(s,:) =  circshift(sensor_data(s,:), -sensors_sampleShift(s)); %dimension for each sensor = 1*num_sensor_points
           shifted_sensor_data(s,end-sensors_sampleShift(s):end) =  0; %dimension for each sensor = 1*num_sensor_points
            
        end
        
        shiftedSum_sensor_data = sum((shifted_sensor_data)); %summing over sensors. Dimension = 1*Nt
        grid_vals(i,j) = abs(shiftedSum_sensor_data(1));
    end
end

figure;
imagesc(kgrid.y_vec*1000, kgrid.x_vec*1000, source_pic);
title("Original image");
xlabel('y pos[mm]');
ylabel('x pos[mm]');
colorbar;
colormap jet;

figure;
imagesc(kgrid.y_vec*1000, kgrid.x_vec*1000, grid_vals);
title("DAS Reconstructed image using "+num_sensor_points+" sensors");
xlabel('y pos[mm]');
ylabel('x pos[mm]');
colorbar;
colormap jet;



%% Q3 & first part of Q5

% create the computational grid
Nx = 250;           % number of grid points in the x (row) direction
Ny = 250;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

medium.sound_speed = water_speed;
medium.density = water_density;


% load the initial pressure distribution from an image and scale the magnitude
p0_magnitude = 3;
source_pic = loadImage('EDITED_example_pr_2D_tr_circular_sensor_01.png');
p0 = p0_magnitude * source_pic;
% resize the image to match the size of the computational grid and assign to the source input structure
resized_source = resize(p0, [Nx, Ny]);
source.p0 = resized_source;


% Create sensor structure
sensor_radius = 8 * 10^(-3); %[m]
sensor_angle = 4*pi/2;

gridCenter_x = round(kgrid.Nx/2) + 1;
gridCenter_y = round(kgrid.Ny/2) + 1;

sensor_radius_grid_points = round(sensor_radius / kgrid.dx);

% assign to sensor structure
num_sensor_points = 400;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points);

sensor.mask = cart_sensor_mask;

% run the simulation
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor); %dimension = num_sensor_points*Nt

dt = kgrid.dt;
Nt = kgrid.Nt;

fc = 2.25e6; %Change it to 6.25e6 for part 1 in question 5
bandw = 0.6*fc;

fs = 1/dt;
filtered_sensor_data = bandFilt(sensor_data2.',fc,bandw,fs).'; %dimension = num_sensor_points*Nt


shifted_sensor_data = zeros(num_sensor_points,Nt);
grid_vals = zeros(Nx,Ny); %Grid points final picture value

sensors_x = sensor.mask(1,:); %dimension = 1*num_sensor_points
sensors_y = sensor.mask(2,:); %dimension = 1*num_sensor_points


for i=1:Nx
    disp(i); %Show the progress
    for j=1:Ny
        point_x = kgrid.x_vec(i);
        point_y = kgrid.y_vec(j);
        
        sensors_dist = sqrt((sensors_x - point_x).^2 + (sensors_y - point_y).^2); %dimension = 1*num_sensor_points
        sensors_timeShift = sensors_dist / medium.sound_speed; %[sec]. Dimension = 1*num_sensor_points
        sensors_sampleShift = round(sensors_timeShift / dt); %[samples]. Dimension = 1*num_sensor_points
        
        for s = 1:num_sensor_points
           shifted_sensor_data(s,:) =  circshift(filtered_sensor_data(s,:), -sensors_sampleShift(s)); %dimension for each sensor = 1*num_sensor_points
           shifted_sensor_data(s,end-sensors_sampleShift(s):end) =  0; %dimension for each sensor = 1*num_sensor_points
            
        end
        
        shiftedSum_sensor_data = sum(shifted_sensor_data); %summing over sensors. Dimension = 1*Nt
        grid_vals(i,j) = abs(shiftedSum_sensor_data(1));
    end
end

figure;
imagesc(kgrid.y_vec*1000, kgrid.x_vec*1000, grid_vals);
title("Recon. image using " + num_sensor_points+" sensors, filter on fc = "+fc*1e-6+" MHz and bandwidth = "+bandw/fc*100+"%");
xlabel('y pos[mm]');
ylabel('x pos[mm]');

colorbar;
colormap jet;







%% Q4
clc;

% create the computational grid
Nx = 270;           % number of grid points in the x (row) direction
Ny = 270;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);




%defining medium
water_speed = 1480;  % [m/s]
water_density = 1000; %[kg/m^3]

medium.sound_speed = water_speed;
medium.density = water_density;



% load the initial pressure distribution from an image and scale the magnitude
p0_magnitude = 3;
source_pic = loadImage('EDITED_example_pr_2D_tr_circular_sensor_01.png');
p0 = p0_magnitude * source_pic;
% resize the image to match the size of the computational grid and assign to the source input structure
resized_source = resize(p0, [Nx, Ny]);
source.p0 = resized_source;


sensor_radius = 9* 10^(-3); %[m]
sensor_angle = 4*pi/2;

gridCenter_x = round(kgrid.Nx/2) + 1;
gridCenter_y = round(kgrid.Ny/2) + 1;

sensor_radius_grid_points = round(sensor_radius / kgrid.dx);

% assign to sensor structure
num_sensor_points = 400;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points);


sensors_x = cart_sensor_mask(1,:); %dimension = 1*num_sensor_points
sensors_y = cart_sensor_mask(2,:); %dimension = 1*num_sensor_points

%NOW: extending each transducer to a length of 6mm
ext_cart_sensor_mask = cart_sensor_mask;
sensors_dx = (cart_sensor_mask(1,:) - (1*dx)); %dimension = 1*num_sensor_points
sensors_dy = (cart_sensor_mask(2,:) - (1*dy)); %dimension = 1*num_sensor_points

sensors_dx(sensors_dx<0) = -dx;
sensors_dx(sensors_dx>0) = dx;
sensors_dx(sensors_dx==0) = 0;

sensors_dy(sensors_dy<0) = -dy;
sensors_dy(sensors_dy>0) = dx;
sensors_dy(sensors_dy==0) = 0;



totalLen = 0; %total length of each transducer [m].

s = 1;
while totalLen < 0.006
    newSensors_x = sensors_x + s * sensors_dx; %dimension = 1*num_sensor_points
    newSensors_y = sensors_y + s * sensors_dy; %dimension = 1*num_sensor_points
    ext_cart_sensor_mask = [ext_cart_sensor_mask, [newSensors_x;newSensors_y]];
    totalLen = s * sqrt(sensors_dx(1)^2 + sensors_dy(1)^2);
    s = s+1;
end

sensor.mask = ext_cart_sensor_mask;

% run the simulation
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor); %dimension = num_sensor_points*Nt

dt = kgrid.dt;
Nt = kgrid.Nt;

fc = 2.25e6; %Change it to 6.25 for part 2 in question 5
bandw = 0.6*fc;

fs = 1/dt;
filtered_sensor_data = bandFilt(sensor_data2.',fc,bandw,fs).'; %dimension = num_sensor_points*Nt


shifted_sensor_data = zeros(num_sensor_points,Nt);
grid_vals = zeros(Nx,Ny); %Grid points final picture value


for i=1:Nx
    disp(i); %Show the progress
    for j=1:Ny
        point_x = kgrid.x_vec(i);
        point_y = kgrid.y_vec(j);
        
        sensors_dist = sqrt((sensors_x - point_x).^2 + (sensors_y - point_y).^2); %dimension = 1*num_sensor_points
        sensors_timeShift = sensors_dist / medium.sound_speed; %[sec]. Dimension = 1*num_sensor_points
        sensors_sampleShift = round(sensors_timeShift / dt); %[samples]. Dimension = 1*num_sensor_points
        
        for s = 1:num_sensor_points
           shifted_sensor_data(s,:) =  circshift(filtered_sensor_data(s,:), -sensors_sampleShift(s)); %dimension for each sensor = 1*num_sensor_points
           shifted_sensor_data(s,end-sensors_sampleShift(s):end) =  0; %dimension for each sensor = 1*num_sensor_points
            
        end
        
        shiftedSum_sensor_data = sum(shifted_sensor_data); %summing over sensors. Dimension = 1*Nt
        grid_vals(i,j) = abs(shiftedSum_sensor_data(1));
    end
end

figure;
imagesc(kgrid.y_vec*1000, kgrid.x_vec*1000, grid_vals);
title("Recon. image with extended sensors, using " + num_sensor_points+" sensors, filter on fc = "+fc*1e-6+" MHz and bandwidth = "+bandw/fc*100+"%");
xlabel('y pos[mm]');
ylabel('x pos[mm]');

colorbar;
colormap jet;







%% functions

function filt_data = bandFilt(data,fc,bandw,fs) %band-pass filter on data
    %assumption: data dimension=samples*sensors
    
    butter_filt=designfilt('bandpassiir','FilterOrder',8,...
        'HalfPowerFrequency1',fc-bandw/2,'HalfPowerFrequency2',fc+bandw/2,'SampleRate',fs,'DesignMethod','butter');
%     fvtool(butter_filt);
    filt_data=filter(butter_filt,data); %dimension=samples*channels
end
