%slow swept laser simulation
clear all;
clc;


%setting some prefixed
mm = 1e-3;
um = 1e-6;
nm = 1e-9;
kHz=1e+3;
MHz=1e+6;

                
load('airPuff.mat'); %read exemplary air puff deformation plot
las_freq=5*kHz; %laser repetition frequency
sam_freq=10*MHz; %sampling frequency MHz means MS/s
las_duty_cycle=0.5; %laser duty cycle
n_samples_per_sweep=las_duty_cycle*(sam_freq/las_freq); %samples per sweep
n_ascan=1000; %number of Ascans - here M-mode is considered

%laser spectral characteristics
lambda_start=1305*nm;
lambda_end=1315*nm;

%recalculation to k-space and creating sweep k-vector
k_start=(2*pi)/lambda_start;
k_end=(2*pi)/lambda_end;
k_vector=linspace(k_end,k_start,n_samples_per_sweep);

k2_start=(2*pi)/lambda_end;
k2_end=(2*pi)/lambda_start;
k2_vector=linspace(k_end,k_start,n_samples_per_sweep);
%creating gaussian envelope to simulate Gaussian shape source
x = 0:n_samples_per_sweep;
%gauss_envelope = gaussmf(x,[n_samples_per_sweep/8 n_samples_per_sweep/2]);

gauss_envelope = gauss_distribution(x,n_samples_per_sweep/2, n_samples_per_sweep/2);
gauss_envelope = gauss_envelope(:,1:end-1);


%some variables usefull to scale the final plot (time and depth)
x_time_max=n_ascan/las_freq; %max value on time scale
z_max=(pi*n_samples_per_sweep)/(2*abs(k_start-k_end)); %max depth range
y_scale=2*z_max/n_samples_per_sweep; %scaling in y axis

%defining the surface position (should have length corresponding to
%n_ascans
Delta_z_ini=10*mm;
Delta_z_max=10*mm;

%Delta=linspace(Delta_z_ini,Delta_z_max,n_ascan);
%Delta=sin(3*mm*k_vector)*2*mm+10*mm; %sinosuidal modulation of the surface

%surface deformation as the position input
% x = (0:1:size(air_puff_deform,1))';
% xi = (0:1:n_ascan)';
% air_puff_deform_inter = interp1q(x,air_puff_deform,xi);
%puff_interp = interp(Delta,n_samples_per_sweep);
x_time_air_puff=0.040; %defining time of the deformation plot to 40ms
displacement_amp=2.5*mm;


x = linspace(0, 1, size(air_puff_deform,1));
n_askan_puff=round(n_ascan*(x_time_air_puff/x_time_max));
x1 = linspace(0, 1, n_askan_puff);
air_puff_deform_inter = interp1(x, air_puff_deform, x1, 'linear','extrap');
x2 = linspace(0, 1, n_ascan-n_askan_puff);
xx = linspace(0, 1, size(air_puff_deform(end-2:end),1));
air_puff_no_deform = interp1(xx, air_puff_deform(end-2:end), x2, 'linear','extrap');


air_puff_deform_inter=[air_puff_deform_inter,air_puff_no_deform];
air_puff_deform_inter=displacement_amp*air_puff_deform_inter/max(air_puff_deform_inter);
Delta=air_puff_deform_inter+3.5*mm;
Deltatemp = interp(Delta,500);

xd = linspace(0, 1, size(Deltatemp,2));
xdd = linspace(0,1,2*size(Deltatemp,2));
Delta3=interp1(xd,Deltatemp,xdd);
%x3 = linspace(0, 1, size(Delta,2)*n_samples_per_sweep);
Delta2 = reshape(Deltatemp,[1000 500]);
%Delta2 = interp1(x3,Delta,x3);
Delta3 = reshape(Delta3,[1000 1000]);

surface_reflectivity=100.0;

n_fft=1024; %fft size
OCT_fringe=zeros(n_ascan,n_samples_per_sweep); %initialize OCT_fringe array
OCT_int=zeros(n_ascan,n_fft); %initialize OCT_fringe array
OCT_fringe2=zeros(n_ascan,n_samples_per_sweep); %initialize OCT_fringe array
OCT_int2=zeros(n_ascan,n_fft); %initialize OCT_fringe array

for ii=1:n_ascan
    noise_array = poissrnd(400,[1 n_samples_per_sweep])/10;
    OCT_fringe2(ii,:)=surface_reflectivity*((cos(2*Delta3(:,ii)'.*k2_vector)+1).*gauss_envelope).*noise_array;
    OCT_fringe(ii,:)=surface_reflectivity*((cos(2*Delta3(:,ii)'.*k_vector)+1).*gauss_envelope).*noise_array;
    OCT_int(ii,:)=20*log10(abs(fft(OCT_fringe(ii,:),n_fft)));
    OCT_int2(ii,:)=20*log10(abs(fft(OCT_fringe2(ii,:),n_fft)));
end
subplot(3,1,1);
imagesc(OCT_int(:,1:n_fft/2).');  
colormap(gray);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/las_freq)
xlabel('Time [s]')

yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', (yt*y_scale)*1000)
ylabel('Depth [mm]')


for iii=1:500
    
    [upperenv, lowerenv] = envelope(OCT_fringe(iii,:));
    OCT_fringe(iii,:)=OCT_fringe(iii,:)-0.5*upperenv;
    %octfft = fft2(OCT_fringe);
    hhh = hilbert(OCT_fringe(iii,:));
    sigphase(iii,:) = (unwrap(angle(hhh)))';
    
end

for iii=1:500
    
    [upperenv, lowerenv] = envelope(OCT_fringe2(iii,:));
    OCT_fringe2(iii,:)=OCT_fringe2(iii,:)-0.5*upperenv;
    %octfft = fft2(OCT_fringe);
    hhh = hilbert(OCT_fringe2(iii,:));
    sigphase2(iii,:) = (unwrap(angle(hhh)))';
    
end
subplot(3,1,2);
plot(OCT_fringe(1,:));

subplot(3,1,3);
imagesc(sigphase.');
colormap(jet);

figure(2);
subplot(2,1,1)
for iii=1:500
    
    xs = 1:size(sigphase(iii,:),2);
    
    pp = polyfit(xs,sigphase(iii,:),1);
    ys(iii,:) = polyval(pp,xs);
    slope(iii)=pp(1);
end
imagesc(OCT_int(:,1:n_fft/2).');
hold on;
plot(200*slope,'LineWidth',4);
subplot(2,1,2)
for iii=1:500
    
    xs = 1:size(sigphase2(iii,:),2);
    
    pp = polyfit(xs,sigphase2(iii,:),1);
    ys(iii,:) = polyval(pp,xs);
    slope2(iii)=pp(1);
end
imagesc(OCT_int2(:,1:n_fft/2).');
hold on;
plot(200*slope2,'LineWidth',4);

figure (3)
subplot(2,1,1);
imagesc(OCT_int(:,1:n_fft/2).');  
colormap(gray);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/las_freq)
xlabel('Time [s]')

yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', (yt*y_scale)*1000)
ylabel('Depth [mm]')

subplot(2,1,2);
imagesc(OCT_int2(:,1:n_fft/2).');  
colormap(gray);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/las_freq)
xlabel('Time [s]')

yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', (yt*y_scale)*1000)
ylabel('Depth [mm]')

