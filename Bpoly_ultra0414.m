%%
clc;
clear;
close all;
path(path,'C:\Users\Classic\Documents\MATLAB\Field_II_PC7');
field_init;

%% 参数设置
% Generate the transducer apertures for send and receive
f0=3e6; % Transducer center frequency [Hz]
fs=120e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wave length [m]
width=lambda; % Width of element
element_height=5/1000; % Height of element [m]
kerf=width/20; % Kerf width of interval between elements [m]
focus=[0 0 50]/1000; % Fixed focal point [m]
N_elements=192; % Number of elements in the transducer
N_active=64; % Active elements in the transducer

%% 定义发射、接受换能器
% Set the sampling frequency of the system
set_sampling(fs);
% Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response 设置换能器脉冲响应，加窗
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
% Set excitation of the emit aperture 为发射换能器设置激励
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);
% Generate aperture for reception
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

%% 添加一个设置attenuation
set_field('att',1.5*100);
set_field ('Freq_att',0.54*100/1e6);
set_field ('att_f0',3e6);
set_field ('use_att',1);

%% 加载phantom进行处理
% Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
% Do linear array imaging
no_lines=N_elements-N_active+1; % Number of A-lines in image
dx=width; % Increment for image
z_focus=50/1000;
% Pre-allocate some storage
image_data=zeros(1,no_lines);
for i=1:no_lines
    % Find position for imaging
    x=(i-1-no_lines/2)*dx;
    % Set the focus for this direction
    xdc_center_focus (emit_aperture, [x 0 0]);
    xdc_focus (emit_aperture, 0, [x 0 z_focus]);
    xdc_center_focus (receive_aperture, [x 0 0]);
    xdc_focus (receive_aperture, 0, [x 0 z_focus]);
    % Set the active elements using the apodization
    apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
    xdc_apodization (emit_aperture, 0, apo);
    xdc_apodization (receive_aperture, 0, apo);
    % Calculate the received response
    [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    % Store the result
    image_data(1:max(size(v)),i)=v;
    times(i) = t1;
end
% Free space for apertures 释放资源
xdc_free (emit_aperture)
xdc_free (receive_aperture)

%% 对齐处理，对数压缩，灰度范围校正后显示图像
% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines
rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
env(1:size(rf_env,1),i)=rf_env;
end
% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;

% env_gray即为要展示的矩阵，x*1000, depth*1000作为x y定义图像在xy轴上的位置
image(x*1000, depth*1000, env_gray)
xlabel('lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('image of cyst phantom (60 db dynamic range)')

figure
image(env_gray);
colormap(gray(128));
xlswrite('C:\Users\Classic\Documents\MATLAB\testpoly.xlsx',[env_gray]);

%% 建立phantom
function [positions, amp] = cyst_phantom(N)
x_size = 40/1000; % Width of phantom [m]
y_size = 10/1000; % Transverse width of phantom [m]
z_size = 50/1000; % Height of phantom [m]
z_start = 30/1000; % Start of phantom surface [m];
% Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;

% %不知道scatter啥意思 画多一点试试看，发现不行
% x = (rand (N,1)-0.5)* 0.4;
% y = (rand (N,1)-0.5)*0.4;
% z = rand (N,1)*0.4;

% amp with a Gaussian distribution
amp=randn(N,1);
% Make the cyst and set the amplitudes to 0.3 inside
%定义五边形
L = linspace(0,2.*pi,6);
xv = 0.008*cos(L)';
zv = 0.008*sin(L)'+ z_start*2;
inornot = inpolygon(x,z,xv,zv);
inside = ( inornot == 1); %判断，内部则返回1
amp = amp.*(1-inside) + 0.15*amp.*inside;

% Place the point scatterers in the phantom
dz = z_size/10;
for i = N-9:N
    x(i) = -15/1000;
    y(i) = 0;
    z(i) = z_start + (i-N+9)*dz;
    amp(i) = 100;
end
%reture the variables
positions = [x y z];
end