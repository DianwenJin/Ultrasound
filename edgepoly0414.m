clc;
clear;
close all;
%读取原图
figure(1)
% poly是五边形 Circle是圆形
img = xlsread("testpoly.xlsx");
image(img);
colorbar
colormap(gray(128));
%插值后长宽相等
figure
imagesqu = imresize(img,[1000,1000]);
% imshow(imagesqu);
image(imagesqu)
colorbar
colormap(gray(128));

%截取cyst部分显示
figure
[row,column] = size(imagesqu);
cystimg = imagesqu(0.25*row:0.6*row,0.3*column:0.7*column);
image(cystimg);
colorbar
colormap(gray(128));
title('原图像');

%画一个五边形，失败
% [m,n]=size(cystimg);
% %withedge = cystimg;
% withedge = zeros(m,n);
% L = linspace(0,2.*pi,6);
% xv = 100*cos(L)'+ 200;
% zv = 100*sin(L)'+ 150;
% % figure
% % plot(xv,zv);
% in = zeros(m*n,1);
% on = zeros(m*n,1);
% for i = 1:m
%     for j =1:n
%     [in(i*j),on(i*j)] = inpolygon(i,j,xv,zv);
%     if on(i*j) == 1
%         withedge(i,j)=128;
%     end
%     end
% end
% figure
% image(withedge);
% colorbar
% colormap(gray(128));

%% 各种灰度算子 
BW1 = edge(cystimg,'canny');
% 其他5种不同提取子，但是结果差不多,都不理想
BW2 = edge(cystimg,'sobel');
BW3 = edge(cystimg,'prewitt');
BW4 = edge(cystimg,'roberts');
BW5 = edge(cystimg,'log');
BW6 = edge(cystimg,'zerocross');
figure
subplot(2,2,1)  
imagesc(BW1);
colormap(gray);
subplot(2,2,2)  
imagesc(BW2);
colormap(gray);
subplot(2,2,3)  
imagesc(BW3);
colormap(gray);
subplot(2,2,4)  
imagesc(BW4);
colormap(gray);


%% 没用的程序 并不能检测
% [m,n]=size(cystimg); 
% LEN= zeros(m,n); %建立矩阵保持对应像素的梯度幅度
% THETA = zeros(m,n); %建立矩阵保持对应像素的角
% JC = zeros(m,n); %建立边缘检测图矩阵
% hv = fspecial('prewitt'); 
% hh = hv.'; 
% gv = abs(imfilter(cystimg,hv,'replicate')); %水平方向梯度
% gh = abs(imfilter(cystimg,hh,'replicate'));  %竖直方向梯度
% LEN = sqrt(gv.^2 + gh.^2); %梯度幅值
% THETA = atan(gh./gv); %保存角度量
% for i=2:m-1  
%     for j=2:n-1 
%         if LEN(i,j) < 127 %建立删选条件
%             JC(i,j) = 0; 
%         else  
%             JC(i,j) = 1; 
%         end
%     end 
% end  
% figure
% imshow(JC); 
% title('prewitt边缘检测图');

%%
% % 转化为二值图像
% % BWcyst = im2bw(cystimg,0.2);
% BWcyst = imbinarize(cystimg);
% figure
% Bcyst = ~BWcyst;
% imagesc(Bcyst);
% colormap(gray);
% 
% %边界追踪
% [m,n]=size(BWcyst);
% edge1=zeros(m,n);        %边界标记图像
% ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %从左上角像素，逆时针搜索
% for i=2:m-1
%     for j=2:n-1
%         if Bcyst(i,j)==1 && edge1(i,j)==0      %当前是没标记的白色像素
%             if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %块内部的白像素不标记
%                 ii=i;         %像素块内部搜寻使用的坐标
%                 jj=j;
%                 edge1(i,j)=2;    %本像素块第一个标记的边界，第一个边界像素为2
%                 
%                 while edge1(ii,jj)~=2    %是否沿着像素块搜寻一圈了。
%                     for k=1:8           %逆时针八邻域搜索
%                         tmpi=ii+ed(k,1);        %八邻域临时坐标
%                         tmpj=jj+ed(k,2);
%                         if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %搜索到新边界，并且没有搜索一圈
%                             ii=tmpi;        %更新内部搜寻坐标，继续搜索
%                             jj=tmpj;
%                             edge1(ii,jj)=1;  %边界标记图像该像素标记，普通边界为1
%                             break;
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% figure;
% edge1=edge1>=1;
% imagesc(edge1);
% colormap(gray);

