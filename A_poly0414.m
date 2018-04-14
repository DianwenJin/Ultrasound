clc;
clear;
close all;
% poly是五边形 Circle是圆形
img = xlsread("testpoly.xlsx");
% img = xlsread("testCircle.xlsx");     %圆形

%截取有信号的部分
figure(1)
[row0,column0] = size(img);
imagesqu = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
image(imagesqu);
colorbar
colormap(gray(128));
title("原图有info的部分");

%截取cyst部分显示
figure(2)
[row,column] = size(imagesqu);
% img_cyst = imagesqu(0.15*row:0.45*row,0.35*column:0.75*column);   %poly2位置
img_cyst = imagesqu(0.45*row:0.75*row,0.3*column:0.8*column); %poly1位置
image(img_cyst);
colorbar
colormap(gray(128));
title('截取cyst部分');

%median filt
Img_filt = medfilt2(img_cyst,[10,10]);
%转化为二值图像
BWcyst = imbinarize(Img_filt);
figure
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);
title("二值化cyst");

%边界追踪
[m,n]=size(BWcyst);
edge1=zeros(m,n);        %边界标记图像
ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %从左上角像素，逆时针搜索
for i=2:m-1
    for j=2:n-1
        if Bcyst(i,j)==1 && edge1(i,j)==0      %当前是没标记的白色像素
            if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %块内部的白像素不标记
                ii=i;         %像素块内部搜寻使用的坐标
                jj=j;
                edge1(i,j)=2;    %本像素块第一个标记的边界，第一个边界像素为2
                
                while edge1(ii,jj)~=2    %是否沿着像素块搜寻一圈了。
                    for k=1:8           %逆时针八邻域搜索
                        tmpi=ii+ed(k,1);        %八邻域临时坐标
                        tmpj=jj+ed(k,2);
                        if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %搜索到新边界，并且没有搜索一圈
                            ii=tmpi;        %更新内部搜寻坐标，继续搜索
                            jj=tmpj;
                            edge1(ii,jj)=1;  %边界标记图像该像素标记，普通边界为1
                            break;
                        end
                    end
                end
                
            end
        end
    end
end
figure
edge1=edge1>=1;
imagesc(edge1);
colormap(gray);

%算不插值面积
white_num = sum(Bcyst(:)==1);
%算理论面积
z_start = 30;
L = linspace(0,2.*pi,6);
xv = 8*cos(L)';
zv = 8*sin(L)'+ z_start*2;
figure
plot(xv,zv);
A = polyarea(xv,zv);

%% 
img_info = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
%插值后长宽相等
figure
imagesqu = imresize(img_info,[1000,1000]);
% imshow(imagesqu);
image(imagesqu)
colorbar
colormap(gray(128));
title("插值后正方形图像");

%截取cyst部分显示
figure
[row,column] = size(imagesqu);
% img_cyst = imagesqu(0.15*row:0.45*row,0.35*column:0.75*column);   %poly2位置
img_cyst = imagesqu(0.45*row:0.75*row,0.3*column:0.8*column); %poly1位置
%img = imagesqu(0.15*row:0.35*row,0.3*column:0.7*column);   %圆形
% img_cyst = imagesqu(0.25*row:0.6*row,0.3*column:0.7*column);
image(img_cyst);
colorbar
colormap(gray(128));
title('插值后截取cyst');

%median filt
Img_filt = medfilt2(img_cyst,[10,10]);

%转化为二值图像
% BWcyst = imbinarize(img_cyst);
BWcyst = imbinarize(Img_filt);
figure
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);

%算面积
white_num2 = sum(Bcyst(:)==1);

%边界追踪
[m,n]=size(BWcyst);
edge1=zeros(m,n);        %边界标记图像
ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %从左上角像素，逆时针搜索
for i=2:m-1
    for j=2:n-1
        if Bcyst(i,j)==1 && edge1(i,j)==0      %当前是没标记的白色像素
            if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %块内部的白像素不标记
                ii=i;         %像素块内部搜寻使用的坐标
                jj=j;
                edge1(i,j)=2;    %本像素块第一个标记的边界，第一个边界像素为2
                
                while edge1(ii,jj)~=2    %是否沿着像素块搜寻一圈了。
                    for k=1:8           %逆时针八邻域搜索
                        tmpi=ii+ed(k,1);        %八邻域临时坐标
                        tmpj=jj+ed(k,2);
                        if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %搜索到新边界，并且没有搜索一圈
                            ii=tmpi;        %更新内部搜寻坐标，继续搜索
                            jj=tmpj;
                            edge1(ii,jj)=1;  %边界标记图像该像素标记，普通边界为1
                            break;
                        end
                    end
                end
                
            end
        end
    end
end
figure
edge1=edge1>=1;
imagesc(edge1);
colormap(gray);
