clc;
clear;
close all;
% poly������� Circle��Բ��
img = xlsread("testpoly.xlsx");
% img = xlsread("testCircle.xlsx");     %Բ��

%��ȡ���źŵĲ���
figure(1)
[row0,column0] = size(img);
imagesqu = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
image(imagesqu);
colorbar
colormap(gray(128));
title("ԭͼ��info�Ĳ���");

%��ȡcyst������ʾ
figure(2)
[row,column] = size(imagesqu);
% img_cyst = imagesqu(0.15*row:0.45*row,0.35*column:0.75*column);   %poly2λ��
img_cyst = imagesqu(0.45*row:0.75*row,0.3*column:0.8*column); %poly1λ��
image(img_cyst);
colorbar
colormap(gray(128));
title('��ȡcyst����');

%median filt
Img_filt = medfilt2(img_cyst,[10,10]);
%ת��Ϊ��ֵͼ��
BWcyst = imbinarize(Img_filt);
figure
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);
title("��ֵ��cyst");

%�߽�׷��
[m,n]=size(BWcyst);
edge1=zeros(m,n);        %�߽���ͼ��
ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %�����Ͻ����أ���ʱ������
for i=2:m-1
    for j=2:n-1
        if Bcyst(i,j)==1 && edge1(i,j)==0      %��ǰ��û��ǵİ�ɫ����
            if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %���ڲ��İ����ز����
                ii=i;         %���ؿ��ڲ���Ѱʹ�õ�����
                jj=j;
                edge1(i,j)=2;    %�����ؿ��һ����ǵı߽磬��һ���߽�����Ϊ2
                
                while edge1(ii,jj)~=2    %�Ƿ��������ؿ���ѰһȦ�ˡ�
                    for k=1:8           %��ʱ�����������
                        tmpi=ii+ed(k,1);        %��������ʱ����
                        tmpj=jj+ed(k,2);
                        if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %�������±߽磬����û������һȦ
                            ii=tmpi;        %�����ڲ���Ѱ���꣬��������
                            jj=tmpj;
                            edge1(ii,jj)=1;  %�߽���ͼ������ر�ǣ���ͨ�߽�Ϊ1
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

%�㲻��ֵ���
white_num = sum(Bcyst(:)==1);
%���������
z_start = 30;
L = linspace(0,2.*pi,6);
xv = 8*cos(L)';
zv = 8*sin(L)'+ z_start*2;
figure
plot(xv,zv);
A = polyarea(xv,zv);

%% 
img_info = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
%��ֵ�󳤿����
figure
imagesqu = imresize(img_info,[1000,1000]);
% imshow(imagesqu);
image(imagesqu)
colorbar
colormap(gray(128));
title("��ֵ��������ͼ��");

%��ȡcyst������ʾ
figure
[row,column] = size(imagesqu);
% img_cyst = imagesqu(0.15*row:0.45*row,0.35*column:0.75*column);   %poly2λ��
img_cyst = imagesqu(0.45*row:0.75*row,0.3*column:0.8*column); %poly1λ��
%img = imagesqu(0.15*row:0.35*row,0.3*column:0.7*column);   %Բ��
% img_cyst = imagesqu(0.25*row:0.6*row,0.3*column:0.7*column);
image(img_cyst);
colorbar
colormap(gray(128));
title('��ֵ���ȡcyst');

%median filt
Img_filt = medfilt2(img_cyst,[10,10]);

%ת��Ϊ��ֵͼ��
% BWcyst = imbinarize(img_cyst);
BWcyst = imbinarize(Img_filt);
figure
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);

%�����
white_num2 = sum(Bcyst(:)==1);

%�߽�׷��
[m,n]=size(BWcyst);
edge1=zeros(m,n);        %�߽���ͼ��
ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %�����Ͻ����أ���ʱ������
for i=2:m-1
    for j=2:n-1
        if Bcyst(i,j)==1 && edge1(i,j)==0      %��ǰ��û��ǵİ�ɫ����
            if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %���ڲ��İ����ز����
                ii=i;         %���ؿ��ڲ���Ѱʹ�õ�����
                jj=j;
                edge1(i,j)=2;    %�����ؿ��һ����ǵı߽磬��һ���߽�����Ϊ2
                
                while edge1(ii,jj)~=2    %�Ƿ��������ؿ���ѰһȦ�ˡ�
                    for k=1:8           %��ʱ�����������
                        tmpi=ii+ed(k,1);        %��������ʱ����
                        tmpj=jj+ed(k,2);
                        if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %�������±߽磬����û������һȦ
                            ii=tmpi;        %�����ڲ���Ѱ���꣬��������
                            jj=tmpj;
                            edge1(ii,jj)=1;  %�߽���ͼ������ر�ǣ���ͨ�߽�Ϊ1
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
