clc;
clear;
close all;
%��ȡԭͼ
figure(1)
% poly������� Circle��Բ��
img = xlsread("testpoly.xlsx");
image(img);
colorbar
colormap(gray(128));
%��ֵ�󳤿����
figure
imagesqu = imresize(img,[1000,1000]);
% imshow(imagesqu);
image(imagesqu)
colorbar
colormap(gray(128));

%��ȡcyst������ʾ
figure
[row,column] = size(imagesqu);
cystimg = imagesqu(0.25*row:0.6*row,0.3*column:0.7*column);
image(cystimg);
colorbar
colormap(gray(128));
title('ԭͼ��');

%��һ������Σ�ʧ��
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

%% ���ֻҶ����� 
BW1 = edge(cystimg,'canny');
% ����5�ֲ�ͬ��ȡ�ӣ����ǽ�����,��������
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


%% û�õĳ��� �����ܼ��
% [m,n]=size(cystimg); 
% LEN= zeros(m,n); %�������󱣳ֶ�Ӧ���ص��ݶȷ���
% THETA = zeros(m,n); %�������󱣳ֶ�Ӧ���صĽ�
% JC = zeros(m,n); %������Ե���ͼ����
% hv = fspecial('prewitt'); 
% hh = hv.'; 
% gv = abs(imfilter(cystimg,hv,'replicate')); %ˮƽ�����ݶ�
% gh = abs(imfilter(cystimg,hh,'replicate'));  %��ֱ�����ݶ�
% LEN = sqrt(gv.^2 + gh.^2); %�ݶȷ�ֵ
% THETA = atan(gh./gv); %����Ƕ���
% for i=2:m-1  
%     for j=2:n-1 
%         if LEN(i,j) < 127 %����ɾѡ����
%             JC(i,j) = 0; 
%         else  
%             JC(i,j) = 1; 
%         end
%     end 
% end  
% figure
% imshow(JC); 
% title('prewitt��Ե���ͼ');

%%
% % ת��Ϊ��ֵͼ��
% % BWcyst = im2bw(cystimg,0.2);
% BWcyst = imbinarize(cystimg);
% figure
% Bcyst = ~BWcyst;
% imagesc(Bcyst);
% colormap(gray);
% 
% %�߽�׷��
% [m,n]=size(BWcyst);
% edge1=zeros(m,n);        %�߽���ͼ��
% ed=[-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0]; %�����Ͻ����أ���ʱ������
% for i=2:m-1
%     for j=2:n-1
%         if Bcyst(i,j)==1 && edge1(i,j)==0      %��ǰ��û��ǵİ�ɫ����
%             if sum(sum(Bcyst(i-1:i+1,j-1:j+1)))~=9    %���ڲ��İ����ز����
%                 ii=i;         %���ؿ��ڲ���Ѱʹ�õ�����
%                 jj=j;
%                 edge1(i,j)=2;    %�����ؿ��һ����ǵı߽磬��һ���߽�����Ϊ2
%                 
%                 while edge1(ii,jj)~=2    %�Ƿ��������ؿ���ѰһȦ�ˡ�
%                     for k=1:8           %��ʱ�����������
%                         tmpi=ii+ed(k,1);        %��������ʱ����
%                         tmpj=jj+ed(k,2);
%                         if Bcyst(tmpi,tmpj)==1 && edge1(tmpi,tmpj)~=2  %�������±߽磬����û������һȦ
%                             ii=tmpi;        %�����ڲ���Ѱ���꣬��������
%                             jj=tmpj;
%                             edge1(ii,jj)=1;  %�߽���ͼ������ر�ǣ���ͨ�߽�Ϊ1
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

