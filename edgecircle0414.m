clc;
clear;
close all;
%��ȡԭͼ
figure(1)
% poly������� Circle��Բ��
img = xlsread("testCircle2.xlsx");
image(img);
colorbar
colormap(gray(128));
title("ԭͼ");

%��ȡ���źŵĲ���
figure(2)
[row0,column0] = size(img);
imagesqu = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
image(imagesqu);
colorbar
colormap(gray(128));
title("ԭͼ��info�Ĳ���");

%��ȡcyst������ʾ
figure(3)
[row,column] = size(img);
% cystimg = img(0.15*row:0.35*row,0.4*column:0.6*column);
cystimg = img(0.3*row:0.55*row,0.35*column:0.65*column);
image(cystimg);
colorbar
colormap(gray(128));
title("��ȡcyst����");

%median filt
Img_filt = medfilt2(cystimg,[10,10]);
BWcyst = imbinarize(Img_filt);

%ת��Ϊ��ֵͼ��
% BWcyst = imbinarize(cystimg);
figure(4)
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);
title("��ֵ��cyst");
%�����
white_num = sum(Bcyst(:)==1);
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
figure(5)
edge1=edge1>=1;
imagesc(edge1);
colormap(gray);
title("δ��ֵ�߽�");

%%
%��ȡ���źŵĲ���
img_info = img(0.035*row0:0.695*row0,0.18*column0:0.83*column0);
%��ֵ�󳤿����
figure(6)
imagesqu = imresize(img_info,[1000,1000]);
% imshow(imagesqu);
image(imagesqu)
colorbar
colormap(gray(128));
title("��ֵ��������ͼ��");

%��ȡcyst������ʾ
figure(7)
[row,column] = size(imagesqu);
% cystimg = imagesqu(0.2*row:0.45*row,0.35*column:0.65*column);
cystimg = imagesqu(0.45*row:0.75*row,0.3*column:0.7*column);
image(cystimg);
colorbar
colormap(gray(128));
title("��ֵ���ȡcyst");

%median filt
Img_filt = medfilt2(cystimg,[10,10]);
BWcyst = imbinarize(Img_filt);

%ת��Ϊ��ֵͼ��
% BWcyst = imbinarize(cystimg);
figure(8)
Bcyst = ~BWcyst;
imagesc(Bcyst);
colormap(gray);
title("��ֵ���ֵ��cyst");

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
figure(9)
edge1=edge1>=1;
imagesc(edge1);
colormap(gray);
title("��ֵ��cyst�߽�");

% % ��һ��Բ��
% realcircle = zeros(151);
% r = 60;
% cx = 83;
% cy = 80;
% for i = 1:1000
%     for j = 1:1000
%         if((cx-i)^2 + (cy-j)^2) <= r^2+10 && ((cx-i)^2 + (cy-j)^2) >= r^2-10
%             realcircle(i,j)=1;
%         end
%     end
% end
% figure
% imagesc(realcircle);
% colormap(gray);

%% 
% 
% %��ȡcyst������ʾ
% figure(2)
% [row,column] = size(img);
% cystimg = img(0.16*row:0.31*row,0.4*column:0.6*column);
% image(cystimg);
% colorbar
% colormap(gray(128));
% %ת��Ϊ��ֵͼ��
% BWcyst = imbinarize(cystimg);
% figure
% Bcyst = ~BWcyst;
% imagesc(Bcyst);
% colormap(gray);
% 
% 
% %���ϵı߽�׷�ٳ���
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
% 
% figure;
% edge1=edge1>=1;
% imagesc(edge1);
% colormap(gray);
% 
% %ͨ����ԭͼ��ȥ�丯ʴͼ������
% % se = strel('square',3); 
% % edge1=Bcyst-imerode(Bcyst,se);    
% % figure;
% % imagesc(edge1);
% % colormap(gray);
% 
