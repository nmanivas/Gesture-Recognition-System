counter=1;
for f=128:144
str=['C:\Users\valli\Dropbox\Photos\project common\train\a\snapx',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                    %if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(1,counter)=1;
counter=counter+1;

end
for f=7:24
str=['C:\Users\valli\Dropbox\Photos\project common\train\d\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                    %if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(2,counter)=1;
counter=counter+1;

end

for f=529:544
str=['C:\Users\valli\Dropbox\Photos\project common\train\g\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                   % if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(3,counter)=1;
counter=counter+1;

end
for f=537:557
str=['C:\Users\valli\Dropbox\Photos\project common\train\m\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                   % if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(4,counter)=1;
counter=counter+1;

end

for f=536:556
str=['C:\Users\valli\Dropbox\Photos\project common\train\t\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(6,counter)=1;
counter=counter+1;

end
for f=523:539
str=['C:\Users\valli\Dropbox\Photos\project common\train\z\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(7,counter)=1;
counter=counter+1;

end

for f=22:43
str=['C:\Users\valli\Dropbox\Photos\project common\train\b\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                    %if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(8,counter)=1;
counter=counter+1;

end

for f=543:561
str=['C:\Users\valli\Dropbox\Photos\project common\train\x\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                   % if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(9,counter)=1;
counter=counter+1;

end

for f=161:174
str=['C:\Users\valli\Dropbox\Photos\project common\train\r\snapx',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
               %     if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                %    end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(10,counter)=1;
counter=counter+1;

end


for f=545:561
str=['C:\Users\valli\Dropbox\Photos\project common\train\p\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(11,counter)=1;
counter=counter+1;

end



for f=545:554
str=['C:\Users\valli\Dropbox\Photos\project common\train\o\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(13,counter)=1;
counter=counter+1;

end


for f=539:550
str=['C:\Users\valli\Dropbox\Photos\project common\train\h\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x(index,counter)=H(index);

end
for i= 1:14
    
        t(i,counter)=0;
    
end
t(14,counter)=1;
counter=counter+1;

end


size(x)
size(t)



%rand('seed', 491218382)

%net = Neural Network object:
 %   architecture:
  %      numInputs: 0
   %     numLayers: 0
  net = network;
        net=init(net);
        net.numInputs=1;
        net.numLayers=3;
        net.inputConnect(1,1) = 1;
%net.inputConnect(2,1) = 1;
%net.inputConnect(2,2) = 1;
         %net.inputConnect: [0 0 0; 1 0 0; 0 1 0]
        net.biasConnect(1)=1;
        net.biasConnect(2)=1;
        net.biasConnect(3)=1;
        net.layerConnect = [0 0 0; 1 0 0; 0 1 0];
        net.outputConnect = [0 0 1];
        net.inputs{1}.size=108;
         net.layers{1}.size=21;
             net.layers{2}.size=21;
  %               net.layers{3}.size=10;
        net.layers{1}.transferFcn='tansig';
        net.layers{1}.initFcn='initnw';
            net.layers{2}.transferFcn='tansig';
%                  net.inputs{2}.size=10;
        net.layers{2}.initFcn='initnw';
 %            net.inputs{3}.size=10;
        net.layers{3}.transferFcn='tansig';
        net.layers{3}.initFcn='initnw';
        
        net.initFcn='initlay';
        net.performFcn = 'mse';
        net.trainFcn = 'trainlm';
        
%net = patternnet(10);
%net = feedforwardnet(25);
%net = cascadeforwardnet(10);
net = train(net,x,t);


[net,tr] = train(net,x,t);
view(net)
counterr=1;
for f=135:145
str=['C:\Users\valli\Dropbox\Photos\project common\train\y\snapx',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(1,counterr)=1;
counterr=counterr+1;

end
for f=531:547
str=['C:\Users\valli\Dropbox\Photos\project common\train\k\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                   % if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(2,counterr)=1;
counterr=counterr+1;

end
for f=541:555
str=['C:\Users\valli\Dropbox\Photos\project common\train\f\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);


b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(3,counterr)=1;
counterr=counterr+1;

end

for f=29:40
str=['C:\Users\valli\Dropbox\Photos\project common\train\e\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(4,counterr)=1;
counterr=counterr+1;

end


for f=163:174
str=['C:\Users\valli\Dropbox\Photos\project common\train\n\snapx',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                   % if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(5,counterr)=1;
counterr=counterr+1;

end

for f=161:177
str=['C:\Users\valli\Dropbox\Photos\project common\train\s\snapx',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                  %  if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                   % end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(6,counterr)=1;
counterr=counterr+1;

end

for f=537:551
str=['C:\Users\valli\Dropbox\Photos\project common\train\q\',num2str(f),'.jpg'];%change path over here.
Ima=imread(str);
Ima=rgb2gray(Ima);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
len=length(b11);
areas=[b11.Area];
[maxArea largestBlobIndex]=max(areas);
x1=b11(largestBlobIndex).BoundingBox(1);
y1=b11(largestBlobIndex).BoundingBox(2);
x2=b11(largestBlobIndex).BoundingBox(3);
y2=b11(largestBlobIndex).BoundingBox(4);
Ic=imcrop(Ima,[x1 y1 x2 y2]);
imshow(Ic);
 
 
 Im=Ic;
  %figure(11),imshow(Im);
 Im=imfill(Im,8);
 Im=medfilt2(Im,[3 3]);
 figure(10),imshow(Im);
 %Im=imfilter(Im,'gaussian');
 Im=imresize(Im,[90,90]);
nwin_x=3;%set here the number of HOG windows per bound box
nwin_y=3;
B=12;%set here the number of hivstogram bins
[L,C]=size(Im); % L num of lines ; C num of columns
H=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
m=sqrt(L/2);
if C==1 % if num of columns==1
    Im=im_recover(Im,m,2*m);%verify the size of image, e.g. 25x50
    L=2*m;
    C=m;
end
Im=double(Im);
step_x=floor(C/(nwin_x));
step_y=floor(L/(nwin_y));
cont=0;
hx = [-1,0,1];
hy = -hx';
gr=0;
grad_xr = imfilter(double(Im),hx);
grad_yu = imfilter(double(Im),hy);
angles=atan2(grad_yu,grad_xr);
magnit=((grad_yu.^2)+(grad_xr.^2)).^.5;
figure(1),imshow(magnit,[]);
hold on;

quiver(grad_yu,grad_xr);
for i=1:90
    for j=1:90
if(magnit(i,j)<40)
    grad_xr(i,j)=0;
    grad_yu(i,j)=0;
end
end
end

%figure(3),quiver(grad_yu,grad_xr);

for n=0:(nwin_y-1)
    for m=0:(nwin_x-1)
        cont=cont+1;
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot(m-1*step_x,n-1*step_y,'ro');
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
        %rectangle('Position',[(m)*step_x+1,(n)*step_y+1,(m+1)*step_x,(n+1)*step_y],'Edgecolor','r');
        %s=[(m)*step_x+1 (n)*step_y+1 (m+1)*step_x (n+1)*step_y];
        %disp(s);
        %pause(1);
        v_angles=angles2(:);    
        v_magnit=magnit2(:);
        K=max(size(v_angles));
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        bin=0;
        H2=zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            bin=bin+1;
            for k=1:K
                if v_angles(k)<ang_lim
                    v_angles(k)=100;
                    %if(v_magnit(k)>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    
                    %end
                end
            end
        end
        gr=gr+1;
        H3(gr,:)=H2;
        disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\newhog424\ah',num2str(f),'.txt'];
%save(str1,'H','-ASCII');
f1=getframe(gca);
pic=frame2im(f1);
%str=['C:\Users\valli\Desktop\wind',num2str(1),'.jpg'];   
%imwrite(pic,str,'jpg');
hold off;
for index=1:108
x1b(index,counterr)=H(index);

end
for i= 1:7
    
        t1(i,counterr)=0;
    
end
t1(7,counterr)=1;
counterr=counterr+1;

end

 net1 = network;
        net1=init(net1);
        net1.numInputs=1;
        net1.numLayers=3;
        net1.inputConnect(1,1) = 1;
%net.inputConnect(2,1) = 1;
%net.inputConnect(2,2) = 1;
         %net.inputConnect: [0 0 0; 1 0 0; 0 1 0]
        net1.biasConnect(1)=1;
        net1.biasConnect(2)=1;
        net1.biasConnect(3)=1;
        net1.layerConnect = [0 0 0; 1 0 0; 0 1 0];
        net1.outputConnect = [0 0 1];
        net1.inputs{1}.size=108;
         net1.layers{1}.size=21;
             net1.layers{2}.size=21;
  %               net.layers{3}.size=10;
        net1.layers{1}.transferFcn='tansig';
        net1.layers{1}.initFcn='initnw';
            net1.layers{2}.transferFcn='tansig';
%                  net.inputs{2}.size=10;
        net1.layers{2}.initFcn='initnw';
 %            net.inputs{3}.size=10;
        net1.layers{3}.transferFcn='tansig';
        net1.layers{3}.initFcn='initnw';
        
        net1.initFcn='initlay';
        net1.performFcn = 'mse';
        net1.trainFcn = 'trainlm';
        
%net = patternnet(10);
%net = feedforwardnet(25);
%net = cascadeforwardnet(10);
net1 = train(net1,x1b,t1);


[net1,tr] = train(net1,x1b,t1);
view(net1)
save net;
save net1;