%for s1=1:52

%{
str=['C:\Users\valli\Desktop\train\a\',num2str(s1),'.jpg'];
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
%imshow(Ic);
%imshow(Ic);
Im=Ic; 
%figure(11),imshow(Im);
Im=imfill(Im,8);
Im=medfilt2(Im,[3 3]);
%figure(10),imshow(Im);
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
for n=0:nwin_y-1
    for m=0:nwin_x-1
        cont=cont+1;
       
        angles2=angles((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y); 
        magnit2=magnit((m)*step_x+1:(m+1)*step_x,(n)*step_y+1:(n+1)*step_y);
        %plot((m+1)*step_x,(n+1)*step_y,'c*');
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
                    %if(v_magnit>40)
                    H2(bin)=H2(bin)+v_magnit(k);
                    %end
                    v_angles(k)=100;
                end
            end
        end
        gr=gr+1;
       
        
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        H((cont-1)*B+1:cont*B,1)=H2;
        H3(gr,:)=H2;
    end
end
%str1=['C:\Users\ajith\Desktop\integrated\hog',num2str(m1),'.txt'];
%save(str1,'H3','-ASCII');
hold off;
%}
%enter two hand code after hog...till hog its already there
countd=1;
for f=178:179
str=['C:\Users\valli\Dropbox\Photos\project common\test\a\snapx',num2str(f),'.jpg'];%change path over here.
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

%quiver(grad_yu,grad_xr);
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
       % disp(H2);
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
pause(5);
load net;
load net1;

y1a=sim(net,H);
display(y1a);
    
    fingerfound=max(y1a);
    disp(fingerfound);
    
y22=sim(net1,H);
display(y22);
    
    fingerfound1=max(y22);
    disp(fingerfound1);
disp('h');
for f11=1:15
  switch fingerfound
   case y1a(1)
      disp('letter Neural a');
      neta=1;
      str=['C:\Users\valli\Dropbox\Photos\project common\train1\a\',num2str(f11),'.jpg'];
         case y1a(2)
            disp('letter Neural d');
            neta=2;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\d\',num2str(f11),'.jpg'];
             case y1a(3)
            disp('letter Neural g');
            neta=3;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\g\',num2str(f11),'.jpg'];
             case y1a(4)
            disp('letter Neural m');
            neta=4;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\m\',num2str(f11),'.jpg'];
             case y1a(5)
            disp('letter Neural s');
            neta=5;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\s\',num2str(f11),'.jpg'];
             case y1a(6)
            disp('letter Neural t');
            neta=6;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\t\',num2str(f11),'.jpg'];
            case y1a(7)
            disp('letter Neural z');
            neta=7;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\z\',num2str(f11),'.jpg'];
                case y1a(8)
            disp('letter Neural b');
                        neta=8;
                        str=['C:\Users\valli\Dropbox\Photos\project common\train1\b\',num2str(f11),'.jpg'];
      case y1a(9)
          disp('letter Neural x');
                      neta=9;
                      str=['C:\Users\valli\Dropbox\Photos\project common\train1\x\',num2str(f11),'.jpg'];
      case y1a(10)
          disp('letter Neural r');
                      neta=10;
                      str=['C:\Users\valli\Dropbox\Photos\project common\train1\r\',num2str(f11),'.jpg'];
      case y1a(11)
          disp('letter Neural p');
                      neta=11;
                      str=['C:\Users\valli\Dropbox\Photos\project common\train1\p\',num2str(f11),'.jpg'];
      case y1a(12)
          disp('letter Neural n');
                      neta=12;
                      str=['C:\Users\valli\Dropbox\Photos\project common\train1\n\',num2str(f11),'.jpg'];
      case y1a(13)
          disp('letter Neural o');
                      neta=13;
                      str=['C:\Users\valli\Dropbox\Photos\project common\train1\o\',num2str(f11),'.jpg'];
      case y1a(14)
          disp('letter Neural h');
                      neta=14;
          str=['C:\Users\valli\Dropbox\Photos\project common\train1\h\',num2str(f11),'.jpg'];
          otherwise
            disp('letter Neural nothing ');
    
  end
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
Hh=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
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

%quiver(grad_yu,grad_xr);
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
       % disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        Hh((cont-1)*B+1:cont*B,1)=H2;
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
x3(index,countd)=Hh(index);

end
for i= 1:2
    
        t3(i,countd)=0;
    
end
t3(1,countd)=1;
countd=countd+1;
end

for f11=1:15

  switch fingerfound1
   case y22(1)
      disp('letter Neural y');
      netb=1;
      str=['C:\Users\valli\Dropbox\Photos\project common\train1\y\',num2str(f11),'.jpg'];
         case y22(2)
            disp('letter Neural k');
            netb=2;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\k\',num2str(f11),'.jpg'];
             case y22(3)
            disp('letter Neural f');
            netb=3;
            str=['C:\Users\valli\Dropbox\Photos\project common\train1\f\',num2str(f11),'.jpg'];
      case y22(4)
          disp('letter Neural e');
          netb=4;
          str=['C:\Users\valli\Dropbox\Photos\project common\train1\e\',num2str(f11),'.jpg'];
      case y22(5)
          disp('letter Neural n');
          netb=5;
          str=['C:\Users\valli\Dropbox\Photos\project common\train1\n\',num2str(f11),'.jpg'];
      case y22(6)
          disp('letter Neural s');
          netb=6;
          str=['C:\Users\valli\Dropbox\Photos\project common\train1\s\',num2str(f11),'.jpg'];
      case y22(7)
          disp('letter Neural q');
          netb=7;
          str=['C:\Users\valli\Dropbox\Photos\project common\train1\q\',num2str(f11),'.jpg'];
          otherwise
            disp('letter Neural nothing ');
    
  end
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
Hh=zeros(nwin_x*nwin_y*B,1); % column ector with zeros
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

%quiver(grad_yu,grad_xr);
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
       % disp(H2);
        %disp('Hello');
       %disp(H2);         
        H2=H2/(norm(H2)+0.01); 
        %disp(H2);
        Hh((cont-1)*B+1:cont*B,1)=H2;
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
x3(index,countd)=Hh(index);

end
for i= 1:2
    
        t3(i,countd)=0;
    
end
t3(2,countd)=1;
countd=countd+1;
    
end



 net3 = network;
        net3=init(net3);
        net3.numInputs=1;
        net3.numLayers=3;
        net3.inputConnect(1,1) = 1;
%net.inputConnect(2,1) = 1;
%net.inputConnect(2,2) = 1;
         %net.inputConnect: [0 0 0; 1 0 0; 0 1 0]
        net3.biasConnect(1)=1;
        net3.biasConnect(2)=1;
        net3.biasConnect(3)=1;
        net3.layerConnect = [0 0 0; 1 0 0; 0 1 0];
        net3.outputConnect = [0 0 1];
        net3.inputs{1}.size=108;
         net3.layers{1}.size=21;
             net3.layers{2}.size=21;
  %               net.layers{3}.size=10;
        net3.layers{1}.transferFcn='tansig';
        net3.layers{1}.initFcn='initnw';
            net3.layers{2}.transferFcn='tansig';
%                  net.inputs{2}.size=10;
        net3.layers{2}.initFcn='initnw';
 %            net.inputs{3}.size=10;
        net3.layers{3}.transferFcn='tansig';
        net3.layers{3}.initFcn='initnw';
        
        net3.initFcn='initlay';
        net3.performFcn = 'mse';
        net3.trainFcn = 'trainlm';
        
%net = patternnet(10);
%net = feedforwardnet(25);
%net = cascadeforwardnet(10);
net3 = train(net3,x3,t3);


[net3,tr] = train(net3,x3,t3);
view(net3)
y223=sim(net3,H);
display(y223);
    
    fingerfound1=max(y223);
    disp(fingerfound1);
    if y223(1)==fingerfound1
        
  switch neta
   case 1
      disp('decision letter Neural a');
      
         case 2
            disp('decision letter Neural d');
           
             case 3
            disp('decision letter Neural g');
           
            
             case 4
            disp('decision letter Neural m');
           
             case 5
            disp('decision letter Neural s');
           
             case 6
            disp('decision letter Neural t');
   
      case 7
            disp('decision letter Neural z');
           
                case 8
            disp('decision letter Neural b');
                       
      case 9
          disp('decision letter Neural x');
                      
      case 10
          disp('decision letter Neural r');
                     
      case 11
          disp('decision letter Neural p');
                     
      case 12
          disp('decision letter Neural n');
                     
      case 13
          disp('decision letter Neural o');
                      
      case 14
          disp('decision letter Neural h');
                    
          
          otherwise
            disp('letter Neural nothing neta');
    
      end

    end
    if y223(2)==fingerfound1
           switch netb
   case 1
      disp('decision letter Neural y');
     
         case 2
            disp('decision letter Neural k');
           
             case 3
            disp('decision letter Neural f');
           
      case 4
          disp('decision letter Neural e');
         
      case 5
          disp('decision letter Neural n');
          
      case 6
          disp('decision letter Neural s');
         
      case 7
          disp('decision letter Neural q');
          
          otherwise
            disp('letter Neural nothing in netb ');
    
  end
        
    end
end
%end