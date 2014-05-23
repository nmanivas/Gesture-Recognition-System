%hand segmentation using skin detection and motion detection
%snapshot1=imread('C:\Users\Sony\Desktop\Image8.jpg');



countreg1=0;
countreg2=0;
X=0;
Y=0;
zz=0;
yy=0;

X1=0;
Y1=0;

Xmax=0;
Ymax=0;
flaag=0;

%[allfing alldeg fingdis tnp arr background left deg]=singlehandtraining();

Movie = VideoReader('C:\Users\valli\Videos\Desktop\New folder\Record Videos\sangiz.avi');
nframes = Movie.NumberOfFrames;
for m1 = 1 :10

snapshot = read(Movie,m1);
%str=['C:\Users\valli\Videos\Desktop\New folder\snap\ss',num2str(i),'.jpg'];   
%imwrite(snapshot,str,'jpg');
%disp('Starting Video Device');
%vid=videoinput('winvideo',1,'YUY2_320x240');
%triggerconfig(vid,'manual');
%start(vid);
%starttime=tic;
%loc=[-1 -1];
%u11=0;

%snap=ycbcr2rgb(snapshot);
snap=snapshot;
snap1=snap;
figure(1),imshow(snap1);

background=snap1;
X=0;
Y=0;
zz=0;
yy=0;

X1=0;
Y1=0;

Xmax=0;
Ymax=0;
background=rgb2gray(background);
end


for m1=160:300
     snapshot2 = read(Movie,m1);
str=['C:\Users\valli\Videos\Desktop\New folder\snap\sss',num2str(m1),'.jpg'];   
imwrite(snapshot2,str,'jpg');
    clearvars -except background vid m1 zz yy X Y X1 Y1 Xmax bw Ymax u11 snapshot2 Movie nframes flaag Ic m1 nettry3 netmn netrmns nethygz netgy netxe netxf
     %snapshot1=ycbcr2rgb(getsnapshot(vid));
     snapshot1=snapshot2;
 snap=snapshot1;   
    snapshot1=imfilter(snapshot1,fspecial('average',[5,5]));
    snapshot1=imadjust(snapshot1,[],[0.1,1]);
    set(0,'currentfigure',1);
 %   figure(m1),imshow(snapshot1);
str=['C:\Users\valli\Videos\Desktop\integrated\snapshot1',num2str(m1),'.jpg'];   
imwrite(snapshot1,str,'jpg');
snapshot1=segmentation(snapshot1,m1,snap,background);  %segmentation call

     figure(10),imshow(snapshot1);
     disp(m1);
hold on;


Ima=rgb2gray(snapshot1);
th=graythresh(Ima);
bw=im2bw(Ima,th);
%imshow(bw);
bz=bwlabel(bw);
b11=regionprops(bz,'all');
if(~isempty(b11))
len=length(b11);
%areas=[b11.Area];
%[maxArea largestBlobIndex]=max(areas);
i2=0;
for i2=1:len
x1=b11(i2).BoundingBox(1);
y1=b11(i2).BoundingBox(2);
x2=b11(i2).BoundingBox(3);
y2=b11(i2).BoundingBox(4);
rectangle('Position',[x1,y1,x2,y2],'Edgecolor','r');
end
hold off;
if(m1==1)
    if(len==2)
     flaag=1;
    end
end

if(len==1)
    if(flaag==1)
        Ic=imcrop(Ima,[x1 y1 x2 y2]);
        disp(m1);
        %twohand(Ic,m1,  net, net1, neta1, netxe, netbf, netgf1, netgk, nettk, nettf, neths1, netak, netdk, netaf, netdf, netxf, netms, netok, netoy, netas, netmn, nette, netzf, netaq, netay, netbn, netbs, netbq, netge, netgy, netgn, netgs, netgq, netoe, neton, netof, netos, netoq, netty, netts, nettn, nettq, netxy, netxk, netxs, netxn, netxq, netze, netzn, netzk, netby, netbk, netbe, netdy, netde, netdn, netds, netdq, netmk, netmy, netme, netmf, netmq, nethf, nethk, nethe, nethn, netpe, netpf, netpq, netps, netrk, netrf, netre, nethy, nethq); %two hand call
        %twohand(Ic,m1); %two hand call
        % twohand(Ic,m1,net,net1,netmn,netms,neta1,netaq,netay,netdk,netdf,netdq,netzq,netzf,netzk,nette,nettf,netrk);
    twohand(Ic,m1,bw,nettry3,netmn,netrmns);
    end
    if(flaag==0)
        disp('singlehand');%singlehand(bw);%enter single hand code after segmentation
    end
end
if(len==2)
     flaag=1;
 end
end
end




    
    


    
    