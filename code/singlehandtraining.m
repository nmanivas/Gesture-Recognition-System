function [allfing alldeg fingdis tnp arr background left deg]=singlehandtraining()
for i=1:10
    tfingdis(i)=0;
end
finger1x=0;

 disp('Starting Video Device');
vid=videoinput('winvideo',1,'YUY2_320x240');
triggerconfig(vid,'manual');
start(vid);
starttime=tic;
loc=[-1 -1];
u11=0;
flaag=0;
%imwrite(Frame, 'F:\image%04d.jpg',x);

 %str=['F:\',num2str(i),'.jpg'];   
%imwrite(snapshot,str,'jpg');
%end
for i=1:10
snap=ycbcr2rgb(getsnapshot(vid));
snap1=snap;
figure(1),imshow(snap1);
end
background=snap1;
background=rgb2gray(background);

X=0;
Y=0;
zz=0;
yy=0;

X1=0;
Y1=0;

Xmax=0;
Ymax=0;
figure(1)
for m1=1:20
    %loc=[-1 -1];
figure(1)
h_mean=130/255;
h_var=0;
k=0;
count=0;
i=1;
cent=0;
cnt1=0;
arr=0;
c1=0;
c2=0;
fing1dis=0;
fing2dis=0;
fing3dis=0;
fing4dis=0;
fing5dis=0;
allfing=0;
left=0;
%TRAINING
%for i=1:30
   % snapshot=ycbcr2rgb(getsnapshot(vid));
  
   % clearvars -except background vid m1 zz yy X Y X1 Y1 Xmax Ymax u11 flaag;
     snapshot1=ycbcr2rgb(getsnapshot(vid));
 snap=snapshot1;   
    snapshot1=imfilter(snapshot1,fspecial('average',[5,5]));
    snapshot1=imadjust(snapshot1,[],[0.1,1]);
    set(0,'currentfigure',1);
 %   figure(m1),imshow(snapshot1);
str=['C:\Users\valli\Videos\Desktop\integrated\snapshot1',num2str(m1),'.jpg'];   
imwrite(snapshot1,str,'jpg');
snapshot1=segmentation(snapshot1,m1,snap,background);  %segmentation call

     figure(10),imshow(snapshot1);
     snapshot1=rgb2gray(snapshot1);
     thres=graythresh(snapshot1);
     snapshot1=im2bw(snapshot1,thres);
     skin=snapshot1;
     skin_l=bwlabel(snapshot1);
    s=skin_l(skin_l~=0);
    [maxareai maxarea]=mode(s(:));
    b=skin_l==maxareai;
  out=loc;
    points=[];
    c=[];
    hc=0;
    if maxareai>0
    if maxarea/numel(skin)>0.05
        hc=1;
        d=bwdist(~b,'euclidean');
    figure(2),imshow(d,[]);
        maxd=max(max(d));
        hold on
        ind=find(maxd==d);
        [cy cx]=ind2sub(size(d),ind);
        c=mean([cx cy],1);
        plot(c(1),c(2),'c*');
         out=-[c(1),c(2)];
        d_thresh=double(uint32(maxd/3));
        ullang=d>d_thresh;
        ullang=imdilate(ullang,strel('disk',d_thresh,0));
        fingers=b-ullang;
%figure(3),imshow(fingers);
fingers_l=bwlabel(fingers);
maxfingers=max(max(fingers_l));
    ullangarea_est=(2*maxd-1)^2;
    minfingerarea=ullangarea_est/16;
    maxfingerarea=ullangarea_est/4;
    fc=0;
     for j=1:maxfingers
      
        finger=fingers_l==j;
        area=sum(sum(finger));
        if(area>minfingerarea && area<maxfingerarea)
            r=regionprops(finger,'Centroid','MajorAxisLength','Orientation');
            cent=r.Centroid;
            majorl=r.MajorAxisLength;
            ang=r.Orientation;
            x=majorl/2*cosd(ang);
            y=majorl/2*sind(ang);
            x1=cent(1)+x;
            x2=cent(1)-x;
            y1=cent(2)-y;
            y2=cent(2)+y;
             dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
             xx2=(dist/3);
            if(sum((c-[x1 y1]).^2) < sum((c-[x2 y2]).^2))
                fc=fc+1;               
                p=[x2 y2 x1 y1];
            else
                            
           fc=fc+1;        
           p=[x1 y1 x2 y2];
            end
            count=fc;
             plot(0,0,'r*');
    plot(0,50,'g*');
    plot(50,0,'b*');
             plot(cent(1),cent(2),'r*');
             plot(p(1),p(2),'y*');
                plot(p(3),p(4),'g*');
                line([x1 x2],[y1 y2],'linewidth',2,'color',[0 .5 1]);
                points=[points; p];
               
                   if(fc==2)
                   
                dist1=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
                %display(dist1);
                
                r1=dist1;
               % f2=f2+1;
                %fing2dis(f2)=dist1;
                
                f2x=p(3);
                f2y=p(4);
                end
               
               
                   if(fc==4) 
                dist1=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
                %display(dist1);f
                
               i1=dist1;
               %f4=f4+1;
               %fing4dis(f4)=dist1;
                f4x=p(3);
                f4y=p(4);
                
                   end
               
                
                   if(fc==1)
                   
                dist3=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
                %display(dist3);
                l1=dist3;
                %f1=f1+1;
               %fing1dis(f1)=dist3;
                f1x=p(3);
                f1y=p(4);
                   end
                   if(fc==3)
                   
                dist3=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
                %display(dist3);
                
                 m1=dist3;
                 %f3=f3+1;
               %fing3dis(f3)=dist3;
                f3x=p(3);
                f3y=p(4);
                   end
                
               
                   if(fc==5)
                dist1=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
                %display(dist1);
                t1=dist1;
                %f5=f5+1;
               %fing5dis(f5)=dist1;
                
                f5x=p(3);
                f5y=p(4);
                   end
                  
              
        end
     end
    end
    end
    if(count==5)
    cnt1=cnt1+1;
   fing1dis=fing1dis+l1; 
   %display(little);
   fing3dis=fing3dis+m1;
   fing4dis=fing4dis+i1; 
   fing5dis=fing5dis+t1;    
   fing2dis=fing2dis+r1;  
   c1=c1+1;
   c2=c2+1;
   finger1x(c1)=f1x;
   finger1y(c2)=f1y;
   finger4x(c1)=f4x;
   finger4y(c2)=f4y;
   finger2x(c1)=f2x;
   finger2y(c2)=f2y;   
   finger5x(c1)=f5x;
   finger5y(c2)=f5y; 
   finger3x(c1)=f3x;
   finger3y(c2)=f3y;
    end
    
end
%stop(vid);
%fing1fin=sum(fing1dis)/cnt1;
%fing2fin=sum(fing2dis)/cnt1;
%fing3fin=sum(fing3dis)/cnt1;
%fing4fin=sum(fing4dis)/cnt1;
%fing5fin=sum(fing5dis)/cnt1;
%fingfin=[fing1fin fing2fin fing3fin fing4fin fing5fin];
f1xval=sum(finger1x)/length(finger1x);
f2xval=sum(finger2x)/length(finger2x);
f3xval=sum(finger3x)/length(finger3x);
f4xval=sum(finger4x)/length(finger4x);
f5xval=sum(finger5x)/length(finger5x);
f1yval=sum(finger1y)/length(finger1y);
f2yval=sum(finger2y)/length(finger2y);
f3yval=sum(finger3y)/length(finger3y);
f4yval=sum(finger4y)/length(finger4y);
f5yval=sum(finger5y)/length(finger5y);
finx=[f1xval f2xval f3xval f4xval f5xval];
finy=[f1yval f2yval f3yval f4yval f5yval];
thumy=max(finy);
v1=[f1xval,f1yval]-[c(1),c(2)]; 
v2=[c(1)+50,c(2)]-[c(1),c(2)];
a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
deg(1)=angleout*(180/pi);
v1=[f2xval,f2yval]-[c(1),c(2)]; 
v2=[c(1)+50,c(2)]-[c(1),c(2)];
a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
deg(2)=angleout*(180/pi);
v1=[f3xval,f3yval]-[c(1),c(2)]; 
v2=[c(1)+50,c(2)]-[c(1),c(2)];
a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
deg(3)=angleout*(180/pi);
v1=[f4xval,f4yval]-[c(1),c(2)]; 
v2=[c(1)+50,c(2)]-[c(1),c(2)];
a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
deg(4)=angleout*(180/pi);
v1=[f5xval,f5yval]-[c(1),c(2)]; 
v2=[c(1)+50,c(2)]-[c(1),c(2)];
a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
deg(5)=angleout*(180/pi);
%{
if(thumy==f1yval)
    display('Finger 1 is thumb'); 
end
if(thumy==f5yval)
    display('Finger 5 is thumb')
end
%}
fingdis(1)=sqrt(((f1xval-f2xval)*(f1xval-f2xval))+((f1yval-f2yval)*(f1yval-f2yval)));
fingdis(2)=sqrt(((f2xval-f3xval)*(f2xval-f3xval))+((f2yval-f3yval)*(f2yval-f3yval)));
fingdis(3)=sqrt(((f3xval-f4xval)*(f3xval-f4xval))+((f3yval-f4yval)*(f3yval-f4yval)));
fingdis(4)=sqrt(((f4xval-f5xval)*(f4xval-f5xval))+((f4yval-f5yval)*(f4yval-f5yval)));
thumd=max(fingdis);
for i=1:4
    if(fingdis(i)==thumd)
        tnp=i;
    end
end
fing1dis=fing1dis/cnt1;
    fing2dis=fing2dis/cnt1;
    fing3dis=fing3dis/cnt1;
    fing4dis=fing4dis/cnt1;
    fing5dis=fing5dis/cnt1;
    arr=[fing1dis fing2dis fing3dis fing4dis fing5dis];
    display(arr);
    display(tnp);
    display(fingdis);
if(tnp==1)
    display('Thumb is finger 1');
    left=0;
   
end
if(tnp==4)
    display('Thumb is finger 5');
left=1;
end
%Contstraint Matrix building
allfing(1,1)=1000;
allfing(1,2)=sqrt(((f1xval-f2xval)*(f1xval-f2xval))+((f1yval-f2yval)*(f1yval-f2yval)));
allfing(1,3)=sqrt(((f1xval-f3xval)*(f1xval-f3xval))+((f1yval-f3yval)*(f1yval-f3yval)));
allfing(1,4)=sqrt(((f1xval-f4xval)*(f1xval-f4xval))+((f1yval-f4yval)*(f1yval-f4yval)));
allfing(1,5)=sqrt(((f1xval-f5xval)*(f1xval-f5xval))+((f1yval-f5yval)*(f1yval-f5yval)));
allfing(2,1)=1000;
allfing(2,2)=1000;
allfing(2,3)=sqrt(((f2xval-f3xval)*(f2xval-f3xval))+((f2yval-f3yval)*(f2yval-f3yval)));
allfing(2,4)=sqrt(((f2xval-f4xval)*(f2xval-f4xval))+((f2yval-f4yval)*(f2yval-f4yval)));
allfing(2,5)=sqrt(((f2xval-f5xval)*(f2xval-f5xval))+((f2yval-f5yval)*(f2yval-f5yval)));
allfing(3,1)=1000;
allfing(3,2)=1000;
allfing(3,3)=1000;
allfing(3,4)=sqrt(((f3xval-f4xval)*(f3xval-f4xval))+((f3yval-f4yval)*(f3yval-f4yval)));
allfing(3,5)=sqrt(((f3xval-f5xval)*(f3xval-f5xval))+((f3yval-f5yval)*(f3yval-f5yval)));
allfing(4,1)=1000;
allfing(4,2)=1000;
allfing(4,3)=1000;
allfing(4,4)=1000;
allfing(4,5)=sqrt(((f4xval-f5xval)*(f4xval-f5xval))+((f4yval-f5yval)*(f4yval-f5yval)));
allfing(5,1)=1000;
allfing(5,2)=1000;
allfing(5,3)=1000;
allfing(5,4)=1000;
allfing(5,5)=1000;
display(allfing);

alldeg(1,1)=1000;
alldeg(1,2)=abs(deg(1)-deg(2));
alldeg(1,3)=abs(deg(1)-deg(3));
alldeg(1,4)=abs(deg(1)-deg(4));
alldeg(1,5)=abs(deg(1)-deg(5));
alldeg(2,1)=1000;
alldeg(2,2)=1000;
alldeg(2,3)=abs(deg(2)-deg(3));
alldeg(2,4)=abs(deg(2)-deg(4));
alldeg(2,5)=abs(deg(2)-deg(5));
alldeg(3,1)=1000;
alldeg(3,2)=1000;
alldeg(3,3)=1000;
alldeg(3,4)=abs(deg(3)-deg(4));
alldeg(3,5)=abs(deg(3)-deg(5));
alldeg(4,1)=1000;
alldeg(4,2)=1000;
alldeg(4,3)=1000;
alldeg(4,4)=1000;
alldeg(4,5)=abs(deg(4)-deg(5));
alldeg(5,1)=1000;
alldeg(5,2)=1000;
alldeg(5,3)=1000;
alldeg(5,4)=1000;
alldeg(5,5)=1000;
   %display(cnt1);
    %display(little);
    %display(index);
    %display(thumb);
    %display(ring);
    %display(middle);
        %start(vid);
        display('display gestures');%end of config phase
        i=1;
        stop(vid);
        end
