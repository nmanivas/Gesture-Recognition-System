function singlehand(allfing,alldeg,fingdis,tnp,arr,loc,snapshot1,m1,background,left,deg)
    out=loc;
    points=[];
    c=[];
    hc=0;
    fc=0;
    arr2=0;
    of=0;
    adfing=0;
    tfx=0;
    tfy=0;
    %tfingdis=0;
    finfdiff=0;
    fbin(1)=0;
    fbin(2)=0;
    fbin(3)=0;
    fbin(4)=0;
    fbin(5)=0;
    fingmat=0;
     degg(1)=0;
    degg(2)=0;
    degg(3)=0;
    degg(4)=0;
    degg(5)=0;
    fl=0;
    feature=0;
    hc=1;
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
 %       str=['C:\Users\Sony\Desktop\codes\latestcodes\ullang',num2str(i),'.jpg'];   
%imwrite(ullang,str,'jpg');
        fingers=b-ullang;
  %      str=['C:\Users\Sony\Desktop\codes\latestcodes\fing',num2str(i),'.jpg'];   
%imwrite(fingers,str,'jpg');
%figure(3),imshow(fingers);
fingers_l=bwlabel(fingers);
maxfingers=max(max(fingers_l));
    ullangarea_est=(2*maxd-1)^2;
    minfingerarea=ullangarea_est/16;
    maxfingerarea=ullangarea_est/4;
    %fc=0;
    for j=1:maxfingers
      
        finger=fingers_l==j;
        area=sum(sum(finger));
        if(area>minfingerarea && area<maxfingerarea)
            r=regionprops(finger,'Centroid','MajorAxisLength','Orientation');
            cent=r.Centroid;
            majorl=r.MajorAxisLength;
            majorp1=majorl*(1/3);
            majorp2=majorl*(2/3);
            
            ang=r.Orientation;
            x=majorl/2*cosd(ang);
            y=majorl/2*sind(ang);
            p1=majorp1/2*cosd(ang);
            q1=majorp1/2*sind(ang);
            p2=majorp2/2*cosd(ang);
            q2=majorp2/2*sind(ang);
            x1=cent(1)+x;
            x2=cent(1)-x ;
            y1=cent(2)-y;
            y2=cent(2)+y;
            x3=cent(1)+p1;
            x4=cent(1)-p1;
            y3=cent(2)-q1;
            y4=cent(2)+q1;
            x5=cent(1)+p2;
            x6=cent(1)-p2;
            y5=cent(2)-q2;
            y6=cent(2)+q2;
             dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
             xx2=(dist/3);
            if(sum((c-[x1 y1]).^2) < sum((c-[x2 y2]).^2))
                    
   fc=fc+1;                 p=[x2 y2 x1 y1];
            else
                
                    
           fc=fc+1;         p=[x1 y1 x2 y2];
            end
    %plot(cent(1),cent(2),'r*');
  

                points=[points; p];
                %display(i);
               
               %display(fc);
               dist2=sqrt((p(3)-p(1))*(p(3)-p(1))+(p(4)-p(2))*(p(4)-p(2)));
              %display(dist2);
              arr1(fc)=majorl;
              if(fc==1)
                  tfx(fc)=p(3);
                  tfy(fc)=p(4);
                  tfingdis(fc)=majorl;
                  of=fc;
                  v1=[c(1),c(2)]-[p(3),p(4)]; 
                  v2=[c(1),c(2)]-[c(1)+50,c(2)];
                  a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
                  angleout = abs((a1>pi/2)*pi-a1);
                  degg(1)=angleout*(180/pi);
              end
              if(fc==2)
                  tfx(fc)=p(3);
                  tfy(fc)=p(4);
                  tfingdis(fc)=majorl;
                  of=fc;
                   v1=[c(1),c(2)]-[p(3),p(4)]; 
                  v2=[c(1),c(2)]-[c(1)+50,c(2)];
                  a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
                  angleout = abs((a1>pi/2)*pi-a1);
                  degg(2)=angleout*(180/pi);
              end
              if(fc==3)
                  tfx(fc)=p(3);
                  tfy(fc)=p(4);
                  tfingdis(fc)=majorl;
                  of=fc;
                   v1=[c(1),c(2)]-[p(3),p(4)]; 
                  v2=[c(1),c(2)]-[c(1)+50,c(2)];
                  a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
                  angleout = abs((a1>pi/2)*pi-a1);
                  degg(3)=angleout*(180/pi);
              end
              if(fc==4)
                  tfx(fc)=p(3);
                  tfy(fc)=p(4);
                  tfingdis(fc)=majorl;
                  of=fc;
                  v1=[c(1),c(2)]-[p(3),p(4)]; 
                  v2=[c(1),c(2)]-[c(1)+50,c(2)];
                  a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
                  angleout = abs((a1>pi/2)*pi-a1);
                  degg(4)=angleout*(180/pi);
              end
              if(fc==5)
                  tfx(fc)=p(3);
                  tfy(fc)=p(4);
                  tfingdis(fc)=majorl;
                  of=fc;
                    v1=[c(1),c(2)]-[p(3),p(4)]; 
                  v2=[c(1),c(2)]-[c(1)+50,c(2)];
                  a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
                  angleout = abs((a1>pi/2)*pi-a1);
                  degg(5)=angleout*(180/pi);
              end
              %{
arr2(fc)=arr1(fc);
              if(fc<=5)
               for k=1:fc
                  arr2(k)=arr2(k)-arr(k);
                  if(arr2(k)<0)
                      arr2(k)=-1*arr2(k);
                  end
              end
              m=min(arr2);
              for k=1:fc
                  if(arr2(k)==m)
                      mc=k;
                 end
                  
              end
              %}
              
        if(fc<=5)   
        if((arr1(fc)<=((3/4)*arr(fc)))||(fc==5))
            plot(cent(1),cent(2),'b*');
            plot(p(1),p(2),'b*');
            plot(p(3),p(4),'g*');
            line([x1 x2],[y1 y2],'linewidth',2,'color',[0 .5 1]); 
            
               
              else
                 plot(x3,y3,'y*');
                plot(x4,y4,'y*');

               plot(p(1),p(2),'b*');
                plot(p(3),p(4),'g*');
                line([x1 x2],[y1 y2],'linewidth',2,'color',[0 .5 1]);  
        end
        end
             
        end
    end
    end
    end
    %{
    str=['C:\Users\valli\books\neural\snap',num2str(i),'.jpg'];   
imwrite(d,str,'jpg');
    disp(i);
    %}
    of=of-1;
    for i=1:of
        %disp(i);
        %disp(of+1);
        adfing(i)=sqrt(((tfx(of+1)-tfx(i))*(tfx(of+1)-tfx(i)))+((tfy(of+1)-tfy(i))*(tfy(of+1)-tfy(i))));
    end
disp(adfing);

%index_range=[1:5];
y=1;
z=1;
x=1;
finger1=1;
finger2=1;
fbin(1)=0;
fbin(2)=0;
fbin(3)=0;
fbin(4)=0;
fbin(5)=0;

fingmat=0;
    for x=1:of
        for y=1:5
            for z=1:5
                %disp(allfing(j,k));
                %disp(fingdis(i));
            fingdiff(y,z)=abs(allfing(y,z)-adfing(x));
            %disp(fingdiff(j,k));
            end
        end
        %disp(fingdiff);
        minval=min(min(fingdiff(:)));
        %disp(minval);
        %disp(index_range);
        %index_range=[k:5];
        for a=1:5
            for b=1:5
                if(minval==fingdiff(a,b))
                    finger1=a;
                    finger2=b;
                end
            end
        end
        %index_range=[finger2:5];
        fbin(finger1)=1;
        fbin(finger2)=1;
        %str1=['finger ',num2str(finger1),' and ',num2str(finger2),' is open'];
        %disp(index_range);
        %disp(finger1);
        %disp(finger2);
        %disp(str1);
       %disp(fingmat);
    end 
    y=0;
    x=1;
    for x=1:length(fbin)
        if(fbin(x)==0)
            fingmat(x)=0; 
        else
            y=y+1;
            fingmat(x)=tfingdis(y);
        end
    end
   %disp(fingmat);
    x=1;
    for x=1:length(fingmat)
        if(fingmat(x)==0)
            str=['finger',num2str(x),'is closed'];
            feature(x)=0;
            fl(x)=0;
            disp(str);
        else
            if(fingmat(x)<=((3/4)*arr(x)))
                str=['finger',num2str(x),'is half open'];
                feature(x)=1;
                fl(x)=1;
                disp(str);
            else
                str=['finger',num2str(x),'is full open'];
                feature(x)=2;
                fl(x)=2;
                disp(str);
        end
    end
    end
   %disp(fl);
    if((fbin(2)==1)&&(fbin(3)==1)&&(fbin(4)==1)&&(fbin(1)==0)&&(fbin(5)==0))
        disp('Letter is W');
        disp(abs(deg(2)-deg(3)))
        feature(7)=abs(deg(2)-deg(3));
        disp(abs(deg(3)-deg(4)))
        feature(8)=abs(deg(3)-deg(4));
    end
    %if(tnp==4)
     if((fbin(3)==1)&&(fbin(4)==1)&&(fbin(2)==0)&&(fbin(1)==0)&&(fbin(5)==0))
         disp(abs(deg(3)-deg(4)));
         feature(8)=abs(deg(3)-deg(4));
        disp('Letter is V');
     end
    %else
      %  if(tnp==1)
        if((fbin(2)==1)&&(fbin(3)==1)&&(fbin(4)==0)&&(fbin(1)==0)&&(fbin(5)==0))
            disp(abs(deg(2)-deg(3)));
            feature(7)=abs(deg(2)-deg(3));
        disp('Letter is V');  
        end
        %end   
    %end
    %if(tnp==4)
     if((fbin(4)==1)&&(fbin(5)==1)&&(fbin(1)==0)&&(fbin(2)==0)&&(fbin(3)==0))
     if((fl(4)==2)&&(fl(5)==2)&&(fl(1)==0)&&(fl(2)==0)&&(fl(3)==0))
         if((abs(deg(1)-deg(2))>=50)||(abs(deg(1)-deg(2))<=60))
                disp(abs(deg(1)-deg(2)));
                feature(6)=abs(deg(1)-deg(2));
                disp('Letter C');
            else
        disp(abs(deg(4)-deg(5)));
        feature(9)=abs(deg(4)-deg(5));
        disp('Letter is L');
         end
         
     elseif((fl(4)==1)&&(fl(5)==2)&&(fl(1)==0)&&(fl(2)==0)&&(fl(3)==0))
        disp(abs(deg(4)-deg(5)));
        feature(9)=abs(deg(4)-deg(5));
        disp('Letter is U');
     end
     end
    %else
     %   if(tnp==1)
        if((fbin(1)==1)&&(fbin(2)==1)&&(fbin(4)==0)&&(fbin(3)==0)&&(fbin(5)==0))
         if((fl(1)==2)&&(fl(2)==2)&&(fl(4)==0)&&(fl(3)==0)&&(fl(5)==0))
             if((abs(deg(1)-deg(2))>=50)||(abs(deg(1)-deg(2))<=60))
                disp(abs(deg(1)-deg(2)));
                 feature(6)=abs(deg(1)-deg(2));
                disp('Letter C');
            else
        disp(abs(deg(1)-deg(2)));
            feature(6)=abs(deg(1)-deg(2));
        disp('Letter is L');
             end
        elseif((fl(1)==2)&&(fl(2)==1)&&(fl(4)==0)&&(fl(3)==0)&&(fl(5)==0))
            
                disp(abs(deg(1)-deg(2)));
                    feature(6)=abs(deg(1)-deg(2));
                disp('Letter is U');
        
         end  
        end
        if(fc==1)
            disp(fc);
                feature(2)=2;
            disp('letter I');
        elseif(fc==0)
            disp(fc);
            disp('letter O');
        end
        
        %end
    %end
    feature(10)=left;
    load nets;
    y1=sim(nets,[feature(1) feature(2) feature(3) feature(4) feature(5) feature(6) feature(7) feature(8) feature(9) feature(10)]');
   
    disp(y1);
    
    fingerfound=max(y1);
    disp(feature);
    
      switch fingerfound
   case y1(1)
      disp('letter Neural L');
      
         case y1(2)
            disp('letter Neural C');
            
             case y1(3)
            disp('letter Neural I');
            
             case y1(4)
            disp('letter Neural U');
            
             case y1(5)
            disp('letter Neural W');
            
             case y1(6)
            disp('letter Neural V');
            
          otherwise
            disp('letter Neural nothing ');
    
      end
        
    
    
end
%stop(vid);
        
        
        