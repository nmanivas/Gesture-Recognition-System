function snapshot1=segmentation(snapshot1,m1,snap,background)
[H S V]=rgb2hsv(snapshot1);
    H=mod(H+128/255,1+1/255);
    H=medfilt2(H,[3 3]);
    skin=H>.4&H<.6 & S>.1&S<.9;
    skin=imfill(skin,'holes');
    skin_l1=bwlabel(skin);
blob1=regionprops(skin_l1);
numb1=size(blob1,1);
s1=skin_l1(skin_l1~=0);
[maxarea1i maxarea1]=mode(s1(:));
myarea1=round(maxarea1./10);
bb1=(skin_l1==maxarea1i);
maxi1(1)=maxarea1i;
k1=2;
flag1=0;
l1=1;
%figure(100),imshow(skin_l1)
for i1=1:numb1
    for j1=1:l1
        if(i1==maxi1(j1))
            flag1=1;
        end
    end  
    if(flag1==0)
    big1=blob1(i1).Area;
    if(big1>=myarea1)
        %disp(big);
        %disp(i);
        
        s1=s1(s1~=maxarea1i);
[maxarea1i maxarea1]=mode(s1(:));
maxi1(k1)=maxarea1i;
k1=k1+1;
l1=length(maxi1);
temp1=(skin_l1==maxarea1i);
bb1=bb1+temp1;
%disp(l);
%disp(maxi);
%disp('Hello');
%disp(maxarea1i);
%disp(maxarea1);


    end
    end
    flag1=0;
end
%figure(5),imshow(bb1);
str=['C:\Users\valli\Videos\Desktop\integrated\seg',num2str(m1),'.jpg'];   
imwrite(bb1,str,'jpg');
snap1=snap;
    Im=snap1;
    snap1=rgb2gray(snap1);
    I=imabsdiff(snap1, background);

 for i=1:240
    for j=1:320
        if(I(i,j)<50)
            I(i,j)=0;
        end
    end
end
%figure(10),imshow(I);

for i=1:240
    for j=1:320
        if(I(i,j)~=0)
            I(i,j)=255;
        end
    end
end

imfill(I,'holes');
%figure(5),imshow(I);

diff=graythresh(I);

im2=im2bw(I,diff);
%figure(4),imshow(I);
%figure(5),imshow(im2);


for i=1:240
    for j=1:320
        if(I(i,j)<25)
         snap(i,j,1)=0;
         snap(i,j,2)=0;
         snap(i,j,3)=0;
        end
    end
end
Im=snap;
%figure(2),imshow(Im);

str=['C:\Users\valli\Videos\Desktop\integrated\',num2str(m1+200),'.jpg'];   
imwrite(Im,str,'jpg');

g=rgb2gray(Im);
h=graythresh(g);
Imb=im2bw(g,0.1);
skin_l=bwlabel(Imb);
blob=regionprops(skin_l);
numb=size(blob,1);
s=skin_l(skin_l~=0);

[maxareai maxarea]=mode(s(:));
%myarea=round(maxarea./10);
myarea=2500;
b=(skin_l==maxareai);
maxi(1)=maxareai;
k=2;
flag=0;
l=1;
%figure(10),imshow(b);
%figure(100),imshow(skin_l)
for d=1:numb
    for e=1:l
        if(d==maxi(e))
            flag=1;
        end
    end  
    if(flag==0)
    big=blob(d).Area;
    if(big>=myarea)
        %disp(big);
        %disp(i);
        
        s=s(s~=maxareai);
[maxareai maxarea]=mode(s(:));
maxi(k)=maxareai;
k=k+1;
l=length(maxi);
temp=(skin_l==maxareai);
b=b+temp;
%disp(l);
%disp(maxi);
%disp('Hello');
%disp(maxareai);
%disp(maxarea);


    end
    end
    flag=0;
    
end
%figure(2),imshow(b);
str=['C:\Users\valli\Videos\Desktop\integrated\',num2str(m1+300),'.jpg'];   
imwrite(b,str,'jpg');

%b=logical(b);
%figure(2),imshow(b);

b=imresize(b,[240,320]);
b1=b;
    %background=snap1;
    q1=1;
    v1=0;
    w=0;
    fla=0;
    
for r=1:240
        for u=1:320
            if(bb1(r,u)==1 && b(r,u)==1)
                if(w~=0)
                for w1=1:w
                    if(v1(w1)==skin_l(r,u))
                        fla=1;
                    end
                end
                end
                if(fla==0)
              v1(q1)=skin_l(r,u);
              q1=q1+1;
              w=length(v1);
                end
            end
            fla=0;
        end
    end
   g2=0;
    for r1=1:240
        for u1=1:320
     for t1=1:w
            if(v1(t1)~=skin_l(r1,u1))
                g2=g2+1;
    
            end
    end
    if(g2==w)
    b(r1,u1)=0;
    end
    g2=0;
        end
    end
    
    q1=1;
    v1=0;
    w=0;
    fla=0;
    
    for r=1:240
        for u=1:320
            if(bb1(r,u)==1 && b1(r,u)==1)
                if(w~=0)
                for w1=1:w
                    if(v1(w1)==skin_l1(r,u))
                        fla=1;
                    end
                end
                end
                if(fla==0)
              v1(q1)=skin_l1(r,u);
              q1=q1+1;
              w=length(v1);
                end
            end
            fla=0;
        end
    end
   g2=0;
    for r1=1:240
        for u1=1:320
     for t1=1:w
            if(v1(t1)~=skin_l1(r1,u1))
                g2=g2+1;
    
            end
    end
    if(g2==w)
    bb1(r1,u1)=0;
    end
    g2=0;
        end
    end
    bb1=imfill(bb1,'holes');
    b=imfill(b,'holes');
    area=0;
     for i=1:240
    for j=1:320
        if(b(i,j)==1)
            area=area+1;
         
        end
    end
     end
     area1=0;
     for i=1:240
    for j=1:320
        if(bb1(i,j)==1)
            area1=area1+1;
         
        end
    end
     end
     
     if(area1>area)
         b=bb1;
     end
    
    
     for i=1:240
    for j=1:320
        if(b(i,j)==0)
            snapshot1(i,j,1)=0;
            snapshot1(i,j,2)=0;
         snapshot1(i,j,3)=0;
         
        end
    end
     end
str=['C:\Users\valli\Videos\Desktop\integrated\',num2str(m1+400),'.jpg'];   
    imwrite(b,str,'jpg');
str=['C:\Users\valli\Videos\Desktop\integrated\',num2str(m1+500),'.jpg'];   
    imwrite(snapshot1,str,'jpg');

end