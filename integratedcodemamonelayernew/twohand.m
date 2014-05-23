function          twohand(Ic,m1,bw,nettry3,netmn,netrmns)

if(m1<3)
    

end
%imshow(Ic);

 Im=Ic;
 BW=bw;
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

%load net;
%load net1;

y1a=sim(nettry3,H);
display(y1a);
    
    fingerfound=max(y1a);
    disp(fingerfound);
    


%for f11=1:15
  switch fingerfound
   case y1a(1)
      disp('letter Neural a');
      htxtins = vision.TextInserter('a');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
      neta=1;
         case y1a(2)
            disp('letter Neural d');
            htxtins = vision.TextInserter('d');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
            neta=2;
             case y1a(3)
            disp('letter Neural g');
            neta=3;
             case y1a(4)
            disp('letter Neural m');
                neta=4; 
             case y1a(5)
            disp('letter Neural t');
            neta=5;
            
            case y1a(6)
            disp('letter Neural z');
            neta=6;
  %           htxtins = vision.TextInserter('z');
%htxtins.Color = [255, 255, 255]; % [red, green, blue]
%htxtins.FontSize = 24;
%htxtins.Location = [15 15]; % [x y] [100 315];
%Ic = step(htxtins, Ic);
 %           imshow(Ic);
      case y1a(7)
          disp('letter Neural b');
          neta=7;
          htxtins = vision.TextInserter('b');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
      case y1a(8)
          disp('letter Neural x');
          neta=8;
          case y1a(9)
          disp('letter Neural r');
          neta=9;
      case y1a(10)
          disp('letter Neural p');
          neta=10;
          htxtins = vision.TextInserter('p');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
          case y1a(11)
          disp('letter Neural h');
          neta=11;
          case y1a(12)
      disp('letter Neural y');
      neta=12;
      %htxtins = vision.TextInserter('y');
%htxtins.Color = [255, 255, 255]; % [red, green, blue]
%htxtins.FontSize = 24;
%htxtins.Location = [15 15]; % [x y] [100 315];
%Ic = step(htxtins, Ic);
 %           imshow(Ic);
         case y1a(13)
            disp('letter Neural k');
            neta=13;
            htxtins = vision.TextInserter('k');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
             case y1a(14)
            disp('letter Neural f');
            neta=14;
            case y1a(15)
          disp('letter Neural e');
          neta=15;
           case y1a(16)
          disp('letter Neural n');
          neta=16;
          case y1a(17)
          disp('letter Neural s');
          neta=17;
          case y1a(18)
          disp('letter Neural q');
          neta=18;
          htxtins = vision.TextInserter('q');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
          otherwise
            disp('letter Neural nothing ');
    
      end

%end




mcount=0;
ncount=0;

if(neta==8||neta==15||neta==14||neta==5)
    %yae=sim(netxe,H);
%display(yae);
    
%BW=rgb2gray(Im);
%BW=im2bw(BW,0.1);
%hy = fspecial('sobel');
%BW = edge(Ima,'sobel');

%hx = hy';
%imshow(BW);
ii=1;jj=1;
counter1=0;
white=0;
icol=0;
black=0;
bflag=0;
wflag=0;
wcounter=0;
fflag=0;
counter=0;
%counter=0;
%Iy = imfilter(double(Ima), hy, 'replicate');
%Ix = imfilter(double(Ima), hx, 'replicate');
%gradmag = sqrt(Ix.^2 + Iy.^2);
 %imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
 %value1=regionprops(BW,'PixelList');
 %disp(value1);
 for ii=1:240
     
    for jj=1:320
 if((BW(ii,jj)==1))
     fflag=fflag+1;
     if(fflag==1)
     icol=ii;
     icol2=ii+20;
            icol3=ii+5;
     end
 end
    end
 end
    % counter=counter+1;
     %bflag=1;
     %wcounter=wcounter+1;
     %if(wcounter>3)
     %wflag=1;
     %end
   %  if(counter==1)
   %  for hk=1:5
         %for hj=1:320
             
            for hj=1:320
                if(BW(icol3,hj)==0&&bflag==0)
                    bflag=1;
                    wcounter=0;
                end
                if(BW(icol3,hj)==1&&bflag==1)
                    wcounter=wcounter+1;
                    if(wcounter>3)
                        wflag=1;
                        counter=counter+1;
                    end
                end
                if(BW(icol3,hj)==0&&bflag==1&&wflag==1)
                    bflag=2;
                    wcounter=0;
                end
                if(BW(icol3,hj)==1&&bflag==2)
                   wcounter=wcounter+1;
                   if(wcounter==1)
                       wflag=2;
                   end
                end
            end
            %BW(i+hk,hj)=1;
        % end
    % end
    
     
     for hk=1:5
         for hj=1:320
             BW(icol3+hk,hj)=1;
         end
     end
     disp(counter);
     if(wflag==2)
         disp('x');
         imshow(BW);
     end
     if(wflag~=2&&bflag==2&&counter>28)
         disp('t');
     else
    
     %if(counter==1)
     if(wflag~=2)
     for hj=1:320
         if(BW(icol2,hj)==1&&black==0)
             white=white+1;
         end
         if(BW(icol2,hj)==0&&white~=0)
             black=1;
             %disp(white);
         end
         
     end
     %line([1 240],[j j],'linewidth',2,'color',[0 .5 1]);
     %disp('1');
 
    
 
 imshow(BW);
 disp(white);
 if(white<25)
     disp('e');
 else
     disp('f');
 end
     end
     end
end
if(neta==3||neta==6||neta==11)
 %   BW=rgb2gray(Im);
%BW=im2bw(BW,0.1);
%hy = fspecial('sobel');
%BW = edge(Ima,'sobel');

%hx = hy';
%imshow(BW);
i=1;j=1;
countj=0;
jcol=0;
x=1;
y=1;
countjcol=0;
wflag=0;
bflag=0;
for x=1:240
    
count(x)=0;
    
end
for y=1:320
    county(y)=0;
end
%Iy = imfilter(double(Ima), hy, 'replicate');
%Ix = imfilter(double(Ima), hx, 'replicate');
%gradmag = sqrt(Ix.^2 + Iy.^2);
 %imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
 %value1=regionprops(BW,'PixelList');
 %disp(value1);
 for i=1:240
     
    for j=1:320
 if((BW(i,j)==1))
   %  fflag=fflag+1;
     count(i)=count(i)+1;
 end
    end
 end
     ele=max(count);
     disp(ele);
     for i=1:240
         if(count(i)==ele)
             maxind=i;
         end
     end
     %for mx=1:5
     %for j=1:320
      %   BW(maxind+mx,j)=1;
     %end
     %end
     imshow(BW);
     if(ele>215)
         disp('z');
    
     else
    for i=1:240
     
    for j=1:320
 if((BW(i,j)==1))
     county(j)=county(j)+1;
 end
    end
    end
    ele1=max(county);
    disp(ele1);
    for j=1:320
        if(county(j)==ele1)
            maxind1=j;
        end
    end
    %for i=1:240
     %   for mx=1:5
      %      BW(i,maxind1+mx)=1;
       % end
    %end
    imshow(BW);
    %if(ele1>130)
     %   disp('g');
    %end
    for j=1:320
        for i=1:240
            
            if(BW(i,j)==1)
               countj=countj+1;
                if(countj==1)
               jcol=j;
                end
            end
        end
    end
     jcol=jcol+30;
    for i=1:240
        if(BW(i,jcol)==0&&bflag==0)
            bflag=1;
        end
        if(BW(i,jcol)==1&&bflag==1)          
            countjcol=countjcol+1;
            if(countjcol>3)
            wflag=1;
            end
        end
        if(BW(i,jcol)==0&&wflag==1)
            bflag=2;
        end
    end
    disp(countjcol);
   
   for i=1:240
        for mx=1:5
            BW(i,jcol+mx)=1;
        end
   end
    imshow(BW);
    if(countjcol>70)
        disp('g');
    else
        disp('h');
    end
     end
end
if(neta==4||neta==9||neta==16||neta==17)
    yae1=sim(netrmns,H);
display(yae1);
    
    fingerfound1=max(yae1);
    disp(fingerfound1);
    switch fingerfound1
        case yae1(1)
            disp('letter decision Neural r');
             htxtins = vision.TextInserter('r');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
        case yae1(2)
            mcount=1;
        case yae1(3)
            ncount=1;
        case yae1(4)
            disp('letter decision neural s');
             htxtins = vision.TextInserter('s');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
    end
    if(mcount==1||ncount==1)
         yae2=sim(netmn,H);
display(yae2);
    
    fingerfound2=max(yae2);
    disp(fingerfound2);
    switch fingerfound2
        case yae2(1)
            disp('letter decision Neural m');
             htxtins = vision.TextInserter('m');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
        case yae2(2)
            disp('letter decision Neural n');
             htxtins = vision.TextInserter('n');
htxtins.Color = [255, 255, 255]; % [red, green, blue]
htxtins.FontSize = 24;
htxtins.Location = [15 15]; % [x y] [100 315];
Ic = step(htxtins, Ic);
            imshow(Ic);
    end
    end
end
str=['C:\Users\valli\Videos\Desktop\New folder\snappp\s',num2str(m1),'.jpg'];   
imwrite(Ic,str,'jpg');
 end