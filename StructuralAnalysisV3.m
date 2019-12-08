eleNodes=zeros(20000,2);
Ncounter=1;
i=1;
while 1
    if Ncounter>2
        disp('When the stucture is finished please Press Enter one more time to proceed.')
    else
    end
a=[num2str(i) '.th member starting point'];
    disp(a)
    prompt='\n' ;
    Inputs=input(prompt);
    if size(Inputs)~=[1 2]          %Structure bitince deðer girmeden enter'a bas
    break
end
        eleNodes(Ncounter,:)=[Inputs];
        
   a=[num2str(i) '.th member ending point'];
    disp(a)
    prompt='\n' ;
    Inputs=input(prompt);   
        eleNodes(Ncounter+1,:)=[Inputs];
        Ncounter=Ncounter+2;
        i=i+1;
end

NOE=(i-1);
eleNodes=eleNodes(1:NOE.*2,:);
Nodes=unique(eleNodes,'rows','stable');
NON=size(Nodes,1);

EI=[];
EA=[];
for i=1:NOE
    a=['Enter the EI of the ' ,num2str(i),'. member in kNm^2'];
    disp(a)
    prompt=('\n');
    ei=input(prompt);
    EI=[EI;ei];
    
     
end
    for i=1:NOE
    a=['Enter the EA of the ' ,num2str(i),'. member in kN'];
    disp(a)
    prompt=('\n');
    ea=input(prompt);
    EA=[EA;ea];
    
     
    end
    
    

i=1;
thetas=zeros(NOE,1);
for k=1:NOE
   thetas(k,1)=atand(( eleNodes(i+1,2)-eleNodes(i,2) ) ./ ( eleNodes(i+1,1)-eleNodes(i,1) ) );
   i=i+2;
end


DOFST=zeros(NON,3);
free=0;
for i=1:NON
    a=['For the ',num2str(i),'. Node'];
    disp(a)
prompt=('Is the movement in x-direction restricted? (Y/N)\n');
str=input(prompt,'s');
if str=='Y'||str=='y'
    DOFST(i,1)=-1;
elseif str=='n'||str=='N'
    DOFST(i,1)=0;
    free=free+1;
else
    disp('Error:Unexpected key')
end

prompt=('Is the movement in y-direction restricted? (Y/N)\n');
str=input(prompt,'s');
if str=='Y'||str=='y'
    DOFST(i,2)=-1;
elseif str=='n'||str=='N'
    DOFST(i,2)=0;
    free=free+1;
    else
    disp('Error:Unexpected key')
end
prompt=('Is the Rotational movement restricted? (Y/N)\n');
str=input(prompt,'s');
if str=='Y'||str=='y'
    DOFST(i,3)=-1;
elseif str=='n'||str=='N'
    DOFST(i,3)=0;
    free=free+1;
    else
    disp('Error:Unexpected key')
end
end

restr=(NON.*3)-free;

DOFS=transpose(DOFST);



counter=1;
    for k=1:NON
    for s=1:3
        if DOFS(s,k)==0
            DOFS(s,k)=counter;
            counter=counter+1;
        else
        end
   
    
    end
    end
        

    for k=1:NON
    for s=1:3
        if DOFS(s,k)==-1
            DOFS(s,k)=counter;
            counter=counter+1;
        else
        end
   
    
    end
    end
    


    DOFST=transpose(DOFS);
[Lia, Locb]=ismember(eleNodes,Nodes, 'rows');

k=1;
    for i=1:NOE

        LM(i,1:6)=[DOFST(Locb(k),:) DOFST(Locb(k+1),:)];
            k=k+2;
    end
           
K=zeros(NON.*3);    
    L=zeros(NOE,1);
    kl=zeros(6.*NOE,6);
    R=zeros(6.*NOE,6);
    for i=1:NOE
        L(i,1)= (    (  ( Nodes(    Locb((2.*i),1),1   )-Nodes(Locb(((2.*i)-1),1),1)).^2   )+ (  ( Nodes(    Locb((2.*i),1),2   )-Nodes(Locb(((2.*i)-1),1),2)).^2   ) ).^(1./2);
    kl((6.*i)-5:6.*i,1:6)=[
    EA(i)./L(i,1) 0 0 -EA(i)./L(i,1) 0 0;
    0 12.*EI(i)./(L(i,1).^3) 6.*EI(i)./(L(i,1).^2) 0 -12.*EI(i)./(L(i,1).^3) 6.*EI(i)./(L(i,1).^2);
    0 6.*EI(i)./(L(i,1).^2) 4.*EI(i)./L(i,1) 0 -6.*EI(i)./(L(i,1).^2) 2.*EI(i)./L(i,1);
    -EA(i)./L(i,1) 0 0 EA(i)./L(i,1) 0 0;
    0 -12.*EI(i)./(L(i,1).^3) -6.*EI(i)./(L(i,1).^2) 0 12.*EI(i)./(L(i,1).^3) -6.*EI(i)./(L(i,1).^2);
    0 6.*EI(i)./(L(i,1).^2) 2.*EI(i)./L(i,1) 0 -6.*EI(i)./(L(i,1).^2) 4.*EI(i)./L(i,1)];
    R((6.*i)-5:6.*i,1:6)=[
    cosd(thetas(i)) sind(thetas(i)) 0 0 0 0;
    -sind(thetas(i)) cosd(thetas(i)) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cosd(thetas(i)) sind(thetas(i)) 0;
    0 0 0 -sind(thetas(i)) cosd(thetas(i)) 0;
    0 0 0 0 0 1];

kg((6.*i)-5:6.*i,1:6)=transpose(R((6.*i)-5:6.*i,1:6))*kl((6.*i)-5:6.*i,1:6)*R((6.*i)-5:6.*i,1:6);

    for j=1:6
        for h=1:6
            K(LM(i,j),LM(i,h))=K(LM(i,j),LM(i,h))+kg(j+6.*(i-1),h);
    end
end
    end  
Lcounter=1;    
LLoads=zeros(6,1);
GLoads=zeros(NON.*3,1);
while 1
    if Lcounter>1
        disp('Press "F" when all loads are completed')
    else
    end
    a=['For the ',num2str(Lcounter),'. load press "N" if the Load is on Node or Press "M" if the Load is on Member']; 
    disp(a)
    prompt='\n';
    str=input(prompt,'s');
    
if str=='M' || str=='m'
    prompt='Choose the loading type. \n Press "D" for Distributed Loading \n Press "P" for Point Load \n';
    str=input(prompt,'s');
    
switch str
    case {'D','d'}

%Member üstünde external moment girmek için external momentin olduðu noktaya Node atanmalý         
        
%Distributed Loads on Member
% i=eleman numarasý
% Bazý node larda local member yükleri üstüste biniceði için her döngüde
% local loadlarý sýfýrlamak lazým
LLoads=zeros(6,1);
prompt = {'Enter the number of the member that is carrying the Load','Enter the magnitude of Distributed Load in kN/m','Enter the direction of the Load ("U" for upwards and "D" for downwards with respect to Local Coordinates of that member)'};
dlgtitle = 'Distributed Load on Member';
answer = inputdlg(prompt,dlgtitle);
dir=cell2mat(answer(3,1));
W=cell2mat(answer(2,1));
W=abs(str2double(W));
i=str2double(cell2mat(answer(1,1)));
if dir=='u'||dir=='U'
LLoads(2:3,1)=LLoads(2:3,1)-[(W.*L(i,1))./2;(W.*(L(i,1)).^2)./12];
LLoads(5:6,1)=LLoads(5:6,1)-[(W.*L(i,1))./2;-(W.*(L(i,1)).^2)./12];
LCLoads=transpose(R((6.*i)-5:6.*i,1:6))*LLoads;
else
LLoads(2:3,1)=LLoads(2:3,1)+[(W.*L(i,1))./2;(W.*(L(i,1)).^2)./12];
LLoads(5:6,1)=LLoads(5:6,1)+[(W.*L(i,1))./2;-(W.*(L(i,1)).^2)./12];
LCLoads=transpose(R((6.*i)-5:6.*i,1:6))*LLoads;
end
for k=1:6
    GLoads(LM(i,k),1)=GLoads(LM(i,k),1)+LCLoads(k,1);
end

case {'P','p'}
%Point Loads on Member
LLoads=zeros(6,1);
prompt = {'Enter the number of the member that is carrying the Load','Enter the magnitude of Point Load in kN','Enter the direction of the Load ("U" for upwards and "D" for downwards with respect to Local Coordinates of that member)','Enter the distance between the Point of Load and starting node of that member in meters'};
dlgtitle = 'Point Load On Member';
answer = inputdlg(prompt,dlgtitle);
dir=cell2mat(answer(3,1));
P=cell2mat(answer(2,1));
P=abs(str2double(P));
i=str2double(cell2mat(answer(1,1)));
a=abs(str2double(cell2mat(answer(4,1))));
b=L(i,1)-a;
d2=(P.*a+((P.*(a.^2).*b)./(L(i,1).^2))-((P.*(b.^2).*a)./(L(i,1).^2)))./(a+b);
d1=P-d2;
if dir=='u'||dir=='U'
    
LLoads(2:3,1)=LLoads(2:3,1)-[d1;((P.*(b.^2).*a)./(L(i,1).^2))];
LLoads(5:6,1)=LLoads(5:6,1)-[d2;-((P.*(a.^2).*b)./(L(i,1).^2))];
LCLoads=transpose(R((6.*i)-5:6.*i,1:6))*LLoads;
else
LLoads(2:3,1)=LLoads(2:3,1)+[d1;((P.*(b.^2).*a)./(L(i,1).^2))];
LLoads(5:6,1)=LLoads(5:6,1)+[d2;-((P.*(a.^2).*b)./(L(i,1).^2))];
LCLoads=transpose(R((6.*i)-5:6.*i,1:6))*LLoads;
end

for k=1:6
    GLoads(LM(i,k),1)=GLoads(LM(i,k),1)+LCLoads(k,1);
end

    otherwise
        disp('Error:Unexpected Key')
end
elseif str=='N'||str=='n'
    %Node Loadings
    prompt='Choose the loading type. \n Press "P" for Point Load \n Press "M" for External Moment \n';
    str=input(prompt,'s');
    %Point Load On Node
    if str=='p'||str=='P'
prompt = {'Enter the number of the Node that is carrying the Load','Enter the magnitude of Point Load in kN','Enter the direction of the Load ("U" for upwards, "D" for downwards "R" for rightwards and "L" for leftwards in Global Coordinates)'};
dlgtitle = 'Point Load on Node';
answer = inputdlg(prompt,dlgtitle);
dir=string(answer(3,1));
P=abs(str2double(cell2mat(answer(2,1))));
i=str2double(cell2mat(answer(1,1)));
            switch dir
                case {'u','U'}
                    GLoads(DOFST(i,2),1)=GLoads(DOFST(i,2),1)-P;
                case {'d','D'}
                    GLoads(DOFST(i,2),1)=GLoads(DOFST(i,2),1)+P;
                case {'L','l'}
                    GLoads(DOFST(i,1),1)=GLoads(DOFST(i,1),1)+P;
                case {'r','R'}
                    GLoads(DOFST(i,1),1)=GLoads(DOFST(i,1),1)-P;
                otherwise
                    disp('Error:Unexpected Key')
            end
      %External Moment On Node               
    elseif str=='M'||str=='m'
prompt = {'Enter the number of the Node that is carrying the External Moment','Enter the magnitude of External Moment','Enter the direction of the Moment ("CCW" for Counter-Clockwise and "CW" for Clockwise)'};
dlgtitle = 'External Moment on Node';
answer = inputdlg(prompt,dlgtitle);
dir=string(answer(3,1));
M=abs(str2double(cell2mat(answer(2,1))));
i=str2double(cell2mat(answer(1,1)));
            switch dir
                case {'ccw','Ccw','CCW'}
                    GLoads(DOFST(i,3),1)=GLoads(DOFST(i,3),1)-M;
                case {'cw','Cw','CW'}
                    GLoads(DOFST(i,3),1)=GLoads(DOFST(i,3),1)+M;
                otherwise
                    disp('Error:Unexpected Key')
            end
    end
elseif str=='F'||str=='f'
    break
else 
    disp('Error:Unexpected key.')  
end 
    Lcounter=Lcounter+1;
end

Ur=zeros(restr,1);
Settlement=0;    
disp('Are there any support settlements? (Y/N)')
prompt=('\n');
str=input(prompt,'s');
if str=='Y'||str=='y'
    disp('At which node does the settlement occur? (Enter just the number of the Node)')
    prompt=('\n');
    NSS=input(prompt);
    disp('Enter the support settlement in centimeters')
    prompt=('\n');
    Settlement=input(prompt);
    LNSS=DOFS(2,NSS)-free;
    Ur(LNSS,1)=-abs((Settlement.*0.01));
elseif str=='n'||str=='N'
else
    Disp('Error:Unexpected Key')
end
    
    Kff=K(1:free,1:free);
    Kfr=K(1:free,free+1:end);
    Krf=K(free+1:end,1:free);
    Krr=K(free+1:end,free+1:end);
    
    
    Ff=zeros(free,1);
    Fef=zeros(restr,1);
    Ff(1:free,1)=-GLoads(1:free,1);
    Fef(1:restr,1)=GLoads(free+1:end,1);
  Uf=inv(K(1:free,1:free))*(Ff-(K(1:free,free+1:end)*Ur)) 
  Support_Reactions=(K(free+1:end,1:free)*Uf+K(free+1:end,free+1:end)*Ur)+Fef
  
for i=1:NON
Node_Number(i,1:6)=[num2str(i),'.Node'];
end
  for i=1:NON
      if DOFST(i,1)<free+1
      x_direction(i,1)=Uf((DOFST(i,1)),1);
      else
      x_direction(i,1)=0;
      end
  end
    for i=1:NON
      if DOFST(i,2)<free+1
      y_direction(i,1)=Uf((DOFST(i,2)),1);
      else
      y_direction(i,1)=0;
      end
    end
    for i=1:NON
      if DOFST(i,3)<free+1
      Rotation(i,1)=Uf((DOFST(i,3)),1);
      else
      Rotation(i,1)=0;
      end
    end
    
    if Settlement~=0
      y_direction(NSS,1)=-abs((Settlement.*0.01));
    else
    end
DISPLACEMENTS=table(Node_Number,x_direction,y_direction,Rotation)
  for i=1:NON
      if DOFST(i,1)>free
      x_Direction(i,1)=Support_Reactions((DOFST(i,1)-free),1);
      else
      x_Direction(i,1)=0;
      end
  end
    for i=1:NON
      if DOFST(i,2)>free
      y_Direction(i,1)=Support_Reactions((DOFST(i,2)-free),1);
      else
      y_Direction(i,1)=0;
      end
    end
    for i=1:NON
      if DOFST(i,3)>free
      Rotational(i,1)=Support_Reactions((DOFST(i,3)-free),1);
      else
      Rotational(i,1)=0;
      end
    end
SUPPORT_REACTIONS=table(Node_Number,x_Direction,y_Direction,Rotational)
disp('Report any bug to doguhanyazar@hotmail.com        Doðuhan Yazar')

