%versione che crea una griglia di esagoni, numerando i nodi in senso
%antiorario, e sposta i vertici in modo random, con eliminazione di
%triangoli sui bordi

clear
close all

x0=0; y0=0;
f1=10; f2=10;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

n=6;
alpha = 2*pi/n; 
suddx = 4;
r_lato = f1/(2*suddx-1); %calcolato in modo che si parta con un vertice
% e si finisce con un lato -> r_lato=(f1+r_lato)/(2*suddx)
r_vertice = r_lato/cos(alpha/2);
suddy = ceil(f2/((3*r_vertice)/2)); %per riempire tutto il dominio

d=0.1; % fattore moltiplicativo per i vertici random

x=xmin;
y=ymax;
indelem=0;
elem_x=cell(suddx*suddy,1);
elem_y=cell(suddx*suddy,1);
elem=cell(suddx*suddy,1);

elimina_riga=0;

for i=1:suddy
    for j=1:suddx
      indelem=indelem+1;
      [xprov, yprov] = createHexagon(x,y,n,suddx,f1,xmin,xmax,ymin,ymax);
      if i==1
        if length(yprov)==4
          yprov(end)= ymax; 
        elseif length(yprov)==6
          yprov(2)= ymax;
          yprov(end)= ymax;
          xprov(1)=[];
          yprov(1)=[];
        end
      end
      if i==suddy
        if length(yprov)==3 || abs(yprov(1)-ymin)<(r_vertice*0.7)
            elimina_riga=1;
            switch length(elem_x{indelem-suddx,:})
                case 4
                    if elem_x{indelem-suddx,:}(1,1)<(xmin+xmax)/2
                        elem_y{indelem-suddx,:}(1,2)=ymin;
                        elem_y{indelem-suddx,:}(1,3)=ymin;
                    else
                        elem_y{indelem-suddx,:}(1,3)=ymin;
                        elem_y{indelem-suddx,:}(1,4)=ymin;
                    end
                case 5
                    if elem_x{indelem-suddx,:}(1,1)<(xmin+xmax)/2
                        temp_x = elem_x{indelem-suddx,:}(1,:);
                        temp_y = elem_y{indelem-suddx,:}(1,:);
                        temp_x(3)=[];
                        temp_y(3)=[];
                        temp_y(end-1)=ymin;
                        elem_x{indelem-suddx,:}=temp_x;
                        elem_y{indelem-suddx,:}=temp_y;
                    else
                        temp_x = elem_x{indelem-suddx,:}(1,:);
                        temp_y = elem_y{indelem-suddx,:}(1,:);
                        temp_x(4)=[];
                        temp_y(4)=[];
                        temp_y(3)=ymin;    
                        elem_x{indelem-suddx,:}=temp_x;
                        elem_y{indelem-suddx,:}=temp_y;
                    end   
                case 6
                    temp_x = elem_x{indelem-suddx,:}(1,:);
                    temp_y = elem_y{indelem-suddx,:}(1,:);
                    temp_y(3)=ymin;
                    temp_y(5)=ymin;
                    temp_x(4)=[];
                    temp_y(4)=[];        
                    elem_x{indelem-suddx,:}=temp_x;
                    elem_y{indelem-suddx,:}=temp_y;
                case 7
                    temp_x = elem_x{indelem-suddx,:}(1,:);
                    temp_y = elem_y{indelem-suddx,:}(1,:);
                    temp_y(3)=ymin;
                    temp_y(6)=ymin;
                    temp_x(4:5)=[];
                    temp_y(4:5)=[];     
                    elem_x{indelem-suddx,:}=temp_x;
                    elem_y{indelem-suddx,:}=temp_y;
            end
        end
      end
      elem_x{indelem,:}=xprov;
      elem_y{indelem,:}=yprov;
      if indelem==1
      nodi_x = xprov;
      nodi_y = yprov;
      else
          nodi_x = [nodi_x xprov];
          nodi_y = [nodi_y yprov];
      end
      x=x+2*r_lato;
    end
    if rem(i,2)==0
        x=xmin;
    else       
        x=xmin+r_lato;
    end
    y=y-(r_vertice+r_vertice/2);
end

if elimina_riga==1
    elem_x=elem_x(1:indelem-suddx);
    elem_y=elem_y(1:indelem-suddx);
    indelem=indelem-suddx;
end

%% NODES ENUMERATION
linee_x=uniquetol(nodi_x); %quicksort + ascending order
linee_y=fliplr(uniquetol(nodi_y)); %descending order
nodi_unici=zeros(1,2); %indices of nodes laying on unique x and y lines
count_nodi_unici=zeros(1,1); %nodes numeration for each vertex
count_nodi_globali=0; %counter for global nodes
xvert=zeros(1,1);
yvert=zeros(1,1);

for s=1:indelem
    for k=1:length(elem_x{s,:}) 
        count_nodi_globali=count_nodi_globali+1;
        %find node's position on unique lines with tolerance
        ind_pos_x=ismembertol(linee_x,elem_x{s,:}(1,k));
        ind_pos_y=ismembertol(linee_y,elem_y{s,:}(1,k));
        if s==1
            count_nodi_unici(end+1)=count_nodi_globali;       
            elem{s,:}(1,k)=count_nodi_unici(end);
            xvert(end+1)= linee_x(ind_pos_x);
            yvert(end+1)= linee_y(ind_pos_y);
        else
            %if finds duplicates of nodes' indices
            if find(ismember(nodi_unici,[find(ind_pos_x==1) find(ind_pos_y==1)],'rows'))
                %find the index of the non repeated "original" node
                index = find(ismember(nodi_unici,[find(ind_pos_x==1) find(ind_pos_y==1)],'rows'));
                %insert the nueration already given to the existing node
                count_nodi_unici(end+1)=count_nodi_unici(index(1));
                %update the elements cell
                elem{s,:}(1,k)=count_nodi_unici(end);
                % decrease global counter
                count_nodi_globali=count_nodi_globali-1;
            else
                count_nodi_unici(end+1)=count_nodi_globali;  
                elem{s,:}(1,k)=count_nodi_unici(end);
                %insert only non repeated coordinates
                xvert(end+1)= linee_x(ind_pos_x);
                yvert(end+1)= linee_y(ind_pos_y);
            end
        end
        nodi_unici(end+1,:)=[find(ind_pos_x==1) find(ind_pos_y==1)]; 
    end
end
%remove first element for initialization
xvert(1)=[];
yvert(1)=[];
%% GRID PLOT
xvert_rand=xvert;
yvert_rand=yvert;

for ss=1:30
% building boundary nodes
nnode=length(xvert);
j=0;
b=zeros(1,1);
for i=1:nnode
    if abs(xvert(i)-xmin)<=1e-10 || abs(xvert(i)-xmax)<=1e-10 || abs(yvert(i)-ymin)<=1e-10 || abs(yvert(i)-ymax)<=1e-10
        j=j+1;
        b(j)=i;
    else
        xvert_rand(i)=xvert(i)+d*2*(rand-0.5)*r_lato;  
        yvert_rand(i)=yvert(i)+d*2*(rand-0.5)*r_lato;
    end
end
griglia.dirichlet=b(:);
griglia.neuman=0;

% draw elements with their nodes

for iel=1:indelem
    xvertici=elem{iel,:};
    xv=xvert_rand(xvertici);
    yv=yvert_rand(xvertici);
    plot([xv, xv(1)],[yv, yv(1)],'k','linewidth',1)
    hold on
    h=text(mean(xv), mean(yv), {num2str(iel)});
    set(h,'color','r')
end
for i=1:length(xvert)
    plot( xvert_rand(i),yvert_rand(i),'o'); text(xvert_rand(i)+0.03,yvert_rand(i)+0.03, num2str(i));   
end
for i=1:j
   plot(xvert(b(i)),yvert(b(i)),'m*')
end
griglia.elements=indelem;
griglia.vertices=[xvert_rand; yvert_rand];
griglia.bordo=b(:);
axis equal
pause(0.01);
clf
end