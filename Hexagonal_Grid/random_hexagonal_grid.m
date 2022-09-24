%versione griglia esagoni con nodi sui bordi equidistanti per i due assi

clear
close all

x0=0; y0=0;
f1=1; f2=1;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

n=6;
suddx=5;
alpha = 2*pi/n; %angolo tra due vertici
theta = pi/2; %angolo di rotazione
d=0.2;
r_lato = (f1/(suddx))/2;
r_vertice = r_lato/cos(alpha/2);
suddy = ceil(f2/((3*r_vertice)/2));

x_centro=xmin+r_lato;
y_centro=ymax-r_vertice/2;

indelem=0;
elem_x=cell(suddx*suddy,1);
elem_y=cell(suddx*suddy,1);
elem=cell(suddx*suddy,1);

ultimo_nodo=((y_centro+r_vertice)-(3*r_vertice/2)*(suddy-1)-r_vertice/2);
if ultimo_nodo-ymin<r_vertice/4 && suddx>1
    suddy=suddy-1;
end

dy=f2/suddy;
%% CREAZIONE GRIGLIA ESAGONI
hold on
for i=1:suddy
    for j=1:suddx 
        indelem=indelem+1;
        %creo l'esagono, theta indica la rotazione
        for k=1:n       
           x(k)=x_centro+r_vertice*cos(alpha*(k-1)+theta);
           y(k)=y_centro+r_vertice*sin(alpha*(k-1)+theta);
        end
        if i==1
            if j==1
                y(3)=ymax-dy;
                if suddx==1
                  y(5)=y(3);
                end
            elseif j==suddx
                y(end-1)=ymax-dy;
            end
            x(1)=[];
            y(1)=[];
        elseif i==suddy  
            if rem(suddy,2)==0            
               x=[x(end), x(1), x(1), x(end)+r_lato, x(end)+r_lato];
               y=[y(end), y(1), ymin, ymin, y(1)]; 
               if j==1
                y(2)=ymax-dy*(i-1);
                if suddx==1
                  y(end)=y(2);
                end
               elseif j==suddx
                   y(end)=ymax-dy*(i-1);
               end
            else  
              y(3)=ymin;
              y(5)=ymin;
              y(4)=[];
              x(4)=[];  
              if j==1
                y(2)=ymax-dy*(i-1);
               elseif j==suddx
                   y(end)=ymax-dy*(i-1);
               end
            end
        else
            if j==1
                if rem(i,2)==0
                  x=[x(1), x(4:end)];
                  y=[ymax-dy*(i-1), ymax-dy*i, y(5:end)];
                else          
                  y=[y(1), ymax-dy*(i-1), ymax-dy*i, y(4:end)];
                end
            elseif j==suddx && rem(i,2)==1
                y(end)=ymax-dy*(i-1); 
                y(end-1)=ymax-dy*i;
            end
        end
        
        %plot([x, x(1)], [y, y(1)])
        %plot(mean(x),mean(y),'m*')
        elem_x{indelem,:}=x;
        elem_y{indelem,:}=y;
        if indelem==1
          nodi_x = x;
          nodi_y = y;
        else
          nodi_x = [nodi_x x];
          nodi_y = [nodi_y y];
        end
        x_centro=x_centro+2*r_lato;
    end
    
    if rem(i,2)==0 && i<suddy
        indelem=indelem+1;
            for k=1:n       
               x(k)=x0+(x_centro-xmin)+r_vertice*cos(alpha*(k-1)+theta);
               y(k)=y0+(y_centro-ymin)+r_vertice*sin(alpha*(k-1)+theta);
            end 
            y(1)=ymax-dy*(i-1); 
            y(4)=ymax-dy*i;
            x=x(1:4);
            y=y(1:4);    
            %plot([x, x(1)], [y, y(1)])
            %plot(mean(x),mean(y),'m*')
            elem_x{indelem,:}=x;
            elem_y{indelem,:}=y;
            nodi_x = [nodi_x x];
            nodi_y = [nodi_y y];
        
        x_centro=xmin+r_lato;
    else          
        x_centro=xmin;
    end
    y_centro=y_centro-(r_vertice+r_vertice/2);
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
pause(0.001);
if ss ~= 30
    clf
end
end