clear
close all

x0=0; y0=0;
f1=1; f2=1;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

suddx=5;
suddy=5;
dx=f1/suddx;
dy=f2/suddy;
if suddx==1
    x_l = [0 0 dx dx];
    y_l = [dy 0 0 dy];
else
    x_l = [0 0 dx dx+dx/2 dx];
    y_l = [dy 0 0 dy/2 dy];
end
x_m = [0 dx/2 0 dx dx+dx/2 dx];
y_m = [dy dy/2 0 0 dy/2 dy];
x_r = [0 dx/2 0 dx dx];
y_r = [dy dy/2 0 0 dy];

x_start = xmin;
y_start = ymax-dy;

indelem=0;
elem_x=cell(suddx*suddy,1);
elem_y=cell(suddx*suddy,1);
elem=cell(suddx*suddy,1);

for i=1:suddy
    for j=1:suddx
        indelem=indelem+1;
        if j==1
            x=x_l+x_start;
            y=y_l+y_start;
        elseif j==suddx
            x=x_r+x_start;
            y=y_r+y_start;
        else
            x=x_m+x_start;
            y=y_m+y_start;
        end
        
        elem_x{indelem,:}=x;
        elem_y{indelem,:}=y;
        if indelem==1
          nodi_x = x;
          nodi_y = y;
        else
          nodi_x = [nodi_x x];
          nodi_y = [nodi_y y];
        end
        
        x_start=x_start+dx;
        %plot([x x(1)], [y y(1)])
    end
    y_start=y_start-dy;
    x_start=xmin;
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
%building boundary nodes
nnode=length(xvert);
j=0;
b=zeros(1,1);
for i=1:nnode
    if abs(xvert(i)-xmin)<=1e-10 || abs(xvert(i)-xmax)<=1e-10 || abs(yvert(i)-ymin)<=1e-10 || abs(yvert(i)-ymax)<=1e-10
        j=j+1;
        b(j)=i;
    end
end
griglia.dirichlet=b(:);
griglia.neuman=0;

%draw elements with their nodes

for iel=1:indelem
    xvertici=elem{iel,:};
    xv=xvert(xvertici);
    yv=yvert(xvertici);
    plot([xv, xv(1)],[yv, yv(1)],'k','linewidth',1)
    hold on
    %h=text(mean(xv), mean(yv), {num2str(iel)});
    %set(h,'color','r')
end

%for i=1:length(xvert)
%    plot( xvert(i),yvert(i),'o'); text(xvert(i)+0.03,yvert(i)+0.03, num2str(i));   
%end

%for i=1:j
%   plot(xvert(b(i)),yvert(b(i)),'m*')
%end

griglia.elements=indelem;
griglia.vertices=[xvert;yvert];
griglia.bordo=b(:);
axis equal

