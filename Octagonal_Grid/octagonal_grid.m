%versione griglia ottagoni con nodi sui bordi equidistanti per i due assi

clear
close all

x0=0; y0=0;
f1=1; f2=1;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

n=8;
suddx=5;

alpha = 2*pi/n; %angolo tra due vertici
theta = pi/2+alpha/2; %angolo di rotazione

if suddx==1
    suddx=2; 
end

a=(f1/(suddx-1))/2;
suddy = ceil(abs(ymax-a-ymin)/(2*a))+1;

r_vertice = a/cos(alpha/2);
r_lato = sin(alpha/2)*2*r_vertice;

x=zeros(1,n);
y=zeros(1,n);

x_centro=xmin;
y_centro=ymax;

indelem=0;
elem_x=cell(suddx*suddy+(suddx-1)*(suddy-1),1);
elem_y=cell(suddx*suddy+(suddx-1)*(suddy-1),1);
elem=cell(suddx*suddy+(suddx-1)*(suddy-1),1);

dx=f1/suddx;
dy=f2/suddy;

hold on
for i=1:suddy
    for j=1:suddx            
        indelem=indelem+1;
        %creo l'ottagono, theta indica la rotazione
        for k=1:n       
           x(k)=x_centro+r_vertice*cos(alpha*(k-1)+theta);
           y(k)=y_centro+r_vertice*sin(alpha*(k-1)+theta);
        end 
        %CREAZIONE ROMBI
        if j<suddx && i<suddy
            x_q=[x(6) x(5) x(6) x(6)+abs(x(6)-x(5))];
            y_q=[y(6) y(5) y(5)-abs(y(6)-y(5)) y(5)];
            elem_x{indelem+suddx,:}=x_q;
            elem_y{indelem+suddx,:}=y_q;
        end
        %MODIFICA OTTAGONI
        if i==1
            if j==1
                x=[xmin xmin x(5) x(6) xmin+dx];
                y=[ymax ymax-dy y(4) y(6) ymax];
            elseif j==suddx
                x=[xmax-dx x(3) x(4) xmax xmax];
                y=[ymax y(3) y(4) ymax-dy ymax];
            else
                x=[xmin+dx*(j-1) x(3:6) xmin+dx*j];
                y=[ymax y(3:6) ymax];
            end          
        elseif i==suddy
            if j==1
                x=[xmin xmin xmin+dx x(7:end)];
                y=[ymin+dy ymin ymin y(7:end)];
            elseif j==suddx
                x=[xmax x(1:2) xmax-dx xmax];
                y=[ymin+dy y(1:2) ymin ymin];
            else
                x=[x(1:2) xmin+dx*(j-1) xmin+dx*j x(7:end)];
                y=[y(1:2) ymin ymin y(7:end)];
            end 
        else
            if j==1
                x=[xmin xmin x(5:end)];
                y=[ymax-dy*(i-1) ymax-dy*i y(5:end)];
            elseif j==suddx
                x=[xmax x(1:4) xmax];
                y=[ymax-dy*(i-1) y(1:4) ymax-dy*i];
            end
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
        if j==suddx && i<suddy
            indelem=indelem+suddx-1;
        end
        
        x_centro=x_centro+2*a;
    end
    x_centro=xmin;
    y_centro=y_centro-2*a;
end
%% NUMERAZIONE NODI
linee_x=uniquetol(nodi_x); %usa quicksort + ordine crescente
linee_y=fliplr(uniquetol(nodi_y)); %in ordine decrescente
nodi_unici=zeros(1,2); %indici dei nodi sulle linee x e y univoche
count_nodi_unici=zeros(1,1); %numerazione dei nodi per ciascun vertice
count_nodi_globali=0; %contatore per l'incremento dei nodi globali
xvert=zeros(1,1);
yvert=zeros(1,1);

for s=1:indelem
    for k=1:length(elem_x{s,:}) 
        count_nodi_globali=count_nodi_globali+1;
        %trova la posizione del nodo nelle linee univoche con tolleranza
        ind_pos_x=ismembertol(linee_x,elem_x{s,:}(1,k));
        ind_pos_y=ismembertol(linee_y,elem_y{s,:}(1,k));
        if s==1
            count_nodi_unici(end+1)=count_nodi_globali;       
            elem{s,:}(1,k)=count_nodi_unici(end);
            xvert(end+1)= linee_x(ind_pos_x);
            yvert(end+1)= linee_y(ind_pos_y);
        else
            %se trova una ripetizione degli indici dei nodi
            if find(ismember(nodi_unici,[find(ind_pos_x==1) find(ind_pos_y==1)],'rows'))
                %trovo l'indice del nodo "originale" non ripetuto
                index = find(ismember(nodi_unici,[find(ind_pos_x==1) find(ind_pos_y==1)],'rows'));
                %inserisco la numerazione che avevo dato al nodo già esistente
                count_nodi_unici(end+1)=count_nodi_unici(index(1));
                %aggiorno la cella degli elementi
                elem{s,:}(1,k)=count_nodi_unici(end);
                %decremento il contatore globale
                count_nodi_globali=count_nodi_globali-1;
            else
                count_nodi_unici(end+1)=count_nodi_globali;  
                elem{s,:}(1,k)=count_nodi_unici(end);
                % metto solo le coordinate che non si ripetono
                xvert(end+1)= linee_x(ind_pos_x);
                yvert(end+1)= linee_y(ind_pos_y);
            end
        end
        nodi_unici(end+1,:)=[find(ind_pos_x==1) find(ind_pos_y==1)]; 
    end
end
%rimuovo il primo valore che si era creato quando li ho inizializzati
xvert(1)=[];
yvert(1)=[];

%% CREAZIONE PLOT
% costruzioni nodi di frontiera
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

% disegno gli elementi con i loro nodi

for iel=1:indelem
    xvertici=elem{iel,:};
    xv=xvert(xvertici);
    yv=yvert(xvertici);
    plot([xv, xv(1)],[yv, yv(1)],'k','linewidth',1)
    hold on
  %  h=text(mean(xv), mean(yv), {num2str(iel)});
   % set(h,'color','r')
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
