%versione griglia pentagoni con nodi sui bordi equidistanti per i due assi

clear
close all

x0=0; y0=0;
f1=1; f2=1;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

n=5;
suddx=6;
if suddx==1
    suddx=2;
end
alpha = 2*pi/n; %angolo tra due vertici
theta = pi/2; %angolo di rotazione

d = f1/(suddx-1); %diagonale pentagono, dimensione più grande 
r_vertice = d/(2*sin(alpha));
a = r_vertice*cos(alpha/2);
suddy = ceil(f2/(a+r_vertice));

x=zeros(1,n);
y=zeros(1,n);

x_centro=xmin;
y_centro=ymax-(r_vertice-tan(alpha/2)*d/2);

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
        %creo il pentagono, theta indica la rotazione
        for k=1:n       
           x(k)=x_centro+r_vertice*cos(alpha*(k-1)+theta);
           y(k)=y_centro+r_vertice*sin(alpha*(k-1)+theta);
        end
        %riordino la posizione dei vertici del pentagono capovolto
        if rem(i,2)==0
            x=[x(4:end) x(1:3)];
            y=[y(4:end) y(1:3)];
        end  
        %CREAZIONE E MODIFICA QUADRATI
        if rem(i,2)==0 %se siamo in una riga pari
            if j<suddx && i<suddy %escludo l'ultimo quadrato a destra fuori dal dominio
                if y(4)>ymin %se il quadrato non è completamente fuori                    
                    if j==1
                        x_q=[x(4), x(3), x(4), x(4)+abs(x(3)-x(4))];
                        y_q=[y(4), ymax-dy*i, y(3)-abs(y(4)-y(3)), y(3)];  
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    elseif j==suddx-1
                        x_q=[x(4), x(3), x(4), x(4)+abs(x(3)-x(4))];
                        y_q=[y(4), y(3), y(3)-abs(y(4)-y(3)), ymax-dy*i];  
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    else
                        x_q=[x(4), x(3), x(4), x(4)+abs(x(3)-x(4))];
                        y_q=[y(4), y(3), y(3)-abs(y(4)-y(3)), y(3)];  
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    end
                end
            end
        else
            if j<suddx && i<suddy
                if y(end)>ymin
                    if i==suddy-1  
                        x_q=[x(end), x(end-1), xmin+dx*j, x(end)+abs(x(end)-x(end-1))];
                        y_q=[y(end), y(end-1), ymin, y(end-1)];
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    elseif i==1
                        x_q=[xmin+dx*j, x(end-1), x(end), x(end)+abs(x(end)-x(end-1))];
                        y_q=[y(end), y(end-1), y(end-1)-abs(y(end)-y(end-1)), y(end-1)];
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    else
                        x_q=[x(end), x(end-1), x(end), x(end)+abs(x(end)-x(end-1))];
                        y_q=[y(end), y(end-1), y(end-1)-abs(y(end)-y(end-1)), y(end-1)];
                        elem_x{indelem+suddx,:}=x_q;
                        elem_y{indelem+suddx,:}=y_q;
                    end
                end
            end
        end
        
        % MODIFICA PENTAGONI
        if rem(i,2)==0
            if i>1 && i<suddy
                if j==1
                    x=[xmin xmin x(4:5)];
                    y=[ymax-dy*(i-1) ymax-dy*i y(4:5)];
                elseif j==suddx
                    x=[x(1:2) xmax xmax];
                    y=[y(1:2) ymax-dy*i ymax-dy*(i-1)];
                end
            elseif i==suddy
                if j==1
                    x=[xmin xmin xmin+dx x(end)];
                    y=[ymin+dy ymin ymin y(end)];
                elseif j==suddx
                    x=[x(1) xmax-dx xmax xmax];
                    y=[y(1) ymin ymin ymin+dy];
                else
                    y(2)=ymin;
                    y(4)=ymin;
                    y(3)=[];
                    x(3)=[];
                    x(2)=xmin+dx*(j-1);
                    x(3)=xmin+dx*j;
                end
            end
        else
            if i==1
                if j==1
                    x=[xmin xmin x(4), xmin+dx];
                    y=[y(2) ymax-dy y(4:end)];
                elseif j==suddx
                    x=[xmax xmax-dx x(3) xmax];
                    y=[ymax y(2) y(3) ymax-dy];
                else   
                    y(1)=[];
                    x=[xmin+dx*(j-1) x(3:end-1) xmin+dx*j];              
                end
            elseif i==suddy
                if j==1
                    x=[x(1) xmin xmin+dx x(end)];
                    y=[ymin+dy ymin ymin y(end)];
                elseif j==suddx
                    x=[xmax x(2) xmax-dx xmax];
                    y=[ymin+dy y(2) ymin ymin];
                else
                    x=[x(1) x(2) xmin+dx*(j-1) xmin+dx*j x(end)];
                    y=[y(1) y(2) ymin ymin y(end)];
                end
            else
                if j==1
                    x=[x(1) xmin x(4) x(end)];
                    y=[ymax-dy*(i-1) ymax-dy*i y(4) y(5)];
                elseif j==suddx
                    x=[x(1:3) xmax];
                    y=[ymax-dy*(i-1) y(2:3) ymax-dy*i];
                end
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
        x_centro=x_centro+d;
    end
    if rem(i,2)==0
        theta=pi/2;
        y_centro=y_centro-2*r_vertice;
    else       
        theta=3*pi/2;
        y_centro=y_centro-2*a;
    end
    x_centro=xmin; 
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
   % h=text(mean(xv), mean(yv), {num2str(iel)});
    %set(h,'color','r')
end
%for i=1:length(xvert)
%    plot( xvert(i),yvert(i),'o'); text(xvert(i)+0.03,yvert(i)+0.03, num2str(i));   
%end

%for i=1:j
%   plot(xvert(b(i)),yvert(b(i)),'m*')
%end

griglia.elements=indelem;
griglia.vertices=[xvert; yvert];
griglia.bordo=b(:);
axis equal

