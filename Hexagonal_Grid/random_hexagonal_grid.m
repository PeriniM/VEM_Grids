%versione che crea una griglia di esagoni, numerando i nodi in senso
%antiorario, e sposta i vertici in modo random, senza eliminazione di
%triangoli sui bordi

clear
close all

x0=0; y0=0;
f1=10; f2=10;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  

n=6;
alpha = 2*pi/n; 
suddx = 3;
r_lato = f1/(2*suddx-1); %calcolato in modo che si parta con un vertice
% e si finisce con un lato -> r_lato=(f1+r_lato)/(2*suddx)
r_vertice = r_lato/cos(alpha/2);
suddy = ceil(f2/((3*r_vertice)/2))+1; %per riempire tutto il dominio

d=0.1; % fattore moltiplicativo per i vertici random

x=xmin+r_lato;
y=ymax+r_vertice+r_vertice/2;
indelem=0;
elem_x=cell(suddx*suddy,1);
elem_y=cell(suddx*suddy,1);
elem=cell(suddx*suddy,1);

for i=1:suddy
    for j=1:suddx
      indelem=indelem+1;
      [xprov, yprov] = creaesagono_v1(x,y,n,suddx,f1,xmin,xmax,ymin,ymax);
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
        x=xmin+r_lato;
    else       
        x=xmin;
    end
    y=y-(r_vertice+r_vertice/2);
end

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

xvert_rand=xvert;
yvert_rand=yvert;

for ss=1:50
% costruzioni nodi di frontiera
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

% disegno gli elementi con i loro nodi

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
griglia.vertices=[xvert_rand yvert_rand];
griglia.bordo=b(:);
axis equal
pause(0.01);
clf
end