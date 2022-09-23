%versione che crea esagono valutando solo il nodo successivo
%senza eliminazione di triangoli
function [xprov, yprov] = createHexagon(x1,y1,n,suddx,f1,xmin,xmax,ymin,ymax)

%---------righe da rimettere se si toglie la function---------
%{
x0=0; y0=0;
f1=5; f2=5;  
xmin=x0; xmax=x0+f1;  
ymin=y0; ymax=y0+f2;  
n=6;
suddx = 5;
%}
alpha = 2*pi/n; 
r_lato = f1/(2*suddx-1); %calcolato in modo che si parta con un vertice
% e si finisce con un lato -> r_lato=(f1+r_lato)/(2*suddx)
r_vertice = r_lato/cos(alpha/2);
%suddy = ceil(f2/((3*r_vertice)/2)); %per riempire tutto il dominio

x=zeros(1,n);
y=zeros(1,n);

x(1)=x1; 
y(1)=y1;

if abs(y(1)-ymin)<1e-10
      y(1)=ymin;
end  
if abs(x(1)-xmax)<1e-10
      x(1)=xmax;
end
if abs(x(1)-xmin)<1e-10
      x(1)=xmin;
end
if abs(y(1)-ymax)<1e-10
      y(1)=ymax;
end
% metto 1 se il nodo è dentro al dominio, 0 se è fuori
nodo_check=zeros(1,n);
ind=0; %tiene il conto di quanti vertici sono dentro al dominio

%verifico che il nodo uno sia dentro il dominio
if x(1)>=xmin && x(1)<=xmax && y(1)>=ymin && y(1)<=ymax
    nodo_check(1)=1;
    ind=ind+1;
else
    nodo_check(1)=0;
end
% trovo le coordinate originali dei vertici del poligono
for i=2:n       
   alpha = pi/n-2*pi/n*(i-1);
   x(i)=x(i-1)-r_vertice*cos(alpha);
   y(i)=y(i-1)+r_vertice*sin(alpha);
   %correggo le imprecisioni dei vertici sui bordi (pigreco approssimato)
   if abs(x(i)-xmin)<1e-10
        x(i)=xmin;
   end
   if abs(y(i)-ymin)<1e-10
        y(i)=ymin;
   end  
   if abs(x(i)-xmax)<1e-10
        x(i)=xmax;
   end
   if abs(y(i)-ymax)<1e-10
        y(i)=ymax;
   end
   %segno la posizione dei vertici dentro al dominio
   if x(i)>=xmin && x(i)<=xmax && y(i)>=ymin && y(i)<=ymax
        nodo_check(i)=1;
        ind=ind+1;
   else
        nodo_check(i)=0;
   end
end

xprov=zeros(1,1); % allocamento "brutale"
yprov=zeros(1,1);
controllo=0; %se ho controllato già la x non controllo anche la y per la
             %costruzione della retta
xvert_aux=Inf;
yvert_aux=Inf;

if ind~=n %se i nodi non sono tutti dentro
    for j=1:n %scansiono ciascun nodo
        %inizio la scansione dall'ultimo vertice e vado avanti
        if j==1
           indice_nodo_attuale=n;
           indice_nodo_successivo=1;
        elseif j==2
           indice_nodo_attuale=1;
           indice_nodo_successivo=2;
        else
           indice_nodo_attuale=indice_nodo_attuale+1;
           indice_nodo_successivo=indice_nodo_attuale+1;
        end
        
        if nodo_check(indice_nodo_attuale)==0 && nodo_check(indice_nodo_successivo)==0
            % non faccio nulla
        elseif nodo_check(indice_nodo_attuale)==0 && nodo_check(indice_nodo_successivo)==1
            % verifico in che parte del dominio è fuori
            if x(indice_nodo_attuale)<xmin
                x_retta=xmin;
                y_retta=y(indice_nodo_attuale)+(y(indice_nodo_successivo)-y(indice_nodo_attuale))...
                    *(xmin-x(indice_nodo_attuale))/(x(indice_nodo_successivo)-x(indice_nodo_attuale));
                %verifico se il punto della retta interseca il dominio
                if y_retta>=ymin && y_retta<=ymax
                      if x_retta~=x(indice_nodo_successivo) || y_retta~=y(indice_nodo_successivo)
                         if indice_nodo_attuale==n
                           xprov(end+1)=x(indice_nodo_successivo);
                           yprov(end+1)=y(indice_nodo_successivo);
                           xvert_aux=x_retta;
                           yvert_aux=y_retta;
                         else
                           xprov=[xprov(1:end) x_retta x(indice_nodo_successivo)];
                           yprov=[yprov(1:end) y_retta y(indice_nodo_successivo)];
                         end 
                      else
                        xprov(end+1)=x(indice_nodo_successivo);
                        yprov(end+1)=y(indice_nodo_successivo);
                      end     
                    controllo=1;
                else
                    controllo=0;
                end
            elseif x(indice_nodo_attuale)>xmax
                x_retta=xmax;
                y_retta=y(indice_nodo_attuale)+(y(indice_nodo_successivo)-y(indice_nodo_attuale))...
                    *(xmax-x(indice_nodo_attuale))/(x(indice_nodo_successivo)-x(indice_nodo_attuale));
                %verifico se il punto della retta interseca il dominio
                if y_retta>=ymin && y_retta<=ymax
                      if x_retta~=x(indice_nodo_successivo) || y_retta~=y(indice_nodo_successivo)
                         if indice_nodo_attuale==n
                           xprov(end+1)=x(indice_nodo_successivo);
                           yprov(end+1)=y(indice_nodo_successivo);
                           xvert_aux=x_retta;
                           yvert_aux=y_retta;
                         else
                           xprov=[xprov(1:end) x_retta x(indice_nodo_successivo)];
                           yprov=[yprov(1:end) y_retta y(indice_nodo_successivo)];
                         end 
                      else
                        xprov(end+1)=x(indice_nodo_successivo);
                        yprov(end+1)=y(indice_nodo_successivo);
                      end    
                    controllo=1;
                else
                    controllo=0;
                end
            end
            if y(indice_nodo_attuale)<ymin && controllo == 0
                y_retta=ymin;
                x_retta=x(indice_nodo_attuale)+(x(indice_nodo_successivo)-x(indice_nodo_attuale))...
                    *(ymin-y(indice_nodo_attuale))/(y(indice_nodo_successivo)-y(indice_nodo_attuale));
                  if x_retta~=x(indice_nodo_successivo) || y_retta~=y(indice_nodo_successivo)
                         if indice_nodo_attuale==n
                           xprov(end+1)=x(indice_nodo_successivo);
                           yprov(end+1)=y(indice_nodo_successivo);
                           xvert_aux=x_retta;
                           yvert_aux=y_retta;
                         else
                           xprov=[xprov(1:end) x_retta x(indice_nodo_successivo)];
                           yprov=[yprov(1:end) y_retta y(indice_nodo_successivo)];
                         end 
                   else
                        xprov(end+1)=x(indice_nodo_successivo);
                        yprov(end+1)=y(indice_nodo_successivo);
                   end
            elseif y(indice_nodo_attuale)>ymax && controllo == 0
                y_retta=ymax;
                x_retta=x(indice_nodo_attuale)+(x(indice_nodo_successivo)-x(indice_nodo_attuale))...
                    *(ymax-y(indice_nodo_attuale))/(y(indice_nodo_successivo)-y(indice_nodo_attuale));
                  if x_retta~=x(indice_nodo_successivo) || y_retta~=y(indice_nodo_successivo)
                         if indice_nodo_attuale==n
                           xprov(end+1)=x(indice_nodo_successivo);
                           yprov(end+1)=y(indice_nodo_successivo);
                           xvert_aux=x_retta;
                           yvert_aux=y_retta;
                         else
                           xprov=[xprov(1:end) x_retta x(indice_nodo_successivo)];
                           yprov=[yprov(1:end) y_retta y(indice_nodo_successivo)];
                         end 
                   else
                        xprov(end+1)=x(indice_nodo_successivo);
                        yprov(end+1)=y(indice_nodo_successivo);
                   end   
            end     
        elseif nodo_check(indice_nodo_attuale)==1 && nodo_check(indice_nodo_successivo)==0 
            % verifico in che parte del dominio è fuori
            if x(indice_nodo_successivo)<xmin
                x_retta=xmin;
                y_retta=y(indice_nodo_attuale)+(y(indice_nodo_successivo)-y(indice_nodo_attuale))...
                    *(xmin-x(indice_nodo_attuale))/(x(indice_nodo_successivo)-x(indice_nodo_attuale));
                %verifico se il punto della retta interseca il dominio
                if y_retta>=ymin && y_retta<=ymax
                      if x_retta~=x(indice_nodo_attuale) || y_retta~=y(indice_nodo_attuale)
                       xprov(end+1)=x_retta;
                       yprov(end+1)=y_retta;
                      else
                        
                      end     
                    controllo=1;
                else
                    controllo=0;
                end
            elseif x(indice_nodo_successivo)>xmax
                x_retta=xmax;
                y_retta=y(indice_nodo_attuale)+(y(indice_nodo_successivo)-y(indice_nodo_attuale))...
                    *(xmax-x(indice_nodo_attuale))/(x(indice_nodo_successivo)-x(indice_nodo_attuale));
                %verifico se il punto della retta interseca il dominio
                if y_retta>=ymin && y_retta<=ymax
                      if x_retta~=x(indice_nodo_attuale) || y_retta~=y(indice_nodo_attuale)
                       xprov(end+1)=x_retta;
                       yprov(end+1)=y_retta;
                      else
                        
                      end     
                    controllo=1;
                else
                    controllo=0;
                end
            end
            if y(indice_nodo_successivo)<ymin && controllo == 0
                y_retta=ymin;
                x_retta=x(indice_nodo_attuale)+(x(indice_nodo_successivo)-x(indice_nodo_attuale))...
                    *(ymin-y(indice_nodo_attuale))/(y(indice_nodo_successivo)-y(indice_nodo_attuale));
                  if x_retta~=x(indice_nodo_attuale) || y_retta~=y(indice_nodo_attuale)
                     xprov(end+1)=x_retta;
                     yprov(end+1)=y_retta;
                  else
                     
                  end   
            elseif y(indice_nodo_successivo)>ymax && controllo == 0
                y_retta=ymax;
                x_retta=x(indice_nodo_attuale)+(x(indice_nodo_successivo)-x(indice_nodo_attuale))...
                    *(ymax-y(indice_nodo_attuale))/(y(indice_nodo_successivo)-y(indice_nodo_attuale));
                  if x_retta~=x(indice_nodo_attuale) || y_retta~=y(indice_nodo_attuale)
                     xprov(end+1)=x_retta;
                     yprov(end+1)=y_retta;
                  else
                    
                  end   
            end  
        elseif nodo_check(indice_nodo_attuale)==1 && nodo_check(indice_nodo_successivo)==1
            xprov(end+1)=x(indice_nodo_successivo);
            yprov(end+1)=y(indice_nodo_successivo);
        end
        controllo=0;
    end
    xprov(1)=[];
    yprov(1)=[];
else
    % il poligono ha tutti i nodi dentro
    xprov=x; yprov=y;
end
% per caso 0-1 ultimo-primo vertice del poligono
if xvert_aux~=Inf
    xprov(end+1)=xvert_aux;
    yprov(end+1)=yvert_aux;
end
j=0;
b=zeros(1,1); %trovo i nodi sul bordo
for i=1:length(xprov)
    if abs(xprov(i)-xmin)<=1e-10 || abs(xprov(i)-xmax)<=1e-10 || abs(yprov(i)-ymin)<=1e-10 || abs(yprov(i)-ymax)<=1e-10
        j=j+1;
        b(j)=i;
    end
end
x_centro=mean(x);
y_centro=mean(y);
cord_spigoli = [xmin, ymin;
                xmax, ymin;
                xmax, ymax;
                xmin, ymax];
%trovo la distanza euclidea dal centro del poligono agli spigoli del dominio
V= vecnorm((cord_spigoli-[x_centro,y_centro])');
[dist_spigolo,ind_spigolo]=min(V);
% se ci sono due nodi sul bordo e hanno tutte le coordinate diverse
if j==2 && ind~=n && dist_spigolo<r_vertice 
     if ind_spigolo==1
          xprov=[xprov(1:b(1)) cord_spigoli(ind_spigolo,1) xprov(b(2):end)];
          yprov=[yprov(1:b(1)) cord_spigoli(ind_spigolo,2) yprov(b(2):end)];
      elseif ind_spigolo==2
          xprov=[xprov(1:b(2)) cord_spigoli(ind_spigolo,1) xprov(b(2)+1:end)];
          yprov=[yprov(1:b(2)) cord_spigoli(ind_spigolo,2) yprov(b(2)+1:end)];
      elseif ind_spigolo==3
          xprov=[xprov(1:b(2)) cord_spigoli(ind_spigolo,1) xprov(b(2)+1:end)];
          yprov=[yprov(1:b(2)) cord_spigoli(ind_spigolo,2) yprov(b(2)+1:end)];
      elseif ind_spigolo==4
          xprov=[xprov(1:b(1)) cord_spigoli(ind_spigolo,1) xprov(b(2):end)];
          yprov=[yprov(1:b(1)) cord_spigoli(ind_spigolo,2) yprov(b(2):end)];
     end
end
hold on
%plot([xprov, xprov(1)],[yprov, yprov(1)],'k','linewidth',1)
plot([xmin, xmax, xmax, xmin, xmin],[ymin, ymin, ymax, ymax, ymin])
axis equal

end