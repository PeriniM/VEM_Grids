%versione che modifica il pentagono eliminando i nodi fuori dal dominio
function [xprov, yprov] = createHexagon(x1,y1,n,suddx,f1,xmin,xmax,ymin,ymax)
%{
%---------righe da rimettere se si toglie la function---------
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

ind_pos = find(nodo_check==1); % mi da l'indice dei vertici dentro al dominio
                               % se tutti i vertici sono dentro al dominio
                               % allora si avrà il poligono

xprov=zeros(1,1); % allocamento "brutale"
yprov=zeros(1,1);
controllo=0; %se ho controllato già la x non controllo anche la y per la
             %costruzione della retta

             
if ind~=n %se i nodi non sono tutti dentro
    for j=1:ind %scansiono ciascun nodo interno al dominio
        if ind_pos(j)==1
            indice_nodo_precedente=n;
            indice_nodo_successivo=2;
        elseif ind_pos(j)==n
            indice_nodo_successivo=1;
            indice_nodo_precedente=n-1;
        else
            indice_nodo_precedente=ind_pos(j)-1;
            indice_nodo_successivo=ind_pos(j)+1;
        end
        
        if nodo_check(indice_nodo_precedente)==0 %se il nodo precedente è fuori
            % verifico in che parte del dominio è fuori
            if x(indice_nodo_precedente)<xmin
                xprov(end+1)=xmin;
                yprov(end+1)=y(ind_pos(j))+(y(indice_nodo_precedente)-y(ind_pos(j)))...
                    *(xmin-x(ind_pos(j)))/(x(indice_nodo_precedente)-x(ind_pos(j)));
                %verifico se il punto della retta interseca il dominio
                if yprov(end)>=ymin && yprov(end)<=ymax        
                      % se il nuovo nodo trovato è uguale al nodo dentro allora non
                      % aggiungo il nodo due volte, altrimenti aggiungo prima il
                      % nodo nuovo e poi quello successivo valido
                      if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                         xprov(end+1)=x(ind_pos(j));
                         yprov(end+1)=y(ind_pos(j));
                         if ind==1
                            xprov1=[x(ind_pos(j)) xprov(end-1)];
                            yprov1=[y(ind_pos(j)) yprov(end-1)];
                         end
                      else
                          if ind==1
                            xprov1=[x(ind_pos(j))];
                            yprov1=[y(ind_pos(j))];
                         end
                      end    
                    controllo=1;
                else
                     xprov(end)=[];
                     yprov(end)=[];
                     controllo=0;
                end      
            elseif x(indice_nodo_precedente)>xmax
                xprov(end+1)=xmax;
                yprov(end+1)=y(ind_pos(j))+(y(indice_nodo_precedente)-y(ind_pos(j)))...
                    *(xmax-x(ind_pos(j)))/(x(indice_nodo_precedente)-x(ind_pos(j)));
                %verifico se il punto della retta interseca il dominio
                if yprov(end)>=ymin && yprov(end)<=ymax
                    % se il nuovo nodo trovato è uguale al nodo dentro allora non
                      % aggiungo il nodo due volte, altrimenti aggiungo prima il
                      % nodo nuovo e poi quello successivo valido
                      if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                         xprov(end+1)=x(ind_pos(j));
                         yprov(end+1)=y(ind_pos(j));
                         if ind==1
                            xprov1=[x(ind_pos(j)) xprov(end-1)];
                            yprov1=[y(ind_pos(j)) yprov(end-1)];
                         end
                      else
                          if ind==1
                            xprov1=[x(ind_pos(j))];
                            yprov1=[y(ind_pos(j))];
                         end
                      end  
                    controllo=1;
                else
                     xprov(end)=[];
                     yprov(end)=[];
                     controllo=0;
                end
            end
            if y(indice_nodo_precedente)<ymin && controllo == 0
                xprov(end+1)=x(ind_pos(j))+(x(indice_nodo_precedente)-x(ind_pos(j)))...
                    *(ymin-y(ind_pos(j)))/(y(indice_nodo_precedente)-y(ind_pos(j)));
                yprov(end+1)=ymin;
                
                  % se il nuovo nodo trovato è uguale al nodo dentro allora non
                  % aggiungo il nodo due volte, altrimenti aggiungo prima il
                  % nodo nuovo e poi quello successivo valido
                  if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                     xprov(end+1)=x(ind_pos(j));
                     yprov(end+1)=y(ind_pos(j));
                     if ind==1
                            xprov1=[x(ind_pos(j)) xprov(end-1)];
                            yprov1=[y(ind_pos(j)) yprov(end-1)];
                     end
                   else
                     if ind==1
                        xprov1=[x(ind_pos(j))];
                        yprov1=[y(ind_pos(j))];
                     end
                  end
            elseif y(indice_nodo_precedente)>ymax && controllo == 0
                xprov(end+1)=x(ind_pos(j))+(x(indice_nodo_precedente)-x(ind_pos(j)))...
                    *(ymax-y(ind_pos(j)))/(y(indice_nodo_precedente)-y(ind_pos(j)));
                yprov(end+1)=ymax;
                % se il nuovo nodo trovato è uguale al nodo dentro allora non
                  % aggiungo il nodo due volte, altrimenti aggiungo prima il
                  % nodo nuovo e poi quello successivo valido
                  if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                     xprov(end+1)=x(ind_pos(j));
                     yprov(end+1)=y(ind_pos(j));
                     if ind==1
                            xprov1=[x(ind_pos(j)) xprov(end-1)];
                            yprov1=[y(ind_pos(j)) yprov(end-1)];
                      end
                   else
                     if ind==1
                       xprov1=[x(ind_pos(j))];
                       yprov1=[y(ind_pos(j))];
                     end
                  end
            end
            controllo=0;
        else %se il nodo precedente è dentro
            %non faccio nulla
        end
        
        if nodo_check(indice_nodo_successivo)==0 %se il nodo successivo è fuori
            % verifico in che parte del dominio è fuori
            if x(indice_nodo_successivo)<xmin
                xprov(end+2)=xmin;
                yprov(end+2)=y(ind_pos(j))+(y(indice_nodo_successivo)-y(ind_pos(j)))...
                    *(xmin-x(ind_pos(j)))/(x(indice_nodo_successivo)-x(ind_pos(j)));
                %verifico se il punto della retta interseca il dominio
                if yprov(end)>=ymin && yprov(end)<=ymax
                      % se il nuovo nodo trovato è uguale al nodo dentro allora non
                      % aggiungo il nodo due volte, altrimenti aggiungo prima il
                      % nodo nuovo e poi quello successivo valido
                      if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                       xprov(end-1)=x(ind_pos(j));
                       yprov(end-1)=y(ind_pos(j));
                       if ind==1
                           xprov1=[xprov1(1) xprov(end) xprov1(2:end)];
                           yprov1=[yprov1(1) yprov(end) yprov1(2:end)];
                       end
                       if xprov(end-1)==xprov(end-2) && yprov(end-1)==yprov(end-2)
                         xprov(end-1)=[];
                         yprov(end-1)=[];
                       end
                      else
                        xprov(end-1)=[];
                        yprov(end-1)=[];
                      end     
                    controllo=1;
                else
                    xprov(end-1:end)=[];
                    yprov(end-1:end)=[];
                    controllo=0;
                end
            elseif x(indice_nodo_successivo)>xmax
                xprov(end+2)=xmax;
                yprov(end+2)=y(ind_pos(j))+(y(indice_nodo_successivo)-y(ind_pos(j)))...
                    *(xmax-x(ind_pos(j)))/(x(indice_nodo_successivo)-x(ind_pos(j)));
                %verifico se il punto della retta interseca il dominio
                if yprov(end)>=ymin && yprov(end)<=ymax           
                      % se il nuovo nodo trovato è uguale al nodo dentro allora non
                      % aggiungo il nodo due volte, altrimenti aggiungo prima il
                      % nodo nuovo e poi quello successivo valido
                      if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                       xprov(end-1)=x(ind_pos(j));
                       yprov(end-1)=y(ind_pos(j));
                       if ind==1
                           xprov1=[xprov1(1) xprov(end) xprov1(2:end)];
                           yprov1=[yprov1(1) yprov(end) yprov1(2:end)];
                       end
                       if xprov(end-1)==xprov(end-2) && yprov(end-1)==yprov(end-2)
                         xprov(end-1)=[];
                         yprov(end-1)=[];
                       end
                      else
                        xprov(end-1)=[];
                        yprov(end-1)=[];
                      end            
                    controllo=1;
                else
                    xprov(end-1:end)=[];
                    yprov(end-1:end)=[];
                    controllo=0;
                end
            end
            if y(indice_nodo_successivo)<ymin && controllo == 0
                yprov(end+2)=ymin;
                xprov(end+2)=x(ind_pos(j))+(x(indice_nodo_successivo)-x(ind_pos(j)))...
                    *(ymin-y(ind_pos(j)))/(y(indice_nodo_successivo)-y(ind_pos(j)));
                  % se il nuovo nodo trovato è uguale al nodo dentro allora non
                  % aggiungo il nodo due volte, altrimenti aggiungo prima il
                  % nodo nuovo e poi quello successivo valido
                  if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                     xprov(end-1)=x(ind_pos(j));
                     yprov(end-1)=y(ind_pos(j));
                     if ind==1
                           xprov1=[xprov1(1) xprov(end) xprov1(2:end)];
                           yprov1=[yprov1(1) yprov(end) yprov1(2:end)];
                     end
                     if xprov(end-1)==xprov(end-2) && yprov(end-1)==yprov(end-2)
                      xprov(end-1)=[];
                      yprov(end-1)=[];
                     end
                  else
                     xprov(end-1)=[];
                     yprov(end-1)=[];
                  end   
            elseif y(indice_nodo_successivo)>ymax && controllo == 0
                yprov(end+2)=ymax;
                xprov(end+2)=x(ind_pos(j))+(x(indice_nodo_successivo)-x(ind_pos(j)))...
                    *(ymax-y(ind_pos(j)))/(y(indice_nodo_successivo)-y(ind_pos(j)));
                  % se il nuovo nodo trovato è uguale al nodo dentro allora non
                  % aggiungo il nodo due volte, altrimenti aggiungo prima il
                  % nodo nuovo e poi quello successivo valido
                  if xprov(end)~=x(ind_pos(j)) || yprov(end)~=y(ind_pos(j))
                     xprov(end-1)=x(ind_pos(j));
                     yprov(end-1)=y(ind_pos(j));
                     if ind==1
                           xprov1=[xprov1(1) xprov(end) xprov1(2:end)];
                           yprov1=[yprov1(1) yprov(end) yprov1(2:end)];
                     end
                     if xprov(end-1)==xprov(end-2) && yprov(end-1)==yprov(end-2)
                      xprov(end-1)=[];
                      yprov(end-1)=[];
                     end
                  else
                     xprov(end-1)=[];
                     yprov(end-1)=[];
                  end    
            end
        else %se il nodo successivo è dentro       
            if ind_pos(j)==1
              xprov(end+1)=x(ind_pos(j));
              yprov(end+1)=y(ind_pos(j));  
            end
            if indice_nodo_successivo ~= 1
                xprov(end+1)=x(indice_nodo_successivo);
                yprov(end+1)=y(indice_nodo_successivo); 
            end
        end    
        controllo=0;
    end
    xprov(1)=[];
    yprov(1)=[];
    if ind==3 && abs(y(ind_pos(1))-ymax)<1e-10
        xprov(end)=[];
        yprov(end)=[];
    end
else
    % il poligono ha tutti i nodi dentro
    xprov=x; yprov=y;
end
% se ha solo un nodo dentro il dominio
if ind==1
    xprov=xprov1;
    yprov=yprov1;
end
% caso poligono contenga il punto (xmin,ymin)
if y(1)>ymin && y(1)-2*r_vertice<ymin && (x(1)==xmin || abs(x(1)-xmin)<1e-15)
    xprov=[xprov(1) xmin xprov(2:end)];
    yprov=[yprov(1) ymin yprov(2:end)];
    
end
% caso poligono contenga il punto (xmax,ymin)
if y(1)>ymin && y(1)-2*r_vertice<ymin && (x(1)==xmax || abs(x(1)-xmax)<1e-15)
    xprov(end+1)=xmax;
    yprov(end+1)=ymin;
    xprov(1)=[];
    yprov(1)=[];
end
% caso poligono contenga il punto (xmax,ymax)
if y(1)>ymax && y(1)-2*r_vertice<ymax && (x(1)==xmax || abs(x(1)-xmax)<1e-15)
    xprov(end+1)=xmax;
    yprov(end+1)=ymax;
    xprov(end-1)=[];
    yprov(end-1)=[];
end
if y(1)<ymax && y(1)-2*r_vertice>ymin && (x(1)==xmax || abs(x(1)-xmax)<1e-15)
    xprov(1)=[];
    yprov(1)=[];
    xprov(end)=[];
    yprov(end)=[];
end
end