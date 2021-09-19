##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ELEMENTOS BIDIMENSIONALES (Tensión plana)
## /_____/_/  |_/___/___/___/     |
##                                |    MEFtensionPlana: script principal
##---------CICLO LECTIVO 2021----------------------------------------------------

## CONFIGURACION

clear

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

syms x y

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

t=2*.762e-3;

b=508e-3;

alto=152.4e-3;

E=70e9;

G=5e9;

#nu=.33;

nu=.1;

P=1000;

h=17.8e-3;

q=P/(t*h);

qUnitaria=q/t;

## GEOMETRIA Y DISCRETIZACION

ABSCISAS=[0 b 2 E nu G];% [Xi Xf Discretizacion E nu G]

ORDENADAS=[0 alto 2];% [Xi Xf Discretizacion]

## CARGAS Y CONDICIONES DE CONTORNO

CCx=[3 0;6 0;9 0]; # [Ni CC]

CCy=[3 0;6 0;9 0]; # [Ni C]C

CARGAx=[]; # [Ni CC]

CARGAy=[]; # [Ni CC]


###---- Caracteristicas del elemento

anchoElemento=b/ABSCISAS(3);

altoElemento=alto/ORDENADAS(3);

###---- Cargas distribuidas

#### [coordenadaElemento OrdenadaElemento yi yf Xo Qy]

% 1  Coordenadas del elemento al que se le esta aplicando la carga
% 2  yi yf --> Rango de integracion de la carga EN COORDENADAS LOCALES DEL ELEMENTO
% 3  Xo --> De que lado se aplica la carga, puede valer 0 o el Ancho/alto DEL ELEMENTO (Recordar que son coordenadas locales)

distribuidaX=[1 1 altoElemento-h/2 altoElemento 0 -q;1 2 0 h/2 0 -q]; #   -------> Tiene que estar en ceros si no se usa

distribuidaY=[0 0 0 0 0 0]; # [xi xf Yo Qx] -------> Tiene que estar en ceros si no se usa

###---- Post-Procesamiento

CoX=1; # Coordenada X del elemento de interes

CoY=1;

CsiX=0; # Coordenada LOCAL del elemento para obtener estado de tension

CsiY=altoElemento;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SCRIPT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#[KG,fq]=quad4Lagrange(ABSCISAS,ORDENADAS,t,distribuidaX,distribuidaY);

[KG,fq]=quad4LagrangeOrto(ABSCISAS,ORDENADAS,t,distribuidaX,distribuidaY);

keyboard

GL=size(KG,1); % Cant. de grados de libertad globales



[P,U]=vectorCargas(GL,CARGAx,CARGAy,CCx,CCy);



				% Guyan

INDEX=linspace(1,GL,GL);

i=1;guyan=[1 1];

CC=sum(isnan(U));

				#while (isempty(guyan)<1 && i<CC+1)
while (i<CC+1)
  
  Null=INDEX(isnan(U));

  Zeros=find(U==0);

  if (Null(i)>Zeros(1))

    guyan=[Null(i) Zeros(1)];

    try

      registroGuyan=[registroGuyan;guyan];

    catch

      registroGuyan=[guyan];

    end_try_catch

    [KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,guyan,1,2);
        
  endif

  i++;
  
endwhile


				% RESOLUCION



K11=KG(1:CC,1:CC);

K12=KG(1:CC,CC+1:GL);

K21=KG(CC+1:GL,1:CC);

K22=KG(CC+1:GL,CC+1:GL);

PII=P(1:CC,:); % Cargas conocidas

Fl=fq(1:CC,:); % Cargas locales asociadas a P conocidas

Fp=fq(CC+1:end,:); % Cargas locales NO asociadas a P conocidas

UI=U(CC+1:GL,:); % Desplazamientos conocidos

determinante=det(K11);

if determinante==0

  UII=zeros(size(PII));

else
  
  UII=inv(K11)*(PII-Fl-K12*UI);  % Desplazamientos desconocidos

endif


PI=K21*(UII)+K22*UI+Fp; % Cargas desconocidas


				% REARMADO DE RESULTADOS

j=1;h=j;i=j;

while (j<=size(P,1))

  
  
  if isnan(P(j))==1

    P(j)=PI(i);

    i++;
    
  endif

  if isnan(U(j))==1

    U(j)=UII(h);

    h++;
    
  endif

  j++;
  
endwhile


[KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,registroGuyan,2,1);

				% VERIFICACION

FX=sum(P(1:2:end))-sum(fq(1:2:end)) % sumatoria en X

FY=sum(P(2:2:end))-sum(fq(2:2:end)) % sumatoria en Y


####### ---------- POST PROCESADO


#[tension]=quad4LagrangePostPro(CoX,CoY,CsiX,CsiY,U,anchoElemento,altoElemento,E,nu,ABSCISAS);

[tension]=quad4LagrangeOrtoPostPro(CoX,CoY,CsiX,CsiY,U,anchoElemento,altoElemento,E,nu,G,ABSCISAS);

aster=dlmread("matrizRigidez.csv")*t;

sort(diag(KG))

sort(diag(aster))

keyboard
