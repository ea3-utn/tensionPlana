##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    ELEMENTOS BIDIMENSIONALES (Tensión plana)
## /_____/_/  |_/___/___/___/     |
##                                |    quad4LagrangePostPro: PostProcesamientoQuadLagrange 
##                                |                     
##---------CICLO LECTIVO 2021----------------------------------------------------------------

function [deformaciones,estadoTension]=quad4LagrangeOrtoPostPro(CoordX,CoordY,CoordSigmaX,CoordSigmaY,U,base,alto,Elasticidad,poisson,shear,ABSCISAS)

########## ------------- DATOS DEL ELEMENTO
  
  syms x y b h E nu G

  N=[(1-x/b)*(1-y/h) 0 (x/b)*(1-y/h) 0 (1-x/b)*(y/h) 0 (x/b)*(y/h) 0 ;0 (1-x/b)*(1-y/h) 0 (x/b)*(1-y/h)  0 (1-x/b)*(y/h) 0 (x/b)*(y/h)];

  B=[diff(N(1,1),x) 0 diff(N(1,3),x) 0 diff(N(1,5),x) 0 diff(N(1,7),x) 0;0 diff(N(2,2),y) 0 diff(N(2,4),y) 0 diff(N(2,6),y) 0 diff(N(2,8),y);diff(N(1,1),y)+diff(N(2,1),x) diff(N(1,2),y)+diff(N(2,2),x) diff(N(1,3),y)+diff(N(2,3),x) diff(N(1,4),y)+diff(N(2,4),x) diff(N(1,5),y)+diff(N(2,5),x) diff(N(1,6),y)+diff(N(2,6),x) diff(N(1,7),y)+diff(N(2,7),x) diff(N(1,8),y)+diff(N(2,8),x)];

  D=[-E/(nu^2-1) -(E*nu)/(nu^2-1) 0;-(E*nu)/(nu^2-1) -E/(nu^2-1) 0;0 0 2*G];

  MAX=ABSCISAS(3);
  
  Mconectividad= @(u,v) [2*u-1+(v-1)*2*(MAX+1) 2*u+(v-1)*2*(MAX+1) 2*(u+1)-1+(v-1)*2*(MAX+1) 2*(u+1)+(v-1)*2*(MAX+1) 2*(MAX+u)+1+(v-1)*2*(MAX+1) 2*(MAX+u)+2+(v-1)*2*(MAX+1) 2*(MAX+u)+3+(v-1)*2*(MAX+1) 2*(MAX+u)+4+(v-1)*2*(MAX+1)];

########## ------------- DEFORMACIONES

  keyboard

  conectividad=Mconectividad(CoordX,CoordY);

  desplazamientosLocales=U(conectividad);

  deformaciones=subs(B,{b,h},[base alto])*desplazamientosLocales;

  tension=function_handle(D*deformaciones);

  estadoTension=tension(Elasticidad,shear,poisson,CoordSigmaX,CoordSigmaY)
  
endfunction


