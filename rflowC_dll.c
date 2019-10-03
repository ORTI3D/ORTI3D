#include <stdio.h>
#include <stdlib.h>
#include "rflowC_dll.h"
#include <math.h>
#define _Pi 3.1415292654

/*
gcc -c -DBUILDING_RFLOWC_DLL rflowC_dll.c
gcc -shared -o rflowC_dll.dll rflowC_dll.o -Wl,--out-implib,rflowC_dll.a
*/

//////////////////////////// Calc Tubes /////////////////////////:

double calcTube(double *data,double *alph0i,int nx,int ny, double *vxi,
	double *vyi,int nt,int np,double *xpi,double *ypi,double *climi,
	double *xpo,double *ypo,double *tpo,double *cu1o, double *cu2o, int *nto)
{
  /* definitions de base */
  int n,it,i,ip,jp,cl;
  double vx1,vx2,vy1,vy2,Ax,Ay,vxp,vyp;
  int sensx,sensy;

  double dc,alph,lsrc,ys2,fexp,step,talph,calph;
  int sa;
  double a,a1,a2,xp1,yp1;
  double dtx,dty,dt,dxc,dyc,f,xp0,yp0;
  double x1,x2,y1,y2,xf,yf,b1,bet1,bet2,al0;
  double xr,yr;
  int lign;

  /* constante utilisée dans les conditions */
  double eps=1.e-6,cst1=5.e-5,tmp;

  double x0,dx,y0,dy;
  double aT,aL,ul;

  x0=data[0];dx=data[1];y0=data[2];dy=data[3];
  aT=data[4];aL=data[5];ul=data[6];
  
  
  double xg[nx+1],yg[ny+1];
  double xp[nt][np],yp[nt][np],tp[nt][np],cu1[nt][np],cu2[nt][np];
  int ptin[np];
  double alph0[np];

  /* transferer depuis python */
  for (i = 0; i < np; i++) {alph0[i] = alph0i[i];}

  for (i = 0; i < np; i++) {
    xp[0][i] = xpi[i];yp[0][i] = ypi[i];
    tp[0][i] = 0.;cu1[0][i] = 0.;cu2[0][i] = 0.;
  }

  /* creer xg et yg */
  for (i = 0; i<(nx+1); i++) {xg[i] = x0+i*dx;}
  for (i = 0; i<(ny+1); i++) {yg[i] = y0+i*dy;}

  /* debut des calculs */
  dc = dx/2.+dy/2.;
  a1 = xp[0][0]-xp[0][np-1];a2 = yp[0][0]-yp[0][np-1];
  lsrc = sqrt(a1*a1+a2*a2);

  ys2=lsrc/sqrt(aT)+1.*exp(-lsrc/sqrt(aT)/10.); /*modif v202 pour petit source*/
  fexp=(1.6+2.5/ys2)*.94;
  step=-60.+25.*exp(1.1*log(ys2/2.))+aL;

   /* Boucle generale */
  for (n=0;n<np;n++) // num de particule
    {
    jp = (int)(floor((xp[0][n]-x0)/dx));
    ip = (int)(floor((yp[0][n]-y0)/dy));
    a1=vxi[ip*nx+jp];a2=vyi[ip*nx+jp];
    b1=a1-a2*tan(alph0[n])*lsrc/2.;sensx = (b1>0.) ? 1:-1;
    b1=a2+a1*tan(alph0[n])*lsrc/2.;sensy = (b1>0.) ? 1:-1;
    lign=0;  // 1 : vertical, -1 horizontal

    for (it=1;it<nt;it++)
      {
      dxc = 0.; dyc = 0.;dt = 0.; tmp = 0.; calph=0.;
      cl = climi[ip*nx+jp];
      ptin[n] =((xp[it-1][n]<x0+dx/2.)||(xp[it-1][n]>xg[nx]-dx/2.)||
	   (yp[it-1][n]<y0+dy/2.)||(yp[it-1][n]>yg[ny]-dy/2.)||(cl==0.))? 0:1; //
      if (it<10) {ptin[n]=1;}      // pour le debut
 
      if (ptin[n]==1)
      {
        xp0=xp[it-1][n];yp0=yp[it-1][n];
          
	/* calculer dispersion en passant par formule ~cleary (sauf sur les bords) */
	alph = alph0[n]*erfc(log10(maxD(cu1[it-1][n]+dc/2.+step,1.))/fexp);
        sa = (alph>0.)? 1:-1;sa = (alph==0.)? 0:sa;
	talph=tan(alph)*lsrc/2.;calph=cos(atan(tan(alph)*lsrc/2.));
          
	// vitesses
        if (lign==1) {
           ip=(int)(floor((yp[it-1][n]-y0)/dy));
           sa=(sensx<0)?-1:0;jp=(int)(floor((xp[it-1][n]-x0)/dx)+sa); }
        else if (lign==-1) {
           jp=(int)(floor((xp[it-1][n]-x0)/dx));
           sa=(sensy<0)?-1:0;ip=(int)(floor((yp[it-1][n]-y0)/dy)+sa); }               

	vx1 = vxi[ip*nx+jp];vx2 = vxi[ip*nx+jp+1];      
	vy1 = vyi[ip*nx+jp];vy2 = vyi[(ip+1)*nx+jp];
        vx1=vx1+eps;vx2=vx2+eps;vy1=vy1+eps;vy2=vy2+eps;//eps pour eviter 0. 
	Ax = (vx2-vx1)/dx;x1=xg[jp]; x2=x1+dx;Ax=(Ax==0.)?1.e-9:Ax;
	Ay = (vy2-vy1)/dy;y1=yg[ip]; y2=y1+dy;Ay=(Ay==0.)?1.e-9:Ay;
 	vxp = vx1+Ax*(xp0-x1);
	vyp = vy1+Ay*(yp0-y1);
        sa=(sensx>0)?1:0;xf=xg[jp+sa];
        sa=(sensy>0)?1:0;yf=yg[ip+sa];

	  /* on differencie les cas */
	if ( (fabs(vy1)+fabs(vy2)) < ((fabs(vx1)+fabs(vx2))*cst1) ) // sens x
	   {
           dtx = fabs((xf-xp0)/vx1);
           dty = fabs((yf-yp0)/maxD(vx1*fabs(talph),eps));
           dt = dtx*(1.+eps); lign=1;
           if ((dty<dtx)&&(dty!=0.)) { dt = dty*(1.+eps);lign=-1;}
           dxc = vx1*dt;dyc = dxc*talph;
	   } 
	else if ( (fabs(vx1)+fabs(vx2)) < ((fabs(vy1)+fabs(vy2))*cst1) ) // sens y
	   {
           dty = fabs((yf-yp0)/vy1);
           dtx = fabs((xf-xp0)/maxD(vy1*fabs(talph),eps));
           dt = dty*(1.+eps); lign=-1;
           if ((dtx<dty)&&(dtx!=0.))  {dt = dtx*(1.+eps);lign=1;}
           dyc = vy1*dt;dxc = -dyc*talph;
           }
	else 
	   {
           al0=sqrt(fabs(Ay*Ay-2*Ax*Ay+4*Ax*Ay*talph*talph+Ax*Ax));
    	   bet1=-vx1/Ax+x1-xp0;bet2=vy1/Ay+yp0-y1;
	   a1=(Ay-Ax)*bet1-2*Ay*talph*bet2;
	   a2=(Ay-Ax)*bet2-2*Ax*talph*bet1;

           dt=minD(dx/fabs(vxp),dy/fabs(vyp));
	   xr=exp((Ax+Ay)*dt/2)*(sinh(al0*dt/2)*a1/al0-bet1*cosh(al0*dt/2))+x1-vx1/Ax;
           yr=exp((Ax+Ay)*dt/2)*(sinh(al0*dt/2)*a2/al0+bet2*cosh(al0*dt/2))+y1-vy1/Ay;
           sa=(xr-xp0)>0. ? 1:0; xf=xg[jp+sa];
           sa=(yr-yp0)>0. ? 1:0; yf=yg[ip+sa];
           sa=(yr-yp0)>0. ? 1:-1;
           if ((sa!=sensy)&&(lign==-1)) {yf=yg[ip+(int)(.5+1.5*sa)];}
           dtx=fabs((xf-xp0)/(xr-xp0)*dt);dty=fabs((yf-yp0)/(yr-yp0)*dt);

	   if (dtx<dty) 
 	     {
             lign=1;dxc=(xf-xp0);
	     xp1=exp((Ax+Ay)*dtx/2)*(sinh(al0*dtx/2)*a1/al0-bet1*cosh(al0*dtx/2))+x1-vx1/Ax;
             dtx=maxD(dtx+(xf-xp1)/(xp1-xr)*(dtx-dt),0.);
	     xp1=exp((Ax+Ay)*dtx/2)*(sinh(al0*dtx/2)*a1/al0-bet1*cosh(al0*dtx/2))+x1-vx1/Ax; // ajout pour SeymBug
             dt=maxD(dtx+(xf-xp1)/(xp1-xr)*(dtx-dt),0.);
             yf=exp((Ax+Ay)*dt/2)*(sinh(al0*dt/2)*a2/al0+bet2*cosh(al0*dt/2))+y1-vy1/Ay;
             dyc=(yf-yp0)*(1+eps);
	     }
           else
	     {
             lign=-1;dyc=(yf-yp0);
             yp1=exp((Ax+Ay)*dty/2)*(sinh(al0*dty/2)*a2/al0+bet2*cosh(al0*dty/2))+y1-vy1/Ay;
             dty=maxD(dty+(yf-yp1)/(yp1-yr)*(dty-dt),0.);
             yp1=exp((Ax+Ay)*dty/2)*(sinh(al0*dty/2)*a2/al0+bet2*cosh(al0*dty/2))+y1-vy1/Ay;
             dt=maxD(dty+(yf-yp1)/(yp1-yr)*(dty-dt),0.);
             xf=exp((Ax+Ay)*dt/2)*(sinh(al0*dt/2)*a1/al0-bet1*cosh(al0*dt/2))+x1-vx1/Ax;
             dxc=(xf-xp0)*(1+eps);
	     }
	   }  // fin des differents cas ds cellule
          sensx=(dxc>0.) ? 1:-1;sensy=(dyc>0.) ? 1:-1;
          tmp = sqrt(dxc*dxc+dyc*dyc);f=1.;
          }  // fin du cas ptin <>0

	  else {f = 0.;} //ntmx=(nt>ntmx)? nt:ntmx;break;}  // cas ou ptin=0

	  /* mis a jour des matrices */
	  cu1[it][n] = cu1[it-1][n] + calph*tmp*f;
	  cu2[it][n] = cu2[it-1][n] + tmp*f;
	  xp[it][n] = xp[it-1][n] + dxc*f;
	  yp[it][n] = yp[it-1][n] + dyc*f;
	  tp[it][n] = tp[it-1][n] + dt*f;


	}  // fin de la boucle sur pas de temps it
    }  // fin de la boucle sur particules (n) 

int j; nto[0]=nt;
for (i=0; i<nt; i++)
   for (j=0; j<np; j++)
      {
      xpo[i*np+j] = xp[i][j];ypo[i*np+j] = yp[i][j];
      cu1o[i*np+j] = cu1[i][j];cu2o[i*np+j] = cu2[i][j];
      tpo[i*np+j] = tp[i][j];
      }
 
}  // fin de la fonction

///////////////////////////////////////////////////////////
//////////////// Pts Ligness pour TUBES ///////////////////

double ptsLigne(int *indi,double *teti, double *dci,
	double *xpi,double *ypi,double *tpi,double *cuai,double *cubi,
	double *xpo,double *ypo,double *tpo,double *cuao, double *cubo)
{
  int i,i2,j,nt,np,ip0,ip1,ntmx,nt1;
  double a,b,d,e,f,tet,dc;
  ntmx=indi[0];nt=indi[1];np=indi[2];
  double xp1[ntmx][np],yp1[ntmx][np],tp1[ntmx][np],cu1[ntmx][np],cu2[ntmx][np];

    for (j=0;j<np;j++) //interpolation lineaire, par particule
    { 
      i2 = 1;ip0=1;dc=dci[j];
      xp1[0][j]=xpi[j];yp1[0][j]=ypi[j];
      tp1[0][j]=tpi[j];cu1[0][j]=cuai[j];cu2[0][j]=cubi[j];
      for (i=1;i<nt;i++)
      {  
        ip1=(int)(ceil(cuai[i*np+j]/dc));
        while (ip1>ip0) {
            b=cuai[(i-1)*np+j];d=dc*ip0-b;
            a=cuai[i*np+j];if(a==b){break;}
		f=a-b;e=f-d;
            xp1[i2][j]=xpi[(i-1)*np+j]*e/f+xpi[i*np+j]*d/f; //
            yp1[i2][j]=ypi[(i-1)*np+j]*e/f+ypi[i*np+j]*d/f;
            tp1[i2][j]=tpi[(i-1)*np+j]*e/f+tpi[i*np+j]*d/f;
            cu1[i2][j]=cuai[(i-1)*np+j]*e/f+cuai[i*np+j]*d/f;
            cu2[i2][j]=cubi[(i-1)*np+j]*e/f+cubi[i*np+j]*d/f;
            i2=(i2>ntmx-2)? ntmx-1:i2+1;
	    ip0=ip0+1;
        }  // fin du while
      }  // fin boucle for i; debut boucle pour remplir fin de l'array

      tet=teti[0]+(j-np/2)*.004; //
      dc=xp1[i2-1][j]-xp1[i2-2][j];
      for (i=i2-1;i<ntmx;i++) 
      { 
          xp1[i][j]=xp1[i-1][j]+dc*cos(tet);
          yp1[i][j]=yp1[i-1][j]+dc*sin(tet);
          tp1[i][j]=tp1[i-1][j];
          cu1[i][j]=cu1[i-1][j]+dc;
	  cu2[i][j]=cu2[i-1][j];
      }   // fin boucle i

    }  // fin boucle j (particules)
for (i=0; i<ntmx; i++)
 for (j=0; j<np ; j++)
      {
	xpo[i*np+j] = xp1[i][j];ypo[i*np+j] = yp1[i][j];tpo[i*np+j] = tp1[i][j];
	cuao[i*np+j] = cu1[i][j];cubo[i*np+j] = cu2[i][j];
      }

}  // fin de la fonction

///////////////////////  CALC PARTICULES    /////////////////////////  
////////////////////////////////////////////////////////////////////////  

double calcPart(int *dimsi, double *datai,double *vxi,double *vyi,double *climi,
	double *xpo,double *ypo,double *tpo,double *cu1o, int *ito)

{
  /* definitions de base */
  int it,i,ip,jp,ptin;
  double vx1,vx2,vy1,vy2,Ax,Ay,vxm,vym,x0m,y0m,vxp0,vyp0;
  int sensx,sensy;
  double dc,lb1,lb2,ax1,ay1,ax2,ay2;
  double dtx,dty,dtx1,dty1,dtx2,dty2,dt,dxc,dyc;

  /* constante utilisée dans les conditions */
  double eps=1.e-9,cst1=5.e-5;

  double x0,y0,xp0,yp0,cl,dx,dy;
  x0=datai[0];dx=datai[1];y0=datai[2];dy=datai[3];xp0=datai[4];yp0=datai[5];
  int nx,ny,ndx,ndy,nt;
  nx=dimsi[0];ny=dimsi[1];ndx=dimsi[2];ndy=dimsi[3];nt=dimsi[4];
  
  double xg[nx+1],yg[ny+1];
  double xp[nt+1],yp[nt+1],tp[nt+1],cu1[nt+1];

  /* creer xg et yg */
  for (i = 0; i<nx+1; i++) {
	// if (ndx==1) {dx=dxi[0];} else {dx=dxi[i];}
	xg[i] = x0+i*dx;
	}
  for (i = 0; i<ny+1; i++) {
	// if (ndy==1) {dy=dyi[0];} else {dy=dyi[i];}
	yg[i] = y0+i*dy;
	}

  /* debut des calculs */
  dc = dx/2; // dx=dxi[0];dy=dyi[0];
  it = 0;xp[0]=xp0;yp[0]=yp0;tp[0]=0.;ptin=1;
  while ((it<nt)&&(ptin==1))
  {
      it++;
      dxc = 0.; dyc = 0.;dt = 0.;
      jp = (int)(floor((xp[it-1]-x0)/dx));
      ip = (int)(floor((yp[it-1]-y0)/dy));
      cl = climi[ip*nx+jp];
      if ((jp<(nx))&&(jp>0)&&(ip<(ny))&&(ip>0)&&(cl>0)||(it<5))
      {ptin=1;}
      else {ptin=0;}
      vx1 = vxi[ip*nx+jp];vx2 = vxi[ip*nx+jp+1];      
      vy1 = vyi[ip*nx+jp];vy2 = vyi[(ip+1)*nx+jp];
      Ax = (vx2-vx1)/dx;vxm = 2*vx1-vx2;
      Ay = (vy2-vy1)/dy;vym = 2*vy1-vy2;
      
      x0m = (xp[it-1]-xg[jp]+dx);
      y0m = (yp[it-1]-yg[ip]+dy);
      vxp0 = vxm+Ax*x0m;
      vyp0 = vym+Ay*y0m;
      sensx = (vxp0>=0.) ? 1:-1;         
      sensy = (vyp0>=0.) ? 1:-1;

      /* on differencie les cas */
      if ( (fabs(vy1)+fabs(vy2)) < ((fabs(vx1)+fabs(vx2))*cst1) ) // sens x
         {
          dt = ((1.5+0.5*sensx)*dx-x0m)/vxp0*(1.+eps);
          dxc = vxp0*dt; dyc = 0;
         } 
      else if ( (fabs(vx1)+fabs(vx2)) < ((fabs(vy1)+fabs(vy2))*cst1) ) // sens y
         {
          dt = ((1.5+0.5*sensy)*dy-y0m)/vyp0*(1.+eps);
          dyc = vyp0*dt; dxc = 0;
         }
      else 
         {
          lb1 = (vxp0-vxm)/x0m;
          lb2 = (vyp0-vym)/y0m;
          ax1 = maxD((lb1*dx+vxm)/vxp0,eps);
          ax2 = maxD((lb1*dx*2+vxm)/vxp0,eps);
          ay1 = maxD((lb2*dy+vym)/vyp0,eps);
          ay2 = maxD((lb2*dy*2+vym)/vyp0,eps);
          dtx1 = log(ax1)/lb1;dtx2 = log(ax2)/lb1;dtx = maxD(dtx1,dtx2);
          dty1 = log(ay1)/lb2;dty2 = log(ay2)/lb2;dty = maxD(dty1,dty2);
          if (dtx<=0) { dtx=1.e5; }
          if (dty<=0) { dty=1.e5; }
          dt = minD(dtx,dty)*(1.+eps);
          dxc = ( vxp0*exp(lb1*dt)-vxm )/lb1-x0m; 
          dyc = ( vyp0*exp(lb2*dt)-vym )/lb2-y0m;
         }
      /* mis a jour des matrices */
      cu1[it] = cu1[it-1] + sqrt(dxc*dxc+dyc*dyc);
      xp[it] = xp[it-1] + dxc;
      yp[it] = yp[it-1] + dyc;
      tp[it] = tp[it-1] + dt;
      
    } // fin du while sur it

ito[0]=it;
for (i=0; i<it; i++)
      {
      xpo[i] = xp[i];ypo[i] = yp[i];cu1o[i] = cu1[i];tpo[i] = tp[i];
      }
  
}  // fin de la fonction calcParticules

///////////////////  Pts Lignes2 pour 1 particule  ////////////////////////

double ptsLineP(int *dimsi,double *xpi,double *ypi,double *tpi,double *dti,
	double *zpi,double *xpo,double *ypo,double *tpo,double *zpo,int *ito)
{
  int i,it,i2,nt,ip0,ip1,ntmx;
  double a,b,d,e,f,tet,dt;
  nt=dimsi[0];ntmx=100;dt=dti[0];
  double xp1[ntmx],yp1[ntmx],tp1[ntmx],zp1[ntmx];
  //interpolation lineaire, sur une particule
      i2 = 1;ip0=1;
      xp1[0]=xpi[0];yp1[0]=ypi[0];tp1[0]=tpi[0];zp1[0]=zpi[0];ito[0]=0;
      for (it=1;it<nt;it++) {
        ip1=(int)(ceil(tpi[it]/dt));
        while (ip1>ip0) {
            b=tpi[it-1];d=dt*ip0-b;
            a=tpi[it];if(a==b){break;}
		f=a-b;e=f-d;
            xp1[i2]=xpi[it-1]*e/f+xpi[it]*d/f;
            yp1[i2]=ypi[it-1]*e/f+ypi[it]*d/f;
            tp1[i2]=tpi[it-1]*e/f+tpi[it]*d/f;
            zp1[i2]=zpi[it-1]*e/f+zpi[it]*d/f;ito[i2]=it-1;
            i2 = minD(i2 +1,ntmx-1);ip0=ip0+1;
        }  // fin du while
      }   // fin boucle it
//ito[0]=i2;
for (i=0; i<i2; i++)
      {xpo[i] = xp1[i];ypo[i] = yp1[i];tpo[i] = tp1[i];zpo[i] = zp1[i];}

}  // fin de la fonction

//////////////////////////  INTERP2D     /////////////////////////////////
double interp2d(int *dimsi, double *datai, 
	double *xpi,double *ypi,double *zpi,double *zpo)

{
  /* definitions de base */
  int i,ip,ix,iy;
  double x0,y0,dx,dy,distx,disty;
  x0=datai[0];y0=datai[1];dx=datai[2];dy=datai[3];
  distx=datai[4];disty=datai[5];
  int nx,ny,np;
  nx=dimsi[0];ny=dimsi[1];np=dimsi[2];
  
  double xg[nx],yg[ny];

  /* creer xg et yg */
  xg[0]=dx/2;yg[0]=dy/2;
  for (i = 1; i<nx; i++) {
	// if (ndx==1) {dx=dxi[0];} else {dx=dxi[i];}
	xg[i] = x0+(i-1)*dx+dx/2;
	}
  for (i = 0; i<ny; i++) {
	// if (ndy==1) {dy=dyi[0];} else {dy=dyi[i];}
	yg[i] = y0+(i-1)*dy+dy/2;
	}

  double xp,yp,zp,dist,count[ny][nx];

// creation mise a zero des variables
  for (ix = 0; ix < nx; ix++) 
    for (iy = 0; iy < ny; iy++) {
        zpo[iy*nx+ix]=0.;count[iy][ix]=0.;
  }

// calcul
  for (ip=0;ip<np;ip++)
   for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
        dist=(xg[ix]-xpi[ip])*(xg[ix]-xpi[ip])/distx;
        dist=dist+(yg[iy]-ypi[ip])*(yg[iy]-ypi[ip])/disty;
        zpo[iy*nx+ix]=zpo[iy*nx+ix]+exp(-sqrt(dist))*zpi[ip];
        count[iy][ix]=count[iy][ix]+exp(-sqrt(dist));
  }
   for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
        zpo[iy*nx+ix]=zpo[iy*nx+ix]/count[iy][ix];
  }

}  // fin de la fonction

///////////////////////////////////////////////////////////////////

double maxD(double a, double b) {return (a > b) ? a : b;}

double minD(double a, double b) {return (a < b) ? a : b;}
