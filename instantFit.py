import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.special import erf,erfc
from scipy.stats import norm
from .config import *
import os,time
#from rflow import *

class instantFit():
    
    def __init__(self,core):
        self.core,self.gui = core,core.gui
        self.K,self.H = None, None
        self.cfg = Config(core)
        #self.rflow = rflow(core)

    def calculate(self,dic_options,obs_value):
        '''make the calcualtion according to data present in dic_options
        and uses obs_value to know from which dialog the info is coming (not used now)'''
        # get current zone and changes its value
        #print obs_value,dic_options
        grd = self.core.dicaddin['Grid'];#print grd
        self.domain = [grd['x0'],grd['x1'],grd['y0'],grd['y1']]
        self.dx, self.dy = float(grd['dx'][0]),float(grd['dy'][0])
        self.nrow,self.ncol = len(grd['dy']),len(grd['dx'])
        self.K = self.core.getValueLong('Modflow','lpf.8',0)[0]*1;#print 'K from core',amin(self.K)  # here its 2D
        self.Darcy()
        if dic_options['type']=='Tracer':
            self.aL,self.aT = dic_options['aL'],dic_options['aT']
            self.Sourcexy = self.core.diczone['Mt3dms'].getValue('btn.13','coords',0)
            self.SourceC = float(self.core.diczone['Mt3dms'].getValue('btn.13','value',0))
            xp1,yp1,tp1,cu1,cu2,V,Large,Cinf = self.calcTubes2D();#print shape(xp1),shape(Cinf)
            #to put in the cneter of the cells
            xp1, yp1 = (xp1[:,1:]+xp1[:,:-1])/2,(yp1[:,1:]+yp1[:,:-1])/2
            cu1, tp1 = (cu1[:,1:]+cu1[:,:-1])/2,(tp1[:,1:]+tp1[:,:-1])/2
            t = float(self.gui.guiShow.getCurrentTime())
            C = self.calcTime(xp1,yp1,tp1,cu1,Cinf,t);#print C
            self.xp1,self.yp1,self.C = xp1,yp1,C
        
    def Darcy(self):
        '''K is the Hyd cond matrix, H_bc is a matrix of initial H values with BC on sides
        size of marices : K (nl, nc), bc nl or nc, H nl-2, nc-2 : no bc cells
        the right side is q (=0), except at the boundaries where it is H_bc
        and K=1 to satisfy the BCs'''
        #pass
        K,H_bc = self.K,self.H_bc;#H_bc contains only value at the fixed heads
        nr,nc = shape(K)
        ncell = nc*nr
        Kr = sqrt(K[:-1,:]*K[1:,:])/self.dy; # K interpolated in y dir top-bottom (rows)
        Kc = sqrt(K[:,:-1]*K[:,1:])/self.dx; # K interpolated in x dir left-right (columns)
        Kc0, Kc1 = c_[zeros((nr,1)),Kc],c_[Kc,zeros((nr,1))]
        Kr0, Kr1 = r_[zeros((1,nc)),Kr],r_[Kr,zeros((1,nc))]
        # if fixed head on BC put just one in sK and 0 elsewhere
        Kr0[H_bc>0],Kr1[H_bc>0],Kc0[H_bc>0],Kc1[H_bc>0] = 0,0,0,0
        sK = Kc0+Kc1+Kr0+Kr1
        sK[H_bc>0]=1
        lKc0=r_[ravel(-Kc0)[1:],zeros(1)]
        lKc1=r_[zeros(1),ravel(-Kc1)[:-1]]
        lKr0=r_[ravel(-Kr0)[nc:],zeros(nc)]
        lKr1=r_[zeros(nc),ravel(-Kr1)[:-nc]]
        lK = [ravel(sK),lKc0,lKc1,lKr0,lKr1]
        diags=[0,-1,1,-nc,nc]
        KM = sp.spdiags(lK,diags,ncell,ncell) # put Kv in first diagonal
        Hout = spsolve(sp.csr_matrix(KM),ravel(H_bc));
        H=reshape(Hout,(nr,nc))
        self.Q = [-sqrt(K[:,:-1]*K[:,1:])*(H[:,1:]-H[:,:-1])/self.dx,
                  -sqrt(K[:-1,:]*K[1:,:])*(H[1:,:]-H[:-1,:])/self.dy]
        self.H = H

    def calcTubes2D(self,opt=0,ls1=None):
        """ methode pour calculer l'avancement des particules, necessite la position
        de la source, les dimensions de la grille et les vitesses, contenues dans model
        et les donnees de erfmat.
        La focntion retourne la position des particules pour un temps infini :
        XP0,YP0 les coordonnees des particules
        TP et CU: temps de la particule en chaque point XP,YP et CU distance depuis la source
        ???on ne tient pas compte de variation epaissseur dans calc vitesses
        donc pas necessaire ensuite pour calc de concentration???"""
        #initialisation des valeurs
    #if 3>2:
        nr,nc = shape(self.K)
        q = self.Q # darcy velocities
        x0,x1,y0,y1 = [float(a) for a in self.domain]
        dx,dy = self.dx,self.dy
        aL,aT  =  self.aL,self.aT
        poro = 0.25
        vx,vy = q[0][:,:-1]/poro, q[1][:-1,:]/poro # to have arrays of the same size
        xy  =  self.Sourcexy;x, y = list(zip(*xy)); x=list(x); y=list(y);
        yp0 = erfc(-0.08*arange(26))-1
        #yp0 = erfc(-0.06*arange(18))-1
        yp0 = yp0*0.999/yp0[-1]
        #yp0=array([0,.08,.2,.3,.4,.5,.58,.65,.7,.75,.8,.85,0.89,.92,.94,.96,0.985,0.999])
        # inverse the source if necessary
        jp = ((x[0]+x[1])/2-x0)/dx;jp =int(jp);jp=max(jp,0)
        ip = ((y[0]+y[1])/2-y0)/dy;ip =int(ip);ip=max(ip,0)
        vx1 = vx[ip,jp]; vy1 = vy[ip,jp]
        if abs(vx1)>abs(vy1) :
            if sign(y[0]-y[1])!=sign(vx1):
                x.reverse(); y.reverse()
        else :
            if sign(x[1]-x[0])!=sign(vy1):
                x.reverse(); y.reverse()
        # assemble data for the calculation
        data = (x0,dx,y0,dy,aT,aL,0.);# print data
        ypi=r_[-yp0[-1:0:-1],yp0]/2.
        xpi=ypi*0.;
        xpin = (x[0]+x[1])/2.+ypi*(x[0]-x[1]) # 
        ypin = (y[0]+y[1])/2.+ypi*(y[0]-y[1]);
        np=len(xpin)
        clm=ones((nr,nc));#clm[1:-1,1:-1]=1
        #tstart=time.time()
        alph0=-sign(ypi)*norm.isf(.5+abs(ypi),0.,.25); # distribution of starting point along the source line
        #alph0[abs(alph0)>.5]=alph0[abs(alph0)>.5]*1.5 # pour avoir conc faibe plus sur le cote
        alph0=(exp(abs(alph0))-1)*sign(alph0) # correction for the sides
        [xp0,yp0,tp0,cua,cub]=self.rflow.calcTube1(data,alph0,vx,vy,xpin,ypin,clm);
        #xp0,yp0,tp0,cua,cub = array(xp0),array(yp0),array(tp0),array(cua),array(cub)
        nt,np = shape(xp0)
        indx = ones(nt).astype('int')
        for it in range(nt):
            if sum((xp0[it,:]-xp0[it-1,:])**2+(yp0[it,:]-yp0[it-1,:])**2)==0:
                indx[it] = 0
        xp0,yp0,tp0,cua,cub = xp0[indx>0,:],yp0[indx>0,:],tp0[indx>0,:],cua[indx>0,:],cub[indx>0,:]        
        #print 'xp,yp,tp,cua,cub',xp0,yp0,tp0,cua,cub
        [xp1,yp1,tp1,cu1,cu2,V,Large,Cinf] = self.calcConc(xp0,yp0,tp0,cua,cub); 
        #dt=time.time()-tstart;#print 'xp1b',xp1[0,:]#print dt
        return [self.smoo(xp1),self.smoo(yp1),tp1,cu1,cu2,V,Large,self.smoo(Cinf)]
    
    def smoo(self,v):
        """ fonction qui permet de lisser une matrice par des moyennes mobiles
        """
        [l,c] = shape(v)
        v1 = concatenate([v[:1,:],(v[:-2,:]+v[1:-1,:]+v[2:,:])/3,v[-1:,:]],axis=0)
        v2 = concatenate([v1[:,:1],(v1[:,:-2]+v[:,1:-1]+v1[:,2:])/3,v1[:,c-1:c]],axis=1)
        return v2
        
    def getCfromReglist(self,ix,iy):
        """ix,iy are list of indices referenced on the regular grid, 
        the routines returns the C value at the closest
        points belonging to the Cgrid"""
        n = len(ix)
        x,y = array(ix)*self.dx,array(iy)*self.dy
        clist = []
        for i in range(n):
            d = sqrt((x[i]-self.xp1)**2+(y[i]-self.yp1)**2)
            iy1,ix1=where(d==amin(d))
            clist.append(self.C[iy1[0],ix1[0]])
        return clist
        
    def calcTime(self,xp1,yp1,tp1,cu1,Cinf,t):
        """starting from concentraitons at an infinite time, calculate conc for a given time
        inj : liste [periodes,conc], k cinet 1o, rf retard
        tf temps de la simulation
        """
        putmask(tp1,tp1<=0,1)
        vmoy=cu1/tp1;putmask(vmoy,vmoy<=0,1)
        t2=t+self.aL/200*(1-exp(-3*t));  #a correciton factor 
        a = sqrt(4*self.aL*vmoy*t);
        c1 = erfc((cu1-vmoy*t2)/a)
        c2 = exp(minimum(cu1/self.aL,500))*erfc((cu1+vmoy*t2)/a)
        c3 = (c1+c2)*Cinf/2.;
        return c3
    
    def calcConc(self,xp0,yp0,tp0,cua,cub):
        '''calculate the concentraiton along the lines, using the distance between these lines
        start from the middle and find the point that is on the next line and perpendicular
        '''
        xp0,yp0 = self.smoo(xp0),self.smoo(yp0)
        cua,tp0 = self.smoo(cua),self.smoo(tp0)
    #if 3>2:       
        self.flg = ((xp0[:-1,:]-xp0[1:,:])**2+(yp0[:-1,:]-yp0[1:,:])**2)<1e-10 # true if pts are the same at teh end of the lines
        nt,np = shape(xp0);
        xp1,yp1,vloc = zeros((nt,np)),zeros((nt,np)),zeros((nt,np))
        vloc[0,:] = (cua[1,:]-cua[0,:])/(tp0[1,:]-tp0[0,:])
        self.xp0,self.yp0,self.cua,self.tp0,self.vloc,self.nt = xp0,yp0,cua, tp0, vloc,nt
        # finds the coefficient for the lines in xp0,yp0
        self.lcoefs= [];a,b = str(xp0),str(yp0)
        for ip in range(np): 
            self.lcoefs.append(self.getLcoefs(list(zip(xp0[:,ip],yp0[:,ip]))))
        # then for each point find the width
        wdth,mdle = zeros((nt,np-1)),int((np-1)/2)
        wdth[0,:] = sqrt((xp0[0,:-1]-xp0[0,1:])**2+(yp0[0,:-1]-yp0[0,1:])**2)
        xmn,xmx,ymn,ymx = self.domain
        for it in range(1,nt-2):
            #print it
            x0,y0,sens = xp0[it,mdle]*1,yp0[it,mdle]*1,1
            xp1[it,mdle],yp1[it,mdle] = x0*1,y0*1
            vloc[it,mdle] = (cua[it+1,mdle]-cua[it-1,mdle])/(tp0[it+1,mdle]-tp0[it-1,mdle])
            for ip in range(mdle+1,np):
                wdth[it,ip-1],(x0,y0),v = self.distAtPoint(xp1,yp1,wdth,vloc,x0,y0,it,ip,sens)
                xp1[it,ip],yp1[it,ip],vloc[it,ip] = x0*1,y0*1,v*1
            x0,y0,sens = xp0[it,mdle]*1,yp0[it,mdle]*1,-1
            for ip in range(mdle-1,-1,-1):
                wdth[it,ip],(x0,y0),v = self.distAtPoint(xp1,yp1,wdth,vloc,x0,y0,it,ip,sens)
                xp1[it,ip],yp1[it,ip],vloc[it,ip] = x0*1,y0*1,v*1
                
        #wd1=sqrt((self.xp1[:,1:]-self.xp1[:,:-1])**2+(self.yp1[:,1:]-self.yp1[:,:-1])**2)
        xc1,yc1 = (xp1[:,1:]+xp1[:,:-1])/2,(yp1[:,1:]+yp1[:,:-1])/2
        # find the local velocity
        vloc = (vloc[:,1:]+vloc[:,:-1])/2
    #if 3>2:
        # calculate the concentration (considering the deviaiton from the flow angle)
        Fo = wdth[0,:]
        angl1=arctan((yp0[0,np/2+1]-yp0[0,np/2])/(xp0[0,np/2+1]-xp0[0,np/2]))
        angl2=arctan((yp0[1,np/2]-yp0[0,np/2])/(xp0[1,np/2]-xp0[0,np/2]))
        Fo=Fo*vloc[0,:]*sin(abs(angl1-angl2));
        # CY global
        F=Fo*(wdth*0.+1.);Cy=array(F/wdth/vloc);sumF = sum(Fo)
        Cy,wdth1 = self.smoo(Cy),wdth*1
        for i in range(5,nt-1):
            Cy[i,:]=minimum(Cy[i,:],Cy[i-1,:])
        '''
        for i in range(3,nt):
            n = 0
            while max(Cy[i,:])>max(Cy[i-1,:]):
                print n
                a = Cy[i,:]*vloc[i,:]
                Cy[i,1:-1] = (a[2:]+a[1:-1]+a[:-2])/3/vloc[i,1:-1]
                #Cy[i,1:-1] = (Cy[i,2:]+Cy[i,1:-1]+Cy[i,:-2])/3
                sumFi = sum(Cy[i,:]*vloc[i,:]*wdth[i,:])
                Cy[i,:] = Cy[i,:]*sumF/sumFi
                n+=1
        '''
        '''
        for i in range(5,nt-1):
            dff = Cy[i,:]-Cy[i-1,:]
            if amax(dff)>0: #print i
                n = 0; #print i
                while (n<5) and (amax(dff))>0:
                    jmx = where(dff==max(dff))[0][0]
                    Cy[i,jmx] = Cy[i-1,jmx]
                    f0 = dff[jmx]*vloc[i,jmx]*wdth[i,jmx]
                    put2next(f0,dff,Cy,vloc,wdth,jmx,np,i)
                    n +=1
        '''
        print('inst 235',amax(Cy))
        return [xp0,yp0,tp0,cua,cub,vloc,wdth,Cy]
        
    def getLcoefs(self,poly):
        '''finds the list of coefficients of the equation of each line 
        segment in a polygon'''
        x,y = list(zip(*poly)); #a,b = str(x),str(y)
        n = len(x)
        lcoefs,a=zeros((2,n-1)),zeros((2,2))
        for i in range(n-1): # lcoefs are the coefficient of the lines ax+by=1
            #print i
            a[:,0]=x[i:i+2];a[:,1]=y[i:i+2];#print i,a
            if sum(a[1,:]-a[0,:])==0:
                lcoefs[:,i] = 1
            else :
                lcoefs[:,i] = solve(a,ones((2,1))[:,0])
        return lcoefs
    
    def dintersect(self,x0,y0,it,ip):
        # it is the first point of the segment
        a,b = self.lcoefs[ip][0,it],self.lcoefs[ip][1,it]
        d1 = abs(a*x0+b*y0-1)/sqrt(a**2+b**2) # distances from pt to segment
        xi = (-y0+b/a*x0+1/b)/(b/a+a/b); yi = (1-a*xi)/b # intersection pt
        t = (xi>=self.xp0[it,ip])&(xi<=self.xp0[it+1,ip])&(yi>=self.yp0[it,ip])&(yi<=self.yp0[it+1,ip])
        return d1,(xi,yi),t
                    
    def distAtPoint(self,xp1,yp1,wdth,vloc,x0,y0,it,ip,sens):
        '''returns the distance from the x0,y0 on the ip line, reminf that size of 
        wdth is smaller than the others'''
        dst = sqrt((x0-self.xp0[:,ip])**2+(y0-self.yp0[:,ip])**2)
        idx = argsort(dst)
        it1 = min(idx[0]-1,self.nt-3)
        if self.flg[it1,ip]: # the closest point is stagnant
            return wdth[it-1,ip-1],(2*x0-xp1[it,ip-sens],2*y0-yp1[it,ip-sens]),vloc[it-1,ip]
        d1,pi,t = self.dintersect(x0,y0,it1,ip)
        if t :
            v = (self.cua[it1+1,ip]-self.cua[it1,ip])/(self.tp0[it1+1,ip]-self.tp0[it1,ip])
            return d1,pi,v
        else : 
            d1,pi,t = self.dintersect(x0,y0,it1+1,ip)
            v = (self.cua[it1+2,ip]-self.cua[it1+1,ip])/(self.tp0[it1+2,ip]-self.tp0[it1+1,ip])
            return d1,pi,v

    def calcConc_old(self,xp0,yp0,tp0,cua,cub):
        """ on part des positions des particules limitant les tubes de flux et on calcule
        la concentration le long de chauqe tube pour un temps infini. On va produire de nouvelles
        positions (XP1, YP1) qui sont en fait le centre des cellules irregulieres cree
        par CalculP. en chaque point on aura une concentration
        """
        def put2next(f0,dff,Cy,vloc,dyp,j0,np,i):
            n=0
            while n<np:
                j = j0+n
                if j<np-1 :
                    if (dff[j]<0) and (f0>0):
                        dC = min(-dff[j],f0/vloc[i,j]/dyp[i,j])
                        Cy[i,j] += dC
                        f0 -= dC*vloc[i,j]*dyp[i,j]
                j = j0-n
                if j>=0 :
                    if (dff[j]<0) and (f0>0):
                        dC = min(-dff[j],f0/vloc[i,j]/dyp[i,j])
                        Cy[i,j] += dC
                        f0 -= dC*vloc[i,j]*dyp[i,j]
                n += 1
        #calcul du rapport entre cu2 (chemin reel) et cu (chemin sans disp lat)                    
    #if 3>2:
        [nt,np]=shape(cua);r=cub[-1,:]/cua[-1,:]; 
        x0,x1,y0,y1 = self.domain
        dx,nx = self.dx,self.ncol
        dy,ny = self.dy,self.nrow
        poro = 0.25
        epais = 10
        dc0=min(dx,dy); #sqrt(dx**2.+dy**2.)*1.01
        ntmx=int(max(cua[-1,:]/dc0))+1; # nbre d'intervalle
        dc=cua[-1,:]/ntmx;#dc=clip(dc,dc0*.98,dc0*1.02)
        tet0=[0.];ind=array((ntmx,nt,np));
        #[xp1,yp1,tp1,cu1,cu2]= self.rflow.ptsLigne(ind,tet0,dc,xp0,yp0,tp0,cua,cub)
        xp1,yp1,tp1,cu1,cu2= xp0,yp0,tp0,cua,cub
        # on va chercher la largeur d'interbandes avec cercles inscrits dans les segments
        # calcul du point d'intersection x0 des segments
        xp1=array(xp1);yp1=array(yp1);tp1=array(tp1);cu1=array(cu1); ## pb de compatiblite ancien et new Array??
        a1=yp1[1:,:]*xp1[:-1,:]-yp1[:-1,:]*xp1[1:,:]
        xa1=xp1[1:,:]-xp1[:-1,:];ya1=yp1[1:,:]-yp1[:-1,:];putmask(xa1,xa1==0,1)
        b1=xa1[:,1:]*a1[:,:-1]-xa1[:,:-1]*a1[:,1:]
        b2=xa1[:,:-1]*ya1[:,1:]-xa1[:,1:]*ya1[:,:-1];putmask(b2,b2==0,1)
        x0=-b1/b2
        y0=(yp1[1:,1:]*(x0-xp1[:-1,1:])-yp1[:-1,1:]*(x0-xp1[1:,1:]))/xa1[:,1:]
         # distances entre pts et x0 et calcul des points finaux
        d1=sqrt((xp1[:-1,:-1]-x0)**2+(yp1[:-1,:-1]-y0)**2)
        d2=sqrt((xp1[1:,:-1] -x0)**2+(yp1[1:,:-1] -y0)**2)
        d3=sqrt((xp1[:-1,1:] -x0)**2+(yp1[:-1,1:] -y0)**2)
        d4=sqrt((xp1[1:,1:]  -x0)**2+(yp1[1:,1:]  -y0)**2)
        dmx=maximum(d1,d3);dmn=minimum(d2,d4);d=dmn/2+dmx/2
        xc1=x0+(xp1[:-1,:-1]-x0)*d/d1;yc1=y0+(yp1[:-1,:-1]-y0)*d/d1
        xc2=x0+(xp1[:-1,1:] -x0)*d/d3;yc2=y0+(yp1[:-1,1:] -y0)*d/d3
        # xp2, tp2 calcule au point centre de la bande size:np-1,nt-1
        xp2=xc1/2+xc2/2;d1b=d2-d1
        yp2=yc1/2+yc2/2;d3b=d4-d3
        putmask(xp2,xp2==0,1);putmask(d1b,d1b==0,1);putmask(d3b,d3b==0,1)
        a1=(d-d1)/d1b;a3=(d-d3)/d3b
        tp2=tp1[:-1,:-1]*(1-a1)+ tp1[1:,:-1]*a1+ tp1[:-1,1:]*(1-a3)+ tp1[1:,1:]*a3;tp2=tp2/2
        cu3=cu1[:-1,:-1]*(1-a1)+ cu1[1:,:-1]*a1+ cu1[:-1,1:]*(1-a3)+ cu1[1:,1:]*a3;cu3=cu3/2
        a=arange(np);a1=compress(a!=round(np/2),a)
        a2=take(tp1[1:,:],a1,axis=1);putmask(tp2,a2==0,0);putmask(cu3,a2==0,0)
        # largeur de la bande
        dyp=sqrt((xc1-xc2)**2+(yc1-yc2)**2)
        # calcul vloc
        [l,c]=shape(xp2);a2=a2[1:,:]
        tp3=tp2[1:l,:]-tp2[:l-1,:];putmask(tp3,tp3<=0,1.)
        vloc=(cu3[1:l,:]-cu3[:l-1,:])/tp3
        vloc=concatenate([vloc,vloc[l-2:l,:]],axis=0);#vloc=vloc/vloc[1,:];
        Fo=sqrt((xp0[0,:-1]-xp0[0,1:])**2+(yp0[0,:-1]-yp0[0,1:])**2)
        angl1=arctan((yp0[0,np/2+1]-yp0[0,np/2])/(xp0[0,np/2+1]-xp0[0,np/2]))
        angl2=arctan((yp0[1,np/2]-yp0[0,np/2])/(xp0[1,np/2]-xp0[0,np/2]))
        Fo=Fo*vloc[1,:] *sin(abs(angl1-angl2));
        # CY global
        F=Fo*(dyp*0.+1.);Cy=F/dyp/vloc
        nt,np=shape(xp1)
        #print F[i,16:18]
        #putmask(Cy,(vloc<=0.)|(dyp==0.)|(Cy<0.),0.)#|(indx2[1:,:]==0)|(Cy>1.)
        #print xp1[-1,:],yp1[-1,:],tp1[:,-1],cu1[-1,:],cu2[-1,:],vloc[-1,:],dyp[-1,:],Cy[-1,:];
        # technique brutale
        #for i in range(5,ntmx-1):
        #    Cy[i,:]=minimum(Cy[i,:],Cy[i-1,:])
        # flux distribution
        for i in range(5,ntmx-1):
            dff = Cy[i,:]-Cy[i-1,:]
            if amax(dff)>0: #print i
                n = 0; #print i
                while (n<5) and (amax(dff))>0:
                    jmx = where(dff==max(dff))[0][0]
                    Cy[i,jmx] = Cy[i-1,jmx]
                    f0 = dff[jmx]*vloc[i,jmx]*dyp[i,jmx]
                    put2next(f0,dff,Cy,vloc,dyp,jmx,np,i)
                    n +=1
        Cy=concatenate((ones((1,np-1)),Cy),axis=0);Cy=clip(Cy,0.,1.)
        vloc=concatenate((vloc[:1,:],vloc),axis=0)
        dyp=concatenate((dyp[:1,:],dyp),axis=0)
        #print 'vdyp',vloc[:,np/2-3:np/2+3],Cy[:,np/2-3:np/2+3]
        for i in range(5): Cy=self.smoo(Cy)
        return [xp1,yp1,tp1,cu1,cu2,vloc,dyp,Cy]

    def plume2reg(self):
        '''creates a regular matrix from the C plume otained in calcTube'''
        M  = zeros((self.nrow,self.ncol))
        ix = around(ravel(self.xp1/self.dx)).astype('int');ix = list(clip(ix,0,self.ncol-1))
        iy = around(ravel(self.yp1/self.dy)).astype('int');iy = list(clip(iy,0,self.nrow-1))
        M[iy,ix]=ravel(self.C)
        return M
