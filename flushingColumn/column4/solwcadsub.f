c
c *** solwcadsub2
c
c solwcadsub2 is a set of routines to calculate the saturation content of
c water and carbon dioxide in silicate liquids, and the composition
c of the co-existing gas phase, for any given silicate liquid composition
c in terms of 10 major oxides, pressure, temperature, and total H2O and CO2 
c in the system.
c
c author: Paolo Papale
c date  of last version: December 2005
c
c Primary reference:
c
c Papale, P., Moretti, R., Barbato, D. (2006) The compositional dependence
c of the saturation surface of H2O+CO2 fluids in silicate melts. Chem.
c Geol. 229: 78-95.
c
c Additional references:
c
c Papale, P. (1999) Modeling of the solubility of a two-component
c H2O + CO2 fluid in silicate liquids. Am. Mineral., 84: 477-492.
c
c Papale, P. (1997) Modeling of the solubility of a one-component 
c H2O or CO2 fluid in silicate liquids. Contrib. Mineral. Petrol.,
c 126: 237-251.
c
c -------------------------------------------------------------------------
c
c It makes use of newt.f from Numerical Recipes in FORTRAN, 2nd Edition,
c Press et al., Cambridge University Press, 1992
c
c -------------------------------------------------------------------------
c
      subroutine solwcad(p,t,ox,xy,igas)
c
      implicit real*8 (a-h,o-z)
!f2py intent(in) :: p
!f2py intent(in) :: T
!f2py intent(in) :: ox
!f2py intent(out) :: xy
c
      parameter (nox=12,n=2,neq=2,tol=1.d-1,tolscl=.9,tresho=5.d-5,
     *           tolest=1.d-12)
      common /xgas/ par(10),to(2),rgas,vr,imode
      common /funcvc/ ctrm(2),ctrm2(2),xtot(2),oxtot(2),
     *                rtrm(2),y(2),xtrm(2),pc,tc,whc,xwo,poco2,igs
      common /compc/ pm(nox),f(nox),fmol(nox)
      common /solwcadc/ xym(4),vref(2),vexc(2),gam(2),act(2),phig(2),
     *                  phi(2),wtfrex
      common /solwcadd/ estxy(2)
      common /scale/pref,xscale,iscale,isc
      common /xtt/ xt(2)
      common /modpar/ wo(nox,nox),wij(20,2),fol(2)
      common /check/ xdissw,xdissc,resid2,xgxg,phic,frefc,
     *               xwoc,ctrmc
      common /twij/ ttrm
      common /fff/ ffref1,ffref2
      dimension wo83(nox,nox)
      dimension ox(nox),ox2(nox),x(nox),x2(nox),xy(4),
     *          h2oco2(n),
     *          fref(2),
     *          xgg(neq),xgg1(neq),fvec(neq)
      dimension f1(nox)
      dimension oxnew(nox),xnew(nox)
      dimension covp(20,2),covpar(10),covfol(2)
      dimension wijmax(20,2),wijmin(20,2),parmax(10),parmin(10),
     *          folmax(2),folmin(2)
      logical check
      logical plow,phigh
c
      data pm /18.02,44.01,60.085,79.899,101.96,159.69,71.846,
     *         70.937,40.304,56.079,61.979,94.195/
      data f1 /1.,1.,.25,.25,.375,.375,.25,.25,.25,.25,.375,.375/
c
c Model parameters ---------------------
c
      data wo83 /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *       0.,0.,0.,
     *     0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *     0.,-122861.0680,-328708.4288,11037.09912,-40292.50576,
     *     23118.10624,-126999.4624,-268060.9304,-308604.7272,
     *     -366503.3376,0.,0.,0.,0.,-281791.1448,-28542.49488,
     *     -19223.76456,-8548.7488,53026.3424,-428617.328,-422893.616,
     *     -170291.7288,0.,0.,0.,0.,0.,5189.49888,-249067.6624,
     *     -8023.866,-203655.3632,-411824.0072,-567413.16,-733563.984,
     *     0.,0.,0.,0.,0.,0.,18930.34064,887.828064,-5343.09352,
     *     6358.88504,-15553.51792,1187.109584,0.,0.,0.,0.,0.,0.,0.,
     *     -2942.77456,-242361.5472,-248343.412,-154666.5808,
     *     -353880.628,0.,0.,0.,0.,0.,0.,0.,0.,-11757.4584,
     *     2925.130632,3264.1476,-254.0696344,0.,0.,0.,0.,0.,0.,
     *     0.,0.,0.,-330220.108,-387486.0976,-188961.5736,0.,0.,0.,
     *     0.,0.,0.,0.,0.,0.,0.,-262671.1016,-116767.072,0.,0.,0.,0.,
     *     0.,0.,0.,0.,0.,0.,0.,-75854.6648,0.,0.,0.,0.,0.,0.,0.,0.,
     *     0.,0.,0.,0./
c
      data to /.1000000000000000d+04,.1400000000000000d+04/
c
      data wij /-0.3409307328E+05,0.0000000000E+00,-0.1891165069E+06,
     *          0.1359353090E+06,-0.1957509655E+06,0.0000000000E+00,
     *          -0.8641805096E+05,-0.2099966494E+06,-0.3222531515E+06,
     *          -0.3497975003E+06,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *          -0.5996216205E+05,0.0000000000E+00,-0.5909571939E+06,
     *          0.4469622832E+07,0.2166575511E+05,0.0000000000E+00,
     *          0.5286613028E+05,-0.3287921030E+06,0.1400344548E+06,
     *          0.3090699699E+06,0.6049208301E+04,0.0000000000E+00,
     *          0.4139537240E+05,-0.5293012049E+06,0.1213717202E+04,
     *          0.0000000000E+00,-0.1344620202E+05,0.1278883196E+05,
     *          -0.3521319385E+05,-0.5800953852E+05/
c
      data fol /0.1979347855E+02,0.2348447513E+02/
      data par /0.2418946564E-05,0.1344868775E-07,0.,0.,
     *          0.2106552365E-14,0.,0.,0.,0.,0./
      data covp /0.6316242175E+03,0.0000000000E+00,0.4801528874E+04,
     *           0.1267547652E+05,0.6130377278E+04,0.0000000000E+00,
     *           0.6097020738E+04,0.3500308057E+04,0.4581663847E+04,
     *           0.6257676676E+04,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *           0.1065551645E+05,0.0000000000E+00,0.1887635505E+06,
     *           0.3898097404E+06,0.1746504717E+06,0.0000000000E+00,
     *           0.1001263771E+06,0.7567161446E+05,0.2041019346E+06,
     *           0.1495171552E+06,0.1162983976E+04,0.0000000000E+00,
     *           0.1936693988E+05,0.4662886265E+05,0.1820384707E+05,
     *           0.0000000000E+00,0.1051294566E+05,0.8091964307E+04,
     *           0.2114478691E+05,0.1594015370E+05/
      data covfol /0.1001737547E-01,0.1095849897E+00/
      data covpar /0.1570083361E-05,0.6432481672E-09,0.,0.,
     *             0.1399168896E-14,0.,0.,0.,0.,0./
c
      m=0
      pc=p
      poco2=0.1013250000E+06
      ittrm=0
      wijmix=0.
c
      do i=1,20
        wijmax(i,1)=wij(i,1)+covp(i,1)
        wijmin(i,1)=wij(i,1)-covp(i,1)
        wijmax(i,2)=wij(i,2)+covp(i,2)
        wijmin(i,2)=wij(i,2)-covp(i,2)
      enddo
      folmax(1)=fol(1)+covfol(1)
      folmin(1)=fol(1)-covfol(1)
      folmax(2)=fol(2)+covfol(2)
      folmin(2)=fol(2)-covfol(2)
      do i=1,10
        parmax(i)=par(i)+covpar(i)
        parmin(i)=par(i)-covpar(i)
      enddo
c
      if(ittrm.eq.1) then
        ttrm=t
       else
        ttrm=1.
      endif
c
      do i=1,nox
        f(i)=f1(i)
      enddo
c
c Data from Ghiorso et al. 1983
      do i=3,12
        do j=i+1,12
          wo(j,i)=wo83(j,i)
          wo(i,j)=wo(j,i)
        enddo
      enddo
c
      if(igas.eq.0) goto 9876
      plow=.false.
      phigh=.false.
      if(p.ge.1.d0) plow=.true.
      if(p.le.1.d12) phigh=.true.
      if(plow.and.phigh) then
        goto 9876
       else
        write(*,*) '**error in calling solwcad'
        write(*,*) 'p = ',p,' out of range'
        stop
      endif
c
 9876 continue
      igs=igas
      mm=0
      imode=2
      ligas=0
      ligs=0
      iscale=0
      poh2o=0.
c      poco2=1.01325d5
      rgas=8.31451d0
      sg=4.1d-6
      expg=.5d0
      eps=1.d-12
      redfac=.8d0
      tc=t
      do i=1,2
        h2oco2(i)=0.
      enddo
      if(igas.ne.3) then
        do i=1,4
          xy(i)=0.
          if(i.le.2) then
            vref(i)=0.
            vexc(i)=0.
            gam(i)=0.
            act(i)=0.
            phig(i)=0.
          endif
        enddo
      endif
      do i=3,nox
        do j=i+1,nox
          wo(i,j)=wo(j,i)
        enddo
      enddo
c
c re-express composition -------------------------
c
   29 call comp(ox,x)
check---------------------------
      do i=1,12
        oxnew(i)=ox(i)
      enddo
      oxnew(1)=xdissw
      oxnew(2)=xdissc
      call comp(oxnew,xnew)
check---------------------------
      do i=1,2
        xtot(i)=x(i)
        oxtot(i)=ox(i)
        xt(i)=xtot(i)
      enddo
      if(igas.eq.0) then
              xy(1)=ox(1)
              xy(2)=ox(2)
              xym(1)=x(1)
              xym(2)=x(2)
      endif
      oxtott=ox(1)+ox(2)
      if(igas.ne.3) oxtott=ox(igas)
c
c compositional terms -------------------------
c
      pln=dlog(p/poco2)
      xwo=0.
      gammxh=-1.d6
      gammxc=-1.d6
      gammnh=1.d6
      gammnc=1.d6
      kgam=0
      kgam1=0
 7934 continue
c
 7935 continue
      do i=1,2
        ctrm(i)=0.
        ctrm2(i)=0.
      enddo
      do i=3,nox
        ctrm(1)=ctrm(1)+x(i)*wij(i-2,1)
        ctrm(2)=ctrm(2)+x(i)*(wij(i-2,2)+pln*wij(i+8,2)*ttrm)
c        write(*,*) i,x(i),pln,wij(i-2,2),wij(i+8,2)
        ctrm2(2)=ctrm2(2)+x(i)*wij(i+8,2)
        if(igas.eq.0) ctrm2(1)=ctrm2(1)+x(i)*wij(i-2,2)
        if(i.eq.nox) goto 11
        do j=i+1,nox
          xwo=xwo+x(i)*x(j)*wo(i,j)
        enddo
   11   continue
      enddo
c
      do i=1,2
        xtrm(i)=(ctrm(i)-xwo)/(rgas*t)
      enddo
      if(igas.eq.0) then
              xtrm(1)=(1.-x(1))*ctrm(1)-
     *                x(2)*ctrm2(1)-
     *                (1.-x(1)-x(2))*xwo
              xtrm(1)=xtrm(1)*(1.-x(1)-x(2))/(rgas*t)
              xtrm(2)=-x(1)*ctrm(1)+
     *                (1.-x(2))*ctrm2(1)-
     *                (1.-x(1)-x(2))*xwo
              xtrm(2)=xtrm(2)*(1.-x(1)-x(2))/(rgas*t)
      endif
c
c reference terms -------------------------
c
      if(igas.eq.0) then
              xgg(1)=.5
              pref=(ox(1)/sg)**(1./expg)
              if(ox(1).lt..01 .and. ox(2)/ox(1).ge..01) pref=pref*10.
c              pref=p
              xgg(2)=1.
              goto 81
      endif
      if(igas.ne.3) goto 21
      do ig=1,2
        call fpt(p,t,ycd,ig,fref,phi)
        rtrm(ig)=fol(ig)+fref(ig)/(rgas*t)
        vref(ig)=vr
      enddo
      goto 22
   21 ycd=0.
      call fpt(p,t,ycd,igas,fref,phi)
      rtrm(igas)=fol(igas)+fref(igas)/(rgas*t)
      vref(igas)=vr
   22 continue
c
c decide if one or two gases -----------------------
c
   30 kigas=0
      if(igas.eq.3) then
        if(m.eq.0) then
          igas=1
          kigas=1
         else
          goto 101
        endif
      endif
c
c
c -------- one gas (or two gases and setting initial guess) --------
c ----------------------------------------------------------------------
check---------------------------
      if(igas.ne.3 .and. kigas.ne.1) then
        bterm=rgas*t*dlog(phi(igas)*p/xnew(igas))
        bterm=bterm+(1.-xnew(igas))**2*xwo
        aterm=(1.-xnew(igas))**2*ctrm(igas)+rgas*t*fol(igas)
        aterm=aterm+fref(igas)
        resid2=bterm-aterm
        xgxg=xnew(igas)
        phic=phi(igas)
        frefc=fref(igas)
        xwoc=xwo
        ctrmc=ctrm(igas)
      endif
check---------------------------
c
c set initial guess for one gas
c
      if(m.eq.1) then
c        write(*,*) 'm=1'
        estxy(1)=xtot(1)
        estxy(2)=xtot(2)
        goto 176
      endif
      ox2(1)=sg*p**expg
      ox2(2)=0.
      do i=3,nox
        ox2(i)=ox(i)
      enddo
      call comp(ox2,x2)
      estxy(1)=x2(1)
      estxy(2)=estxy(1)*.1
   31 continue
c
c start iteration -----------------------------------
c
      ah=.8*estxy(igas)
      bh=1.2*estxy(igas)
   39 if(bh.ge.1.) bh=1.
      k=0
c
   40 call nb02ad(k,ah,bh,eps,xg,yh)
      goto(61,62,63),k
c
c set yh(xg) --------------
c
   61 continue
      yh=dlog(phi(igas)*p/xg)-(1.-xg)**2*xtrm(igas)-rtrm(igas)
      goto 40
c
c change boundaries -----------
c
   63 ah=.5*ah
      bh=1.5*bh
      if(ah.lt.1.d-15 .and. bh.gt.1.) then
        write(*,*) '*** error from solwcad2'
        write(*,*) 'at P = ',p,'  and T = ',t
        write(*,*) 'the solution cannot be bracketed'
        return
       else
        goto 39
      endif
c
c xg = x(igas) has been found ------------------
c
   62 do i=1,4
        if(i.le.2) x(i)=0.
        xy(i)=0.
        xym(i)=0.
      enddo
      if(kigas.ne.0) goto 76
      if(xg.ge.xtot(igas)) then
        ligs=1
        x(igas)=xtot(igas)
        xym(igas)=xtot(igas)
        xym(igas+2)=1.
       else
        x(igas)=xg
        xym(igas)=xg
        xym(igas+2)=1.
      endif
      write(*,*) 'xym ',xym(3),xym(4)
      if(kigas.eq.0) goto 70
c
c set initial estimate for two gases ----------------
c
   76 estxy(1)=xg
      estxy(2)=estxy(1)*.1
  176 continue
      do i=1,2
        if(estxy(i).ge.xtot(i)) estxy(i)=xtot(i)-1.d-6
      enddo
      igas=3
      m=1
      goto 30
c
c ------------------------------ two gases -----------------------------
c ----------------------------------------------------------------------
c
  101 continue
c
      pc=p
      tc=t
c
c scaling unknowns ---------
c
      kspur=0
   91 do i=1,neq
        inot=2
        if(i.eq.inot) inot=1
        if(xtot(i).lt.tresho) then
          igas=inot
          ox(i)=0.
          oxtott=ox(igas)
          m=0
          ligas=1
          goto 29 
        endif
        if(estxy(i)/estxy(inot).lt.tolscl) then
          iscale=1
          isc=i
          xscale=estxy(i)/estxy(inot)
          xt(i)=xt(i)/xscale
          estxy(i)=estxy(inot)
        endif
      enddo
c
      if(kspur.ne.0) goto 81
      do i=1,neq
        xgg(i)=estxy(i)
        xgg1(i)=xgg(i)
      enddo
c
c start numerical procedure -------------------------
c
   81 continue
c
      call newt(xgg,neq,check)
c
      if(igas.eq.0) goto 5546
      call funcv(neq,xgg,fvec)
c
      if(check) then
        write(*,*) 'convergency problems in newt'
      endif
      do i=1,2
        if(dabs(fvec(i)).gt.tol) then
          kspur=kspur+1
          m=0
c undersaturated conditions ------------
          do j=1,neq
            xgg(j)=xgg1(j)*.8/kspur
          enddo
          if(kspur.ge.3) then
            do ispur=1,2
              xgg(ispur)=xtot(ispur)
            enddo
            goto 82
          endif
          goto 91
        endif
   82   if(iscale.eq.1 .and. i.eq.isc) xgg(i)=xgg(i)*xscale
        x(i)=xgg(i)
        xym(i)=xgg(i)
      enddo
      xym(3)=y(1)
      xym(4)=1.-xym(3)
      goto 70
c
 5546 continue
      xym(3)=xgg(1)
      xym(4)=1.-xym(3)
      p=xgg(2)*pref
      pln=dlog(p/poco2)
c
   70 continue
c
c ----------------------------------------------------------------------
c set properties and xy(1 -> 4) ----------------------------------------
c ----------------------------------------------------------------------
c
      xvolo=1.-x(1)-x(2)
c H2O
      if(igas.ne.2) then
        vexc(1)=(-x(2)*xvolo*ctrm2(2)+(1.-x(1))*x(2)*h2oco2(2))/p
        gam(1)=(1.-x(1))*(xvolo*ctrm(1)+x(2)*(h2oco2(1)+pln*h2oco2(2)))
     *         -x(2)*xvolo*ctrm(2)-xvolo**2*xwo
        gam(1)=dexp(gam(1)/(rgas*t))
        act(1)=gam(1)*x(1)
        phig(1)=phi(1)
      endif
c
c CO2
      if(igas.ne.1) then
        vexc(2)=(1.-x(2))*(xvolo*ctrm2(2)+x(1)*h2oco2(2))/p
        gam(2)=(1.-x(2))*(xvolo*ctrm(2)+x(1)*(h2oco2(1)+pln*h2oco2(2)))
     *         -x(1)*xvolo*ctrm(1)-xvolo**2*xwo
        gam(2)=dexp(gam(2)/(rgas*t))
        act(2)=gam(2)*x(2)
        phig(2)=phi(2)
      endif
c
      do i=1,2
        estxy(i)=x(i)
      enddo
c
      xyws=xym(3)*pm(1)+xym(4)*pm(2)
      do i=1,nox
        if(i.ge.3) x(i)=x(i)*xvolo
        fmol(i)=x(i)/f(i)
      enddo
      fmol(3)=fmol(3)+.5*(fmol(7)+fmol(8)+fmol(9)+fmol(10))
     *       +fmol(11)+fmol(12)
      oxs=0.
      do i=1,nox
        oxn=fmol(i)*pm(i)
        oxs=oxs+oxn
        if(i.le.2) xy(i)=oxn
      enddo
c
      xy(1)=xy(1)/oxs
      xy(2)=xy(2)/oxs
      xy(3)=xym(3)*pm(1)/xyws
      xy(4)=1.-xy(3)
      if(igas.eq.0) goto 3435
c
c check mass conservation --------------
c
      wtfrex=(oxtott-xy(1)-xy(2))/(1.-xy(1)-xy(2))
c
c undersaturated conditions ------------------------
c
      if(wtfrex.lt.1.d-8) then
        wtfrex=0.
        if(igas.eq.3) then
            do i=1,4
              if(i.le.2) then
                xy(i)=oxtot(i)
               else
                xy(i)=0.
                xym(i)=0.
              endif
            enddo
        endif
      endif
c
 3435 continue
      close(14)
      close(16)
      if(ligas.eq.1) igas=3
c
c normalize volatile-free composition
c
      oxsn=0.
      do i=3,12
        oxsn=oxsn+ox(i)
      enddo
      do i=3,12
        ox(i)=ox(i)/oxsn
      enddo
c
      return
      end
c
c ***********************************************************************
c
      subroutine fpt(p,t,ycd,igas,fref,phi)
c
c returns:
c          fref = liquid components reference terms
c          phi = gas components fugacities
c
      implicit real*8 (a-h,o-z)
c
      common /xgas/ par(10),to(2),rgas,vr,imode
      common /vm1/ c(2,2),d(2,2),e(2,2),b(2),bm,cm,dm,em,ym
      common /qa05dd/ ipr
      common /fintc/ pc,tc,kint
      common /fluid2/vfluid
      common /gasvol/vfld,vmolh
      common /ff/frefcm
c
      dimension bd(10),fref(2),phi(2),x(2)
      dimension xsfl(13),fgsfl(13),actsfl(13)
      external fint
c
      data bd /9.144d-6,3.685d-9,1.168d-11,-1.99d-15,1.22d-16,
     *         -1.945d-17,-1.58d-21,4.68d-24,1.144d-26,-3.96d-33/
c
c calculates phi(i), i = 1 -> 2 (two-component gas phase) ------------
c
      tp=t**1.5
      v=vm(p,t,ycd,igas)
c      vfluid=v
      vfld=v
      vr=0.
      if(igas.eq.1) then
        vr=bd(1)+bd(2)*t+bd(3)*t**2+bd(4)*t**3
     *     +p*(bd(5)+bd(6)*t+bd(7)*t**2)
     *     +p**2*(bd(8)+bd(9)*t)
     *     +p**3*bd(10)
        vmolh=vr
       else
        if(igas.eq.2)
     *  vr=par(1)+par(2)*t+par(3)*t**2+par(4)*t**3
     *     +p*(par(5)+par(6)*t+par(7)*t**2)
     *     +p**2*(par(8)+par(9)*t)
     *     +p**3*par(10)
c     *     +par(10)*p/t
      endif
      vln=dlog((v+bm)/v)
c
      if(igas.eq.3) then
        x(1)=1.-ycd
        x(2)=ycd
        do i=1,2
          phex=(4.*ym-3.*ym**2)/((1.-ym)**2)+b(i)/bm*(4.*ym-2.*ym**2)/
     *         ((1.-ym)**3)-(2.*c(i,i)*x(i)+2.*(1.-x(i))*c(1,2))/
     *         (rgas*tp*bm)*vln-cm*b(i)/(rgas*tp*bm*(v+bm))+
     *         cm*b(i)/(rgas*tp*bm**2)*vln-
     *         (2.*d(i,i)*x(i)+2.*(1.-x(i))*d(1,2)+dm)/(rgas*tp*bm*v)+
     *         (2.*d(i,i)*x(i)+2.*(1.-x(i))*d(1,2)+dm)/(rgas*tp*bm**2)*
     *         vln+b(i)*dm/(rgas*tp*v*bm*(v+bm))+
     *         2.*b(i)*dm/(rgas*tp*bm**2*(v+bm))-2.*b(i)*dm/
     *         (rgas*tp*bm**3)*
     *         vln-2.*(e(i,i)*x(i)+(1.-x(i))*e(1,2)+em)/
     *         (rgas*tp*2.*bm*v**2)+2.*(e(i,i)*x(i)+
     *         e(1,2)*(1.-x(i))+em)/
     *         (rgas*tp*bm**2*v)-2.*(e(i,i)*x(i)+e(1,2)*(1.-x(i))+em)/
     *         (rgas*tp*bm**3)*vln+em*b(i)/(rgas*tp*2.*bm*v**2*(v+bm))-
     *         3.*em*b(i)/(rgas*tp*2.*bm**2*v*(v+bm))+3.*em*b(i)/
     *         (rgas*tp*bm**4)*vln-3.*em*b(i)/(rgas*tp*bm**3*(v+bm))-
     *         dlog(p*v/(rgas*t))
          phi(i)=dexp(phex)
        enddo
        return
      endif
c      stop
c
c calculates phi(igas) (one-component gas phase) --------------
c
      phexp=(8.*ym-9.*ym**2+3.*ym**3)/((1.-ym)**3)-dlog(p*v/(rgas*t))-
     *      cm/(rgas*tp*(v+bm))-dm/(rgas*tp*v*(v+bm))-
     *      em/(rgas*tp*v**2*(v+bm))-cm/(rgas*tp*bm)*vln-
     *      dm/(rgas*tp*bm*v)+dm/(rgas*tp*bm**2)*vln-
     *      em/(rgas*tp*2.*bm*v**2)+em/(rgas*tp*bm**2*v)-
     *      em/(rgas*tp*bm**3)*vln
      phi(igas)=dexp(phexp)
c
c calculates fref(igas) ---------------------------------
c
      if(igas.eq.2) goto 5
c H2O
      t1=2.-t/to(1)
      t2=2.*t-to(1)
      t3=2.*t**2-to(1)**2
      fref(1)=p*(bd(1)*t1+bd(2)*t+bd(3)*t*t2+bd(4)*t*t3)+
     *        p**2*(bd(5)*t1+bd(6)*t+bd(7)*t*t2)*.5+
     *        p**3*(bd(8)*t1+bd(9)*t)/3.+
     *        p**4*bd(10)*t1*.25
      return
c
c CO2
    5 continue
      if(imode.eq.2) then
        t1=2.-t/to(2)
        t2=2.*t-to(2)
        t3=2.*t**2-to(2)**2
        fref(2)=p*(par(1)*t1+par(2)*t+par(3)*t*t2+par(4)*t*t3)+
     *          p**2*(par(5)*t1+par(6)*t+par(7)*t*t2)*.5+
     *          p**3*(par(8)*t1+par(9)*t)/3.+
     *          p**4*par(10)*t1*.25
        frefcm=fref(2)
        return
      endif
c
      ig=2
      po=1.01325d5
      aerr=1.d-12
      rerr=1.d-6
      level=2
      ipr=15
      pc=p
      tc=t
c
      kint=1
c      call qa05ad(frefp,fint,po,p,aerr,rerr,level,error,iflag)
c
      kint=2
c      call qa05ad(freft,fint,po,p,aerr,rerr,level,error,iflag)
c
      fref(2)=frefp-(t-to(2))/to(2)*freft
c
      return
      end
c
c ******************************************************************
c
      subroutine comp(ox,x)
c
c returns the mole fraction composition from the oxides wt distribution
c
      implicit real*8 (a-h,o-z)
c
      common /compc/ pm(12),f(12),fmol(12)
      dimension ox(12),x(12)
c
      xs=0.
      do i=3,12
        xs=xs+ox(i)
      enddo
      do i=1,12
        if(i.ge.3) ox(i)=ox(i)*(1.-ox(1)-ox(2))/xs
        fmol(i)=ox(i)/pm(i)
      enddo
      xsum=0.
      do i=1,12
        if(i.eq.3) fmol(i)=fmol(i)-.5*(fmol(7)+fmol(8)+fmol(9)+
     *                     fmol(10))-fmol(11)-fmol(12)
        x(i)=fmol(i)*f(i)
        xsum=xsum+x(i)
      enddo
      do i=1,12
        x(i)=x(i)/xsum
        if(i.ge.3) x(i)=x(i)/(1.-x(1)-x(2))
      enddo
c
      return
      end
c
c ******************************************************************
c
      real*8 function fint(x)
c
c called from qa05ad
c returns the function to be integrated
c
      implicit real*8 (a-h,o-z)
c
      common /fintc/ pc,tc,kint
c
      ig=2
      p=pc
      t=tc
c
c first integral (kint = 1) ------------------------
c
      v=vm(x,t,ycd,ig)
      if(kint.eq.1) then
        fint=v
        return
      endif
c
c second integral (kint = 2) -----------------------
c
      n=7
      ts=0.
      vs=0.
      tvs=0.
      t2s=0.
      dt=200.
c
      do i=1,n
        if(i.eq.1) then
          t1=t
        else
          t1=700.+dt*(i-1)
        endif
        v=vm(x,t1,ycd,ig)
        if(i.eq.1) v1=v
        ts=ts+t1
        vs=vs+v
        tvs=tvs+t1*v
        t2s=t2s+t1**2
      enddo
      dvdt=(n*tvs-ts*vs)/(n*t2s-ts**2)
      fint=v1-t*dvdt
c
      return
      end
c
c ******************************************************************
c
      real*8 function vm(p,t,ycd,igas)
c
      implicit real*8 (a-h,o-z)
      common /vm1/ c(2,2),d(2,2),e(2,2),b(2),bm,cm,dm,em,ym
      common /xvm/ ivm
      dimension x(2)
c
      r=8.31451d0
c
      if(igas.eq.2) goto 10
c
      b(1)=2.9d-5
      c(1,1)=(290.78-.30276*t+1.4774d-4*t**2)*1.d-1
      d(1,1)=(-8374.+19.437*t-8.148d-3*t**2)*1.d-7
      e(1,1)=(76600.-133.9*t+.1071*t**2)*1.d-13
      bm=b(1)
      cm=c(1,1)
      dm=d(1,1)
      em=e(1,1)
      if(igas.ne.3) goto 20
c
   10 continue
      b(2)=5.8d-5
      c(2,2)=(28.31+.10721*t-8.81d-6*t**2)*1.d-1
      d(2,2)=(9380.-8.53*t+1.189d-3*t**2)*1.d-7
      e(2,2)=(-368654.+715.9*t+.1534*t**2)*1.d-13
      bm=b(2)
      cm=c(2,2)
      dm=d(2,2)
      em=e(2,2)
      if(igas.ne.3) goto 20
c
      x(1)=1.-ycd
      x(2)=ycd
      bm=x(1)*b(1)+x(2)*b(2)
      if(c(1,1)*c(2,2).le.0.) then
            c(1,2)=0.
          else
            c(1,2)=dsqrt(c(1,1)*c(2,2))
      endif
      if(d(1,1)*d(2,2).le.0.) then
            d(1,2)=0.
          else
            d(1,2)=dsqrt(d(1,1)*d(2,2))
      endif
      if(e(1,1)*e(2,2).le.0.) then
            e(1,2)=0.
          else
            e(1,2)=dsqrt(e(1,1)*e(2,2))
      endif
      c(2,1)=c(1,2)
      d(2,1)=d(1,2)
      e(2,1)=e(1,2)
      cm=0.
      dm=0.
      em=0.
      do 15 i=1,2
      do 16 j=1,2
      cm=cm+x(i)*x(j)*c(i,j)
      dm=dm+x(i)*x(j)*d(i,j)
      em=em+x(i)*x(j)*e(i,j)
   16 continue
   15 continue
c
   20 continue
c
      ahfac=.1
      fac=.2
      eps=1.d-16
      kerr=0
      kerr2=0
      kerr3=0
      klim=20
c
c set boundaries ----------------------------------
c
      if(igas.eq.2) then
            ah=2.d-5
            goto 5
      endif
      ah=1.9901d-5-4.6125d-15*p+5.5523d-25*p**2
      if(igas.ne.1) then
            ah=ah+1.3539d-5*ycd
            if(t.le.1300.) ah=ah*.8
      endif
   98 ahst=ah
c
    5 bh=ah*(1.+fac)
c
c start iterating ----------------------------------
c
      k=0
c
    4 call nb02ad(k,ah,bh,eps,v,yh)
      goto (1,2,3),k
c
    1 a=cm+dm/v+em/(v**2)
      ym=.25*bm/v      
      yh=p-r*t*(1.+ym+ym**2-ym**3)/(v*(1.-ym)**3)+
     *   a/(dsqrt(t)*v*(v+bm))
      goto 4
c
c change boundaries ----------------------------------
c
    3 kerr=kerr+1
      if(kerr.eq.klim) then
            kerr=0
            kerr2=kerr2+1
c            write(19,*) 'kerr2 = ',kerr2
            if(kerr2.eq.100) then
                  kerr=0
                  kerr2=0
                  kerr3=kerr3+1
                  ah=ahst
                  ahfac=-ahfac
                  if(kerr3.eq.2) then
c                     write(*,*) '***igas = ',igas
                     ivm=1
                     return
                  endif
            endif
            fac=.1
            ah=ah*(1.+ahfac)
            goto 6
      endif
    6 fac=2.*fac
      goto 5
c
    2 vm=v
      yh=p-r*t*(1.+ym+ym**2-ym**3)/(v*(1.-ym)**3)+
     *   a/(dsqrt(t)*v*(v+bm))
c
      return
      end
c
c ***************************************************************
c
      subroutine funcv(n,x,f)
c
c called from newt
c restitutes the functions f(x) for any given value of x
c
      implicit real*8 (a-h,o-z)
c
      parameter(tin=1.d-16)
      parameter(nox=12)
      common /xgas/ par(10),to(2),rgas,vr,imode
      common /scale/pref,xscale,iscale,isc
      common /funcvc/ ctrm(2),ctrm2(2),xtot(2),oxtot(2),
     *                rtrm(2),y(2),xtrm(2),pc,tc,whc,xwo,poco2,igs
      common /modpar/ wo(nox,nox),wij(20,2),fol(2)
      common /solwcadc/ xym(4),vref(2),vexc(2),gam(2),act(2),phig(2),
     *                  phi(2),wtfrex
      common mm
      common /twij/ ttrm
      dimension x(n),f(n),fref(2),gtrm(2),xtrmx(2),x2(2)
c
      igas=igs
      if(igas.eq.0) goto 8878
c
      mm=0
      do i=1,2
        if(x(i).le.0.) x(i)=tin
      enddo
      x2(1)=x(1)
      x2(2)=x(2)
      do i=1,2
        if(iscale.eq.1 .and. i.eq.isc) x2(i)=x(i)*xscale
      enddo
      y(1)=(xtot(1)*(1.-x2(2))-x2(1)*(1.-xtot(2)))
     *     /(xtot(2)-x2(2)+xtot(1)-x2(1))
      if(y(1).le.tin) y(1)=tin
      if(y(1).ge.1.-tin) y(1)=1.-tin
      y(2)=1.-y(1)
      do i=1,2
        inot=2
        if(i.eq.inot) inot=1
        if(y(i).le.0.) then
          y(i)=tin
          y(inot)=1.-y(i)
        endif
        if(y(i).ge.1.) then
          y(i)=1.-tin
          y(inot)=1.-y(i)
        endif
      enddo
      ycd=y(2)
c      write(*,*) x(1),x(2),ycd
      p=pc
      t=tc
      ig=3
      mm=mm+1
      goto 8879
c
c ---- igas = 0 ----start
c DETERMINE P SATURATION AND GAS COMPOSITION FROM DISSOLVED H2O-CO2 PAIR
 8878 continue 
      t=tc
c
      ycd=1.-x(1)
      p=x(2)*pref
      xh=xym(1)
      xc=xym(2)
      do ig=1,2
        call fpt(p,t,ycd,ig,fref,phi)
        rtrm(ig)=fol(ig)+fref(ig)/(rgas*t)
      enddo
      ig=3
      call fpt(p,t,ycd,ig,fref,phi)
c
      trm1=(1.-xh-xc)*ttrm*dlog(p/poco2)*ctrm2(2)/(rgas*t)
      xtrmx(1)=xtrm(1)-xc*trm1
      xtrmx(2)=xtrm(2)+(1.-xc)*trm1
c
      f(1)=dlog(phi(1)*p*(1.-ycd)/xh)-xtrmx(1)-rtrm(1)
     *     +(1.-xh)*xc*whc
      f(2)=dlog(phi(2)*p*ycd/xc)-xtrmx(2)-rtrm(2)
     *     +(1.-xc)*xh*whc
c
      return
c
c ---- igas = 0 ----end
c
 8879 continue
c
      call fpt(p,t,ycd,ig,fref,phi)
c
      do i=1,2
        inot=2
        if(i.eq.inot) inot=1
        gtrm(i)=dlog(phi(i)*y(i)*p/x2(i))
        xtrm(i)=(1.-x2(1)-x2(2))*((1.-x2(i))*ctrm(i)-
     *          x2(inot)*ctrm(inot))
     *          -(1.-x2(1)-x2(2))**2*xwo+(1.-x2(i))*x2(inot)*whc
        xtrm(i)=xtrm(i)/(rgas*t)
        f(i)=gtrm(i)-xtrm(i)-rtrm(i)
        fugg=xtrm(i)+rtrm(i)+dlog(x2(i))
        fugg=dexp(fugg)
      enddo
c
      return
      end
c
c ***************************************************************
c
c
      subroutine nb02ad(k,az,bz,e2,x,y)
c
c solve non-linear equations by the secant method
c
      implicit real*8 (a-h,o-z)
      common /nb02b/a,b,ya,yb,it,j1,ytest,x1,y1
c
      if(k) 125,10,50
c
c     calculate y(x) at x=az.
   10 a = az
      b = bz
      x = a
      j1 = 1
      it = 1
      k = 1
      go to 20
c
c     branch depending on the value of j1.
   50 go to (60,120,170,170,170),j1
c
c     calculate y(x) at x=bz.
   60 ya = y
      x = b
      j1 = 2
      it=2
      go to 20
c
c     error return because no bracket
   70 k=3
      go to 20
c
c     get y for x = x2
   90 it = it+1
      x = x2
      go to 20
c
c     set the first bracket
  120 b = x
      yb = y
c     make test at 200 fail if ya or yb zero
  125 y=dmax1(dabs(ya),dabs(yb))*2.0
      j1 = 3
      if(ya*yb)126,200,70
  126 if(k.gt.0)go to 100
      k=1
      it=0
c
c     test wether bracket is small enough
  100 if((b-a).le.e2)go to 200
      go to (10,10,130,150,160),j1
c
c     calculate the next x by the secant method based on the bracket.
  130 if(dabs(ya) .le.dabs(yb))go to 140
      x1 = a
      y1 = ya
      x = b
      y = yb
      go to 150
  140 x1 = b
      y1 = yb
      x = a
      y = ya
c
c     use the secant method based on the function values y1 and y.
  150 u=y*(x-x1)/(y-y1)
  155 x2=x-u
      if(x2.eq.x)go to 110
      x1 = x
      y1 = y
      ytest =.50*dmin1(dabs(ya),dabs(yb))
c
c     check that x2 is inside the interval (a,b).
      if((x2-a)*(x2-b) .lt. 0.0)go to 90
c
c     calculate the next value of x by bisection.
  160 x2 = 0.50*(a+b)
      ytest = 0.0
c
c     check wether the maximum accuracy has been achieved.
      if((x2-a)*(x2-b))90,200,200
c
c     move away from fixed point to get close bracket
  110 if(u.eq.0.0)go to 160
      u=u+u
      go to 155
c
c     revise the bracket (a,b).
  170 if(y.eq.0.0)go to 195
      if(ya*y.gt.0.0)go to 180
      b = x
      yb = y
      go to 190
  180 a = x
      ya = y
c
c     use ytest to decide the method for the next value of x.
  190 j1=4
      if(dabs(y).gt.ytest)j1=5
      if(ytest .le.0.0)j1=3
      go to 100
c
c     y = 0 - set closest bracket
  195 if(dabs(x-a).lt.dabs(x-b))go to 196
      a=x
      ya=y
      go to 220
  196 b=x
      yb=y
      go to 220
  200 if(dabs(y).le.dabs(yb))go to 210
      x=b
      y=yb
  210 if(dabs(y).le.dabs(ya))go to 220
      x=a
      y=ya
  220 k=2
c
   20 continue
      return
      end
c
c ************************************************************************
c
      subroutine newt(x,n,check)
c
c from Numerical Recipes in FORTRAN
c Solver for systems of non-linear equations
c ---partially modified---
c
      implicit real*8 (a-h,o-z)
      parameter (neq=2,np=40,maxits=1000,tolf=1.d-12,
     *           tolmin=1.d-12,tolx=1.d-13,
     *stpmx=100.)
      parameter(tolstp=1.d-15)
      logical check
      common /newtv/ fvec(np),nn
      common /temp/ its
      common /step/ tolst
      common /xtt/ xtot(2)
CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      dimension x(n),fjac(np,np),g(np),p(np),xold(np),indx(neq)
      external fmin
      tolst=tolstp
      nn=n
      f=fmin(x)
      test=0.
      do 11 i=1,n
        if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
11    continue
      if(test.lt..01*tolf) then
        return
      endif
      sum=0.
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=stpmx*dmax1(dsqrt(sum),dble(n))
c *********** stpmax definition *************
      stpmax=x(1)*(1.-tolstp)
      xless=x(1)
      do istp=2,n
        if(x(istp).lt.xless) then
          xless=x(istp)
          stpmax=x(istp)*(1.-tolstp)
        endif
      enddo
c *********** stpmax definition *************
      do 21 its=1,maxits
        call fdjac(n,x,fvec,np,fjac)
        do 14 i=1,n
          sum=0.
          do 13 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
13        continue
          g(i)=sum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue
        call ludcmp(fjac,n,np,indx,d)
        call lubksb(fjac,n,np,indx,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.
        do 17 i=1,n
          if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
17      continue
        if(test.lt.tolf)then
          check=.false.
          return
        endif
        if(check)then
          test=0.
          den=dmax1(f,.5d0*n)
          do 18 i=1,n
            temp=dabs(g(i))*dmax1(dabs(x(i)),1.d0)/den
            if(temp.gt.test)test=temp
18        continue
          if(test.lt.tolmin)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.
        do 19 i=1,n
          temp=(dabs(x(i)-xold(i)))/dmax1(dabs(x(i)),1.d0)
          if(temp.gt.test)test=temp
19      continue
        if(test.lt.tolx) then
c          write(*,*) 'return from 4'
          return
        endif
21    continue
      write(*,*) 'ERROR in newt - maxits exceeded in newt'
      return
      end
      subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      implicit real*8 (a-h,o-z)
      parameter (neq=2,alf=1.d-6,tolx=1.d-8)
      common /step/ tolst
      logical check
      external func
CU    USES func
      dimension g(neq),p(neq),x(neq),xold(neq)
      check=.false.
      sum=0.
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=dsqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=tolx/test
      alam=1.
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+alf*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              tmplam=(-b+dsqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=dmax1(tmplam,.1*alam)
      goto 1
      end
      real*8 function fmin(x)
      implicit real*8 (a-h,o-z)
      parameter (np=40)
      common /newtv/ fvec(np),n
CU    USES funcv
      dimension x(*)
      call funcv(n,x,fvec)
      sum=0.
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5*sum
      return
      end
c
      subroutine fdjac(n,x,fvec,np,df)
      implicit real*8 (a-h,o-z)
      parameter (neq=2,nmax=40,eps=1.d-8)
CU    USES funcv
      dimension df(np,np),fvec(neq),x(neq),f(nmax)
      do 12 j=1,n
        temp=x(j)
        h=eps*dabs(temp)
        if(h.eq.0.)h=eps
c ORIGINAL
        x(j)=temp+h
c ORIGINAL
c VolatileCalc
c        x(j)=temp+eps
c VolatileCalc
        h=x(j)-temp
        temp2=x(j)
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
c ORIGINAL
          df(i,j)=(f(i)-fvec(i))/h
c ORIGINAL
c VolatileCalc
c          df(i,j)=(f(i)-fvec(i))/eps
c VolatileCalc
11      continue
12    continue
      return
      end
c
      subroutine lubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      parameter(neq=2)
      dimension a(np,np),b(neq),indx(neq)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end
c
      subroutine ludcmp(a,n,np,indx,d)
c
      implicit real*8 (a-h,o-z)
      parameter (neq=2,nmax=500,tiny=1.0d-20)
      common /temp/ its
      dimension a(np,np),vv(nmax),indx(neq)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.) then
           write(*,*) 'ERROR in newt - aamax.eq.0.'
           return
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=tiny
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      end
