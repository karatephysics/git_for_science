do a=1,int((int(DOSpara(1))-1)*1/2)
   do nky=a,int(3/4*(int(DOSpara(1))-1)-a/2)
         write(10,*) nky
      do nkx=nky,int(3/2*(int(DOSpara(1))-1)-a-nky)!.and.nkx<(DOSpara(1)-1)
