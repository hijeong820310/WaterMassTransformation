begin
ocn_fillvalue = 1e+30

DATA_diri= "/home/hijeong/SouthernOcean/IBS_data/WMT/densf/PD1000db/"
DATA_diro= "/home/hijeong/SouthernOcean/IBS_data/WMT/TimeVari/PD1000db/"
DATA_bins= "/home/hijeong/SouthernOcean/IBS_data/Transport/"
DATA_ht  = "/home/hijeong/SouthernOcean/IBS_data/"
 f_rho   = addfile(DATA_bins+"SIGMA1_bins.nc","r")
 Densbin = f_rho->sigma1_bins
 nbins    = dimsizes(Densbin)

do YR=1,70 
   YYYY = sprinti("%0.4i",YR)

  DATA_diroy = DATA_diro + YYYY
  system("mkdir -p "+DATA_diroy)

do MO=1,12 
   MM =sprinti("%2.2i",MO)

  fin   = addfile(DATA_diri+"PD.densf."+YYYY+"-"+MM+".nc","r")
  fht   = addfile(DATA_ht+"HT.90S-50S.nc","r")
  HT    = fht->HT(0:428,:) ; ~60S

  dens  = fin->rho(0:428,:)
  dens  = dens - 1000.

  io    = fin->io(0:428,:)
  brine = fin->brine(0:428,:)
  melt  = fin->melt(0:428,:)
  ao    = fin->ao(0:428,:)
  ioao  = fin->ioao(0:428,:)
  hap   = fin->hap(0:428,:)

  dens  = where(HT.ge.100000, 0, dens)
  io    = where(HT.ge.100000, 0, io)
  brine = where(HT.ge.100000, 0, brine)
  melt  = where(HT.ge.100000, 0, melt)
  ao    = where(HT.ge.100000, 0, ao)
  ioao  = where(HT.ge.100000, 0, ioao)
  hap   = where(HT.ge.100000, 0, hap) 

  dens@_FillValue = ocn_fillvalue

  bio    = io
  bbrine = io
  bmelt  = io
  bao    = io
  bioao  = io
  bhap   = io

  do di = 0, nbins-2
  bio     = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),io,0)
  bbrine  = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),brine,0)
  bmelt   = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),melt,0)
  bao     = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),ao,0)
  bioao   = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),ioao,0)
  bhap    = where(dens.ge.Densbin(di).and.dens.lt.Densbin(di+1),hap,0) 

  bio!0   = "lat"
  bio!1   = "lon"
  copy_VarCoords(bio,bbrine)
  copy_VarCoords(bio,bmelt) 
  copy_VarCoords(bio,bao)
  copy_VarCoords(bio,bioao) 
  copy_VarCoords(bio,bhap)

  system("rm -rf "+DATA_diroy+"/PD1000db.WMT."+YYYY+"-"+MM+".bin."+di+".nc")
  fout=addfile(DATA_diroy+"/PD1000db.WMT."+YYYY+"-"+MM+".bin."+di+".nc","c")
  fout->io=bio
  fout->brine=bbrine
  fout->melt=bmelt
  fout->ao=bao  
  fout->ioao=bioao 
  fout->hap=bhap

  delete(bio)
  delete(bbrine)
  delete(bmelt)
  delete(bao)
  delete(bioao)
  delete(bhap) 
  end do

   delete(fin)
end do
end do
end
