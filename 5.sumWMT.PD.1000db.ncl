begin

   ocn_fillvalue = 1e+30

DATA_diri= "/home/hijeong/SouthernOcean/IBS_data/WMT/TimeVari/PD1000db/"
DATA_diro= "/home/hijeong/SouthernOcean/IBS_data/WMT/WMT/PD1000db/"
DATA_dirr= "/home/hijeong/SouthernOcean/IBS_data/"
DATA_bins= "/home/hijeong/SouthernOcean/IBS_data/Transport/"

 DATA_rr = addfile(DATA_dirr+"TAREA.90S-50S.nc","r")
 areaCell = DATA_rr->TAREA
 areaCell = areaCell/1e+4 ; cm2 to m2
 RAC      = areaCell
 RAC!0    = "lat"
 RAC!1    = "lon"
 
 f_rho   = addfile(DATA_bins+"SIGMA1_bins.nc","r")
 Denbins = f_rho->sigma1_bins

 drho    = new((/227/),"float")
 do di=0,226
 drho(di)=Denbins(di+1) - Denbins(di)
 end do

;======================
; Input Files
;======================

  do YR = 1,70 
  YYYY  = sprinti("%0.4i",YR)
  DATA_diriy = DATA_diri + YYYY  + "/"

  sumWTR   = new((/12,227/),"double",ocn_fillvalue)
  sumWTR!0 = "time"
  sumWTR!1 = "lev"
  sumio    = sumWTR
  sumbrine = sumWTR
  summelt  = sumWTR
  sumao    = sumWTR
  sumioao  = sumWTR
  sumhap   = sumWTR


  do MO=1,12
  MM =sprinti("%2.2i",MO)
  do di=0,226 

  fin    = addfile(DATA_diriy+"PD1000db.WMT."+YYYY+"-"+MM+".bin."+di+".nc","r")
  io     = fin->io
  brine  = fin->brine
  melt   = fin->melt
  ao     = fin->ao  
  ioao   = fin->ioao 
  hap    = fin->hap 

  io     = io*RAC/1e+6
  brine  = brine*RAC/1e+6
  melt   = melt*RAC/1e+6
  ao     = ao*RAC/1e+6
  ioao   = ioao*RAC/1e+6
  hap    = hap*RAC/1e+6

  iosum    =dim_sum(dim_sum(io))/drho(di)
  brinesum =dim_sum(dim_sum(brine))/drho(di)
  meltsum  =dim_sum(dim_sum(melt))/drho(di)
  aosum    =dim_sum(dim_sum(ao))/drho(di)
  ioaosum  =dim_sum(dim_sum(ioao))/drho(di)
  hapsum   =dim_sum(dim_sum(hap))/drho(di)

  sumio(MO-1,di)    = iosum
  sumbrine(MO-1,di) = brinesum
  summelt(MO-1,di)  = meltsum
  sumao(MO-1,di)    = aosum
  sumioao(MO-1,di)  = ioaosum
  sumhap(MO-1,di)   = hapsum

  delete(io)
  delete(brine)
  delete(melt)
  delete(ao)  
  delete(ioao)
  delete(hap) 
  end do
  end do

  system("rm -rf "+DATA_diro+"/PD1000db.sumWTR."+YYYY+".nc")
  fout=addfile(DATA_diro+"/PD1000db.sumWTR."+YYYY+".nc","c")
  fout->io=sumio
  fout->brine=sumbrine
  fout->melt=summelt
  fout->ao=sumao  
  fout->ioao=sumioao
  fout->hap=sumhap 

  end do
end 
