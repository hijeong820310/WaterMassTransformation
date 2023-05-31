load "./Sub.EOS.drhodt.1000db.ncl"
load "./Sub.EOS.drhods.1000db.ncl"
load "./Sub.EOS.rho_1000db.ncl"

begin

   ocn_fillvalue = 1e+30

DATA_dir_fw = "/home/hijeong/SouthernOcean/IBS_data/DATA/PD/" 
DATA_dir_ts = "/home/hijeong/SouthernOcean/IBS_data/DATA/PD/"
DATA_diro   = "./densf/PD1000db/"

  DATA_t    = addfile(DATA_dir_ts+"PD.SST.year.month.nc","r")
  DATA_s    = addfile(DATA_dir_ts+"PD.SSS.year.month.nc","r")

  DATA_fw   = addfile(DATA_dir_fw+"PD.fw.year.month.nc","r")
  DATA_ao   = addfile(DATA_dir_fw+"PD.ao.year.month.nc","r")
  DATA_so   = addfile(DATA_dir_fw+"PD.so.year.month.nc","r")
  DATA_br   = addfile(DATA_dir_fw+"PD.siform.year.month.nc","r")
  DATA_mt   = addfile(DATA_dir_fw+"PD.simelt.year.month.nc","r")
  DATA_hf   = addfile(DATA_dir_fw+"PD.shf.year.month.nc","r")
  DATA_aice = addfile(DATA_dir_fw+"PD.aice.year.month.nc","r")

  do yy = 0,69
  do mm = 0, 11

   brine = DATA_br->siform(yy,mm,:,:)
   melt  = DATA_mt->simelt(yy,mm,:,:)
   IOAO  = DATA_fw->fw(yy,mm,:,:)
   AO    = DATA_ao->ao(yy,mm,:,:)
   IO    = DATA_so->so(yy,mm,:,:)
   shf   = DATA_hf->shf(yy,mm,:,:)
   temp  = DATA_t->TEMP(yy,mm,:,:)
   salt  = DATA_s->SALT(yy,mm,:,:)
   aice  = DATA_aice->aice(yy,mm,:,:)

   brine = brine*(-1)
   melt  = melt*(-1)
   AO    = AO*(-1)
   IO    = IO*(-1)
   IOAO  = IOAO*(-1)

  ;================================================
  ; Variable Setting
  ;================================================

   rUnit2mass      = 1.035e+03
   HeatCapacity_cp = 3.994e+03

   fluxfac_T  = 1.0/rUnit2mass/HeatCapacity_cp
   fluxfac_S  = 1.0/rUnit2mass

   Dens    = eosstat_1000_rho(temp,salt)
   drhodt  = eosstat_1000_drhodt(temp,salt)
   drhods  = eosstat_1000_drhods(temp,salt)

   Dens   = where(Dens.lt.1000.5,ocn_fillvalue,Dens)
   salt1  = where(Dens.lt.1000.5,ocn_fillvalue,salt)
   brine1 = where(Dens.lt.1000.5,ocn_fillvalue,brine)
   melt1  = where(Dens.lt.1000.5,ocn_fillvalue,melt)
   ao1    = where(Dens.lt.1000.5,ocn_fillvalue,AO)
   io1    = where(Dens.lt.1000.5,ocn_fillvalue,IO)
   ioao1  = where(Dens.lt.1000.5,ocn_fillvalue,IOAO)
   hap1   = where(Dens.lt.1000.5,ocn_fillvalue,shf)

   dens_brine= new((/603,3600/),"double",ocn_fillvalue)
   dens_melt = dens_brine
   dens_ioao = dens_brine
   dens_ao   = dens_brine
   dens_io   = dens_brine
   dens_hap  = dens_brine

   dens_io   = drhods * fluxfac_S * salt1 * io1
   dens_brine= drhods * fluxfac_S * salt1 * brine1
   dens_melt = drhods * fluxfac_S * salt1 * melt1
   dens_ao   = drhods * fluxfac_S * salt1 * ao1  
   dens_ioao = drhods * fluxfac_S * salt1 * ioao1
   dens_hap  = drhodt * fluxfac_T * hap1

   dens_io!0    = "lat"
   dens_io!1    = "lon"
   copy_VarCoords(dens_io,dens_brine)
   copy_VarCoords(dens_io,dens_melt) 
   copy_VarCoords(dens_io,dens_ao) 
   copy_VarCoords(dens_io,dens_ioao) 
   copy_VarCoords(dens_io,dens_hap)  
   copy_VarCoords(dens_io,Dens)  

   YYYY = sprinti("%4.4i",yy + 1)
   MM   = sprinti("%2.2i",mm + 1)

   system("rm -rf "+DATA_diro+"PD.densf."+YYYY+"-"+MM+".nc")
   fout=addfile(DATA_diro+"PD.densf."+YYYY+"-"+MM+".nc","c")
   fout->rho=Dens
   fout->io=dens_io
   fout->brine=dens_brine
   fout->melt=dens_melt
   fout->ao=dens_ao
   fout->ioao=dens_ioao
   fout->hap=dens_hap

   delete(Dens)
   delete(salt1)
   delete(io1)
   delete(brine1)
   delete(melt1)
   delete(ao1)
   delete(ioao1)
   delete(hap1)
   delete(drhodt)
   delete(drhods)
   delete(brine)
   delete(melt)
   delete(IOAO)
   delete(AO)
   delete(IO)
   delete(shf)
   delete(temp)
   delete(salt)
   delete(aice)
end do
end do
end 
