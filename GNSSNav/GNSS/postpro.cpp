#include "GNSS/gnss_pro.h"
#include "BaseFunction/basefunction.h"
#include "GNSS/PosModel/position.h"
/* post-processing server class ------------------------------------------------------------------- */
/* Constructor -------------------------------------------------------------------- */
gnss_postsvr_c::gnss_postsvr_c() {
    nav=new gnss_nav_c;
    obsBuffIndex[0]=obsBuffIndex[1]=0;
    obs[0]=NULL; obs[1]=NULL;
    epochObs[0]=epochObs[1]=gnss_obs_c();
    rbStaInf[0]=gnss_sta_c(); rbStaInf[0]=gnss_sta_c();
    solopt[0]=solopt[1]=gnss_solopt_c();
    solstream[0]=solstream[1]=file_c();
}
gnss_postsvr_c::~gnss_postsvr_c() {
    if (gnss_pro) { delete gnss_pro; gnss_pro=NULL; }
    if (nav) { delete nav; nav=NULL; }
    prcopt=NULL;
}
/* Implementation functions ------------------------------------------------------- */
/* set sat\rcv antenna information ------------------------------------------------ */
void gnss_postsvr_c::setpcv(in_gnss_atx_c *atx,int satrcv) {
    /* update satellite antenna to gnss_pro->nav */
    if (satrcv==0) {
        for (int i=0; i<MAXSAT; i++) {
            if (!(satsys(i+1,NULL)&gnss_pro->opt->navsys)) continue;
            for (int np=0; np<atx->satpcv.size(); np++) {
                if (atx->satpcv[np].sat!=i+1) continue; //test prn
                if (atx->satpcv[np].te.time>0&&atx->satpcv[np].te.timediff(pstopt->time_start)<0.0)
                    continue;
                nav->pcvs[i]=atx->satpcv[np];
                break;
            }
        }
    }
    /* update receiver antenna to gnss_pro->opt */
    else {
        gnss_pro->opt->pcvr[0]=atx->recpcv[0];
        gnss_pro->opt->pcvr[1]=atx->recpcv[1];
    }
}
/* new gnss_pro according to opt -------------------------------------------------- */
void gnss_postsvr_c::ini_gnss_pro(gnss_prcopt_c *Prcopt,gnss_filopt_t *Filopt) {
    if (gnss_pro) { delete gnss_pro; gnss_pro=NULL; }
    if (Prcopt->mode==PMODE_SINGLE) gnss_pro=new gnss_single_c;
    else if (Prcopt->mode>=PMODE_DGPS) gnss_pro=new gnss_relative_c;
    else  gnss_pro=new gnss_ppp_c;
    gnss_pro->opt=Prcopt;
    /* set base station position */
    for (int i=0; i<6; i++) {
        gnss_pro->rb[i]=i<3 ? Prcopt->rb[i] : 0.0;
    }

    /* set frequency priority according to configuration */
    nav->freq.change_freq_priority(Prcopt);

    /* initialize navigation and observation pointer */
    //gnss_pro->nav=nav; gnss_pro->obsr=obs[0]; gnss_pro->obsb=obs[1];
    gnss_pro->nav=nav; gnss_pro->obsr=&epochObs[0]; gnss_pro->obsb=&epochObs[1];
    /* read erp file */
    if (Filopt->erp.length()>10) {
        in_gnss_erp_c erp(Filopt->erp);
        erp.readerp(&nav->erp);
        gnss_pro->tidefunc.init_erp(nav);
        gnss_pro->satantfunc.init_erp(nav);
    }
    /* satellite antenna information file*/
    sat_ant=Filopt->satantp;
    /* receiver antenna information file */
    rec_ant=Filopt->rcvantp;
    /* station ocean-loading tide file */
    blq_sta=Filopt->blq;
    /* station SNX solution file path */
    snx_sta=Filopt->snx;

    /* read DCB file */
    if (Filopt->dcb.length()>10) {
        in_gnss_bias_c dcb(Filopt->dcb);
        dcb.readBias(nav);
    }

    /* read ionex tec grid file */
    if (Filopt->iono.length()>10) {
        in_gnss_ionex_c ion(Filopt->iono);
        ion.readIonex(nav);
    }

    /* read ztd file */
    if (Filopt->ztd.length()>10) {
        int iztd=-1;
        str2int(Filopt->ztd.substr(Filopt->ztd.size()-1),iztd);
        Filopt->ztd.pop_back();
        in_gnss_ztd_c ztd(iztd,Filopt->ztd);
        ztd.readZTD(nav);
    }
    
    gnss_pro->log_path=Filopt->test;
}
/* update station information to prcopt ------------------------------------------- */
void gnss_postsvr_c::updatesta() {
    /* rover station information */
    prcopt->name[0]=rbStaInf[0].name.substr(0,4);
    prcopt->anttype[0]=rbStaInf[0].antdes;
    for (int i=0; i<3; i++) prcopt->antdel[0][i]=rbStaInf[0].del[i];
    /* station position */
    if ( prcopt->rovpos==POSOPT_RINEX && norm(rbStaInf[0].pos,3)>1E5 ) {
        prcopt->ru[0]=rbStaInf[0].pos[0];
        prcopt->ru[1]=rbStaInf[0].pos[1];
        prcopt->ru[2]=rbStaInf[0].pos[2];
    }
    /* base station information */
    if (prcopt->mode>=PMODE_DGPS) {
        prcopt->name[1]=rbStaInf[1].name.substr(0,4);
        prcopt->anttype[1]=rbStaInf[1].antdes;
        for (int i=0; i<3; i++) prcopt->antdel[1][i]=rbStaInf[1].del[i];
        /* station position */
        if ( prcopt->baspos==POSOPT_RINEX && norm(rbStaInf[1].pos,3)>1E5 ) {
            prcopt->rb[0]=rbStaInf[1].pos[0];
            prcopt->rb[1]=rbStaInf[1].pos[1];
            prcopt->rb[2]=rbStaInf[1].pos[2];
        }
    }
    epochObs[0].sta=rbStaInf[0]; epochObs[1].sta=rbStaInf[1];

    /* receiver antenna information and ocean-loading correction */
    if (rec_ant.length()>10) {
        in_gnss_atx_c rcv(rec_ant);
        rcv.read_recatx(prcopt,nav);
        setpcv(&rcv,1);
    }

    /* read blq file */
    if (blq_sta.length()>10&&prcopt->tidecorr&2) {
        in_gnss_blq_c blq(blq_sta);
        blq.readblq(prcopt->name[0],nav->ocean_par[0]);
        if (prcopt->mode>=PMODE_DGPS)
            blq.readblq(prcopt->name[1],nav->ocean_par[1]);
    }

    /* read snx file */
    if ( snx_sta.length()>10&& ( prcopt->rovpos==POSOPT_SNX || prcopt->baspos==POSOPT_SNX ) ) {
        in_gnss_snx_c snx(snx_sta);
        snx.readSNX(prcopt,rbStaInf);
    }

    /* tidal displacement functions */
    gnss_pro->tidefunc.init_otl(gnss_pro->opt,nav);
}
/* initialize read stream --------------------------------------------------------- */
int gnss_postsvr_c::ini_Read_Stream() {
    /* initialize rover_obs and navigation read stream */
    rover.ini_ReadO(pstopt->rover_obs,"",pstopt->time_start,pstopt->time_end,pstopt->time_inter,0);
    readN.ini_ReadN(pstopt->nav,"");
    /* read head of rover observation file */
    rover.Head(nav,&rbStaInf[0]);
    /* update post option start or end time if == 0 */
    if (pstopt->time_start.time<=0) {
        pstopt->time_start=rover.epochFirst;
        pstopt->time_start.ep[3]=pstopt->time_start.ep[4]=pstopt->time_start.ep[5]=0;
        pstopt->time_start.epoch2time(pstopt->time_start.ep);
        pstopt->time_start.time2str(3);
    }
    if (pstopt->time_end.time<=0) {
        pstopt->time_end=rover.epochFirst;
        pstopt->time_end.ep[3]=pstopt->time_end.ep[4]=pstopt->time_end.ep[5]=0;
        pstopt->time_end.epoch2time(pstopt->time_end.ep);
        pstopt->time_end.timeadd(86400.0);
        pstopt->time_end.time2str(3);
    }
    /* read navigation file */
    readN.Head(nav); readN.Body(nav,prcopt); readN.closeF();

    /* open stat */
    int stat=rover.test_open();

    /* initailize base_obs read stream */
    if ((prcopt->mode>=PMODE_DGPS)&&pstopt->base_obs.length()>=12) {
        base.ini_ReadO(pstopt->base_obs,"",pstopt->time_start,pstopt->time_end,pstopt->time_inter,1);
        /* read head of base observation file */
        base.Head(nav,&rbStaInf[1]);
        stat=stat&&base.test_open();
        if (stat) {
            if (prcopt->freqPriority) {
                nav->freq.change_gnss_freq(SYS_GPS,&rover.index[iGPS],&base.index[iGPS]);
                nav->freq.change_gnss_freq(SYS_GLO,&rover.index[iGLO],&base.index[iGLO]);
                nav->freq.change_gnss_freq(SYS_GAL,&rover.index[iGAL],&base.index[iGAL]);
                nav->freq.change_gnss_freq(SYS_BDS,&rover.index[iBDS],&base.index[iBDS]);
            }
            rover.set_freq_index(nav);
            base.set_freq_index(nav);
        }
    }
    else if (stat) {
        if (prcopt->freqPriority) {
            nav->freq.change_gnss_freq(SYS_GPS,&rover.index[iGPS],NULL);
            nav->freq.change_gnss_freq(SYS_GLO,&rover.index[iGLO],NULL);
            nav->freq.change_gnss_freq(SYS_GAL,&rover.index[iGAL],NULL);
            nav->freq.change_gnss_freq(SYS_BDS,&rover.index[iBDS],NULL);
        }
        rover.set_freq_index(nav);
    }

    /* update station information */
    updatesta();
    /* update satellite carrier wave lengths */
    nav->update_sat_lambda(gnss_pro->ssat);

    /* read sat antenna information file */
    if (sat_ant.length()>10) {
        in_gnss_atx_c sat(sat_ant);
        sat.read_satatx(nav);
        setpcv(&sat,0);
    }

    /* read precise ephemeris */
    if (prcopt->sateph==EPHOPT_PREC&&pstopt->prseph.length()>=12) {
        readEph.ini_readEph(pstopt->prseph,pstopt->predict);
        readEph.readsp3(nav,prcopt);
    }

    /* read precise clocks */
    if (prcopt->sateph==EPHOPT_PREC&&pstopt->prsclk.length()>=12) {
        readC.ini_ReadC(pstopt->prsclk,"");
        readC.readclk(nav,prcopt);
    }

    /* solution types */
    int year,doy;
    string yer_str,doy_str,suf=satsys_flag(prcopt->navsys);
    size_t fileDot;
    year=int(pstopt->time_start.ep[0]); doy=int(pstopt->time_start.time2doy());
    int2str(4,"0",year,yer_str); int2str(3,"0",doy,doy_str);
    for (int i=0; i<2; i++) {
        if ((fileDot=pstopt->output[i].find_last_of("."))!=string::npos) {
            pstopt->output[i].insert(fileDot,"_"+suf+yer_str+doy_str);
        }
        if (solopt[i].outopt&&!solstream[i].StreamOpen(pstopt->output[i].c_str(),STR_FILE,STR_MODE_W)) {
            for (i--; i>=0; i--) solstream[i].StreamClose();
            return 0;
        }
    }
    /* open test file */
    if (gnss_pro->log_path.length()>5) {
        if ((fileDot=gnss_pro->log_path.find_last_of("."))!=string::npos) {
            gnss_pro->log_path.insert(fileDot,"_"+suf+yer_str+doy_str);
        }
        gnss_pro->log_stream.open(gnss_pro->log_path,ios::out);
    }

    /* wirte solution file header */
    writesolhead();

    return stat;
}
/* rover/base observation synchronization ----------------------------------------- */
int gnss_postsvr_c::obs_synchron() {
    if (!base.test_open()) return 0;

    /* test this base observation time */
    if (epochObs[1].n>0
        &&fabs(epochObs[0].data[0].time.timediff(epochObs[1].data[0].time))<DTTOL) 
        return 1;
    else if (epochObs[1].n>0
        &&epochObs[0].data[0].time.timediff(epochObs[1].data[0].time)<-DTTOL) 
        return 0;
    /* read new base observation */
    //while (Read_One_Epoch(1)) {
    while(base.Read_One_Epoch(epochObs+1,prcopt,pstopt)){
        // test time sychronization of rover and base
        if ( epochObs[1].n>0
            &&fabs(epochObs[0].data[0].time.timediff(epochObs[1].data[0].time))<DTTOL ) 
            return 1;
        else if ( epochObs[1].n>0
            &&epochObs[0].data[0].time.timediff(epochObs[1].data[0].time)<-DTTOL ) 
            break;
    }

    return 0;
}
/* write solution header to output stream ----------------------------------------- */
void gnss_postsvr_c::writesolhead() {
    unsigned char buff1[8192]={ 0 };
    unsigned char buff2[8192]={ 0 };
    int n;

    if (solopt[0].outopt) {
        n=solopt[0].outsolheads(prcopt,buff1);
        solstream[0].StreamWrite(buff1,n);
    }
    if (solopt[1].outopt) {
        n=solopt[1].outsolheads(prcopt,buff2);
        solstream[1].StreamWrite(buff2,n);
    }
}
/* write solution to output stream ------------------------------------------------ */
void gnss_postsvr_c::writesol() {
    if (solopt[0].outopt) writesolstr(0);
    if (solopt[1].outopt) writesolstr(1);
}
/* write solution to each solstream ----------------------------------------------- */
void gnss_postsvr_c::writesolstr(int index) {
    unsigned char buff[MAXSOLMSG+1]={ 0 };
    char *p=(char *)buff;

    /* write solution to buff */
    /* [1] write solution time */
    string soltime=gnss_pro->sol.back().fortime(solopt+index);
    p+=sprintf(p,"%s%s",soltime.c_str(),solopt[index].sep.c_str());

    /* [2] position solution */
    string strpv;
    if ( gnss_pro->sol.back().stat==SOLQ_NONE ) return;
    else if (solopt[index].posf==SOLF_ENU) {
        if ( solopt[index].refsta==REFSAT_BASE || prcopt->mode==PMODE_MOVEB ) {
            for (int i=0; i<3; i++) {
                gnss_pro->sol.back().refpos[i]=gnss_pro->rb[i];
            }
        }
        else {
            for (int i=0; i<3; i++) {
                gnss_pro->sol.back().refpos[i]=prcopt->ru[i];
            }
        }
        if ( norm(gnss_pro->sol.back().refpos,3)<1E5 ) return;
    }
    strpv=gnss_pro->sol.back().forposvel(solopt+index,gnss_pro);
    p+=sprintf(p,"%s\n",strpv.c_str());

    /* write solution to stream[index+3] */
    solstream[index].StreamWrite(buff,p-(char *)buff);
}
/* read observation file from opbsBuff -------------------------------------------- */
int gnss_postsvr_c::Read_One_Epoch(int rovbas) {
    /* write obsBuff[rovbas] from observation file */
    if (obsBuff->size()<=0||obsBuffIndex[rovbas]>=obsBuff[rovbas].size()) {
        if (rovbas==0) { 
            if (rover.Read2_Obs_Buff(obsBuff[0],prcopt,pstopt)<=0) return 0;
            obsBuffIndex[0]=0; 
        }
        else { 
            if (base.Read2_Obs_Buff(obsBuff[1],prcopt,pstopt)<=0) return 0; 
            obsBuffIndex[1]=0; 
        }
    }
    /* make obs[rovbas] point to the current epoch in obsBuff[rovbas] */
    obs[rovbas]=&obsBuff[rovbas][obsBuffIndex[rovbas]++];
    obs[rovbas]->sta=rbStaInf[rovbas];

    return obs[rovbas]->n; 
}
/* initialize postsvr ------------------------------------------------------------- */
int gnss_postsvr_c::postsvrini(all_option_c *option) {
    /* option pointor */
    prcopt=&option->prcopt;
    pstopt=&option->pstopt;

    /* solution types */
    for (int i=0; i<2; i++) {
        solopt[i]=option->solopt[i];
    }
    /* initialize gnss_pro */
    ini_gnss_pro(prcopt,&option->filopt);
    gnss_pro->gnss_pro_init();

    return 1;
}
/* reset postsvr ------------------------------------------------------------------ */
void gnss_postsvr_c::postreset() {
    /* close output files */
    for (int i=0; i<2; i++) {
        solstream[i].StreamClose();
    }
    if (gnss_pro) gnss_pro->log_stream.close();
    /* close input files */
    rover.closeF(); base.closeF();
    readN.closeF(); readEph.closeF(); readC.closeF();
}
/* read Navigation file ----------------------------------------------------------- */
int gnss_postsvr_c::readNav() {
    /* initialize navigation file stream */
    readN.ini_ReadN(pstopt->nav,"");
    /* read navigation file */
    readN.Head(nav); readN.Body(nav,prcopt); readN.closeF();
    return 1;
}
/* post-position epoch by epoch --------------------------------------------------- */
int gnss_postsvr_c::Post_Position_Epoch() {
    /* read observation and navigation file */
    if (!ini_Read_Stream()) { postreset(); return 0; }

    /* read observation file and post-position epoch by epoch */
    //while (Read_One_Epoch(0)) {
    while(rover.Read_One_Epoch(epochObs,prcopt,pstopt)){
        if ((prcopt->mode>=PMODE_DGPS)&&!obs_synchron()) continue;

        /* debug */
        if (!epochObs[0].data[0].time.sep.compare("2022/04/24 04:40:40.000"))
            gnss_pro->obsr=&epochObs[0];

        /* SPP for base station */
        if (gnss_pro->opt->mode>=PMODE_DGPS) {
            if (!gnss_pro->basepos()) continue;
            if (gnss_pro->opt->mode==PMODE_MOVEB) {
                for (int i=0; i<3; i++) {
                    gnss_pro->rb[i]=gnss_pro->b_sol.back().xpos[i];
                }
            }
        }

        /* position process */
        if (gnss_pro->gnss_pos()) writesol();
    }

    postreset();

    return 1;
}
