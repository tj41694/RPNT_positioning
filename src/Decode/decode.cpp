#include "Decode/decode.h"

/* decode data for kinds of formats ------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------- */
decode_data_c::decode_data_c() {
    format=-1;
    ephsat=0; dgps=NULL; Svr=NULL;
    /* class */
    obs.n=0; obs.data.assign(MAXOBS,gnss_obsd_c());
    nav.n=MAXSAT; nav.ng=MAXPRNGLO;
    nav.eph.assign(MAXSAT,gnss_eph_c()); nav.geph.assign(MAXPRNGLO,gnss_geph_c());

    for (int i=0; i<MAXSAT; i++) ssr[i]=gnss_ssr_c();
}
decode_data_c::~decode_data_c() {
    dgps=NULL; Svr=NULL;
}
int decode_data_c::decode(unsigned char data) {
    return 0;
}