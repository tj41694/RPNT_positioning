#ifndef STREAM_H
#define STREAM_H

#include "BaseFunction/timesys.h"

/* macros ----------------------------------------------------------------------------------------- */

#ifdef WIN32
#define dev_t               HANDLE
#define socket_t            SOCKET
typedef int socklen_t;
#else
#define dev_t               int
#define socket_t            int
#define closesocket         close
#endif

/* proto types for static functions ------------------------------------------*/
/* TCP and NTRIPT date structure */
typedef struct {                        /* tcp control type */
    int state;                          /* state (0:close,1:wait,2:connect) */
    string saddr;                       /* address string */
    int port;                           /* port */
    struct sockaddr_in addr;            /* address resolved */
    socket_t sock;                      /* socket descriptor */
    int tcon;                           /* reconnect time (ms) (-1:never,0:now) */
    unsigned int tact;                  /* data active tick */
    unsigned int tdis;                  /* disconnect tick */
} tcp_t;
/* ntrip caster connection type */
typedef struct {
    int state;                          /* state (0:close,1:connect) */
    string mntpnt;                      /* mountpoint (256) */
    string str;                         /* mountpoint string for server (256) */
    int nb;                             /* request buffer size */
    unsigned char buff[32768];          /* request buffer = NTRIP_MAXRSP in stream.cpp */
} ntripc_con_t;

/* generate binary raw message functions ---------------------------------------------------------- */
int gen_ubx(const char *msg,unsigned char *buff);
int gen_stq(const char *msg,unsigned char *buff);
int gen_nvs(const char *msg,unsigned char *buff);
int gen_lexr(const char *msg,unsigned char *buff);
int gen_hex(const char *msg,unsigned char *buff);

/* stream class ----------------------------------------------------------------------------------- */
class stream_c {
/* Constructor */
public:
    stream_c();
    virtual ~stream_c();
/* Implementation functions */
protected:
    /* lock/unlock stream ----------------------------------------------- */
    void strlock();
    void strunlock();
    /* set bitrate ------------------------------------------------------------ */
    int set_brate(int brate);
    /* TCP functions ---------------------------------------------------------- */
    /* decode tcp/ntrip path (path=[user[:passwd]@]addr[:port][/mntpnt[:str]])- */
    void decodetcppath(const string path,string *Addr,string *Port,string *User,
        string *Passwd,string *Mntpnt,string *Strp);
    /* set socket option ------------------------------------------------------ */
    int setsock(socket_t sock);
    /* non-block accept ------------------------------------------------------- */
    socket_t accept_nb(socket_t sock,struct sockaddr *Saddr,socklen_t *len);
    /* non-block connect ------------------------------------------------------ */
    int connect_nb(socket_t sock,struct sockaddr *Saddr,socklen_t len);
    /* non-block receive ------------------------------------------------------ */
    int recv_nb(socket_t sock,unsigned char *buff,int n);
    /* non-block send --------------------------------------------------------- */
    int send_nb(socket_t sock,unsigned char *buff,int n);
    /* generate tcp socket ---------------------------------------------- */
    int gentcp(tcp_t &tcp,int type);
    /* disconnect tcp --------------------------------------------------- */
    void discontcp(tcp_t &tcp,int tcon);


public:
    /* special functions ------------------------------------------------------ */
    /* get stream from file_c time (only for file_c) -------------------- */
    virtual gtime_c strgettime();
    /* sync file_c stream with another (only for file_c) ---------------- */
    virtual void strsync(void *str_file);

    /* common functions ------------------------------------------------- */
    /* send receiver command -------------------------------------------- */
    void SendCmd(const char *cmd);
    /* open stream from kinds of format --------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close stream from kinds of format -------------------------------- */
    virtual void StreamClose();
    /* read stream from kinds of format --------------------------------- */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write stream from kinds of format -------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);

/* Components */
public:
    int Stype;                          /* type (STR_???) */
    int Smode;                          /* mode (STR_MODE_?) */
    int Sstate;                         /* state (-1:error,0:close,1:open) */
    unsigned int inb,inr;               /* input bytes/rate */
    unsigned int outb,outr;             /* output bytes/rate */
    unsigned int tick_i;                /* input tick tick */
    unsigned int tick_o;                /* output tick */
    unsigned int tact;                  /* active tick */
    unsigned int inbt,outbt;            /* input/output bytes at tick */
    lock_t Slock;                       /* lock flag */
    string path;                        /* stream path */
    string msg;                         /* stream message */
};

/* type definition -------------------------------------------------------------------------------- */
/* stream :: file control class ------------------------------------------------------------------- */
class file_c : public stream_c {
/* Constructor */
public:
    file_c();
    ~file_c();
/* Implementation Functions */
protected:
    /* open file -------------------------------------------------------- */
    int openfile_(gtime_c ttime);
    /* open new swap file ----------------------------------------------- */
    void swapfile(gtime_c ttime);
    /* close old swap file ---------------------------------------------- */
    void swapclose();
public:
    /* special functions ------------------------------------------------------ */
    /* get stream from file_c time (only for file_c) -------------------------- */
    virtual gtime_c strgettime();
    /* sync file_c stream with another (only for file_c) ---------------------- */
    virtual void strsync(void *str_file);

    /* common functions ------------------------------------------------------- */
    /* open from file download ------------------------------------------------ */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from file download ----------------------------------------------- */
    virtual void StreamClose();
    /* read from file download ------------------------------------------------ */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from file download ----------------------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    FILE *fp;                           /* file pointer */
    FILE *fp_tag;                       /* file pointer of tag file */
    FILE *fp_tmp;                       /* temporary file pointer for swap */
    FILE *fp_tag_tmp;                   /* temporary file pointer of tag file for swap */
    string fpath;                       /* file path */
    string openpath;                    /* open file path */
    int mode;                           /* file mode */
    int timetag;                        /* time tag flag (0:off,1:on) */
    int repmode;                        /* replay mode (0:master,1:slave) */
    int offset;                         /* time offset (ms) for slave */
    gtime_c time;                       /* start time */
    gtime_c wtime;                      /* write time */
    unsigned int tick;                  /* start tick */
    unsigned int tick_f;                /* start tick in file */
    unsigned int fpos;                  /* current file position */
    double start;                       /* start offset (s) */
    double speed;                       /* replay speed (time factor) */
    double swapintv;                    /* swap interval (hr) (0: no swap) */
    lock_t lock;                        /* lock flag */
};
/* stream :: tcp server class --------------------------------------------------------------------- */
class tcpsvr_c : public virtual stream_c {
/* Constructor */
public:
    tcpsvr_c();
    virtual ~tcpsvr_c();
/* Implementation Functions */
protected:
    /* update tcp server ------------------------------------------------ */
    void updatetcpsvr();
    /* accept client connection ----------------------------------------- */
    int accsock();
    /* wait socket accept ----------------------------------------------- */
    int waittcpsvr();
public:
    /* open from tcp server download ------------------------------------------ */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from tcp server download ----------------------------------------- */
    virtual void StreamClose();
    /* read from tcp server download ------------------------------------------ */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from tcp server download ----------------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    tcp_t svr;                          /* tcp server control */
    tcp_t cli[32];                      /* tcp client controls  =MAXCLI in stream.cpp */
};
/* stream :: tcp cilent class --------------------------------------------------------------------- */
class tcpcli_c : public virtual stream_c {
/* Constructor */
public:
    tcpcli_c();
    virtual ~tcpcli_c();
/* Implementation Functions */
protected:
    /* connect server --------------------------------------------------- */
    int consock();
    /* wait socket connect ---------------------------------------------- */
    int waittcpcli();
public:
    /* open from tcp cilent download ------------------------------------------ */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from tcp cilent download ----------------------------------------- */
    virtual void StreamClose();
    /* read from tcp cilent download ------------------------------------------ */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from tcp cilent download ----------------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    tcp_t svr;                          /* tcp server control */
    int toinact;                        /* inactive timeout (ms) (0:no timeout) */
    int tirecon;                        /* reconnect interval (ms) (0:no reconnect) */
};
/* stream :: tcpsvr_c :: serial control class ----------------------------------------------------- */
class serial_c : public tcpsvr_c {
/* Constructor */
public:
    serial_c();
    ~serial_c();
/* Implementation Functions */
protected:
    /* write serial data ------------------------------------------------------ */
    int writeseribuff(unsigned char *bbuff,int n);
public:
    /* read serial data ------------------------------------------------------- */
    int readseribuff(unsigned char *rbuff,int nmax);
    /* open from serial download ---------------------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from serial download --------------------------------------------- */
    virtual void StreamClose();
    /* read from serial download ---------------------------------------------- */
    virtual int StreamRead(unsigned char *bbuff,int n);
    /* write from serial download --------------------------------------------- */
    virtual int StreamWrite(unsigned char *bbuff,int n);
/* Components */
public:
    dev_t dev;                          /* serial device */
    int error;                          /* error state */
    int state,wp,rp;                    /* state,write/read pointer */
    int buffsize;                       /* write buffer size (bytes) */
    lock_t lock;                        /* lock flag */
    unsigned char *buff;                /* write buffer */
#ifdef WIN32
    HANDLE thread;                      /* write thread */
#endif
};
/* stream :: tcpcli_c :: ntrip control class ------------------------------------------------------ */
class ntrip_c : public tcpcli_c {
/* Constructor */
public:
    ntrip_c();
    ~ntrip_c();
/* Implementation Functions */
protected:
    /* NTRIP functions -------------------------------------------------- */
    /* send ntrip server request ---------------------------------------- */
    int reqntrip_s();
    /* send ntrip client request ---------------------------------------- */
    int reqntrip_c();
    /* send ntrip server response --------------------------------------- */
    int rspntrip_s();
    /* send ntrip client response --------------------------------------- */
    int rspntrip_c();
    /* wait ntrip request/response -------------------------------------- */
    int waitntrip();
public:
    /* open from ntrip download ----------------------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from ntrip download ---------------------------------------------- */
    virtual void StreamClose();
    /* read from ntrip download ----------------------------------------------- */
    virtual int StreamRead(unsigned char *bbuff,int n);
    /* write from ntrip download ---------------------------------------------- */
    virtual int StreamWrite(unsigned char *bbuff,int n);
/* Components */
public:
    int state;                          /* state (0:close,1:wait,2:connect) */
    int type;                           /* type (0:server,1:client) */
    int nb;                             /* response buffer size */
    string url;                         /* url for proxy (256) */
    string mntpnt;                      /* mountpoint (256) */
    string user;                        /* user (256) */
    string passwd;                      /* password (256) */
    string str;                         /* mountpoint string for server (256) */
    unsigned char buff[32768];          /* response buffer = NTRIP_MAXRSP in stream.cpp */
};

/* stream :: tcpsvr_c :: ntrip caster control class ----------------------------------------------- */
class ntripc_c : public tcpsvr_c {
/* Constructor */
public:
    ntripc_c();
    ~ntripc_c();
/* Implementation Functions */
protected:
    /* disconnect ntrip-caster connection ------------------------------- */
    void discon_ntripc(int i);
    /* test mountpoint in source table ---------------------------------- */
    int test_mntpnt(const char *mmntpnt);
    /* send ntrip source table ------------------------------------------ */
    void send_srctbl(socket_t sock);
    /* test ntrip-caster client request --------------------------------- */
    void rsp_ntripc_c(int i);
    /* test ntrip-caster client request --------------------------------- */
    void rsp_ntripc_s(int i);
    /* handle ntrip-caster connect request ------------------------------ */
    void wait_ntripc();
public:
    /* open from ntrip caster download ---------------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from ntrip caster download --------------------------------------- */
    virtual void StreamClose();
    /* read from ntrip caster download ---------------------------------------- */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from ntrip caster download --------------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    int state;                          /* state (0:close,1:wait,2:connect) */
    int type;                           /* type (0:server,1:client) */
    string mntpnt;                      /* selected mountpoint (256) */
    string user;                        /* user (256) */
    string passwd;                      /* password (256) */
    char *srctbl;                       /* source table */
    lock_t lock_srctbl;                 /* lock flag for source table */
    ntripc_con_t con[32];               /* ntrip caster connections =MAXCLI in stream.cpp */
};
/* stream :: udp class ---------------------------------------------------------------------------- */
class udp_c : public stream_c {
/* Constructor */
public:
    udp_c();
    ~udp_c();
/* Implementation Functions */
protected:
    /* generate udp socket ---------------------------------------------------- */
    int genudp(int udptype,int strport,string Saddr);
public:
    /* open from udp download ------------------------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from udp download ------------------------------------------------ */
    virtual void StreamClose();
    /* read from udp download ------------------------------------------------- */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from udp download ------------------------------------------------ */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    int state;                          /* state (0:close,1:open) */
    int type;                           /* type (0:server,1:client) */
    int port;                           /* port */
    string saddr;                       /* address (server:filter,client:server) (256) */
    struct sockaddr_in addr;            /* address resolved */
    socket_t sock;                      /* socket descriptor */
};
/* stream :: ftp download control class ----------------------------------------------------------- */
class ftp_c : public stream_c {
/* Constructor */
public:
    ftp_c();
    ~ftp_c();
/* Implementation Functions */
protected:
    /* decode ftp path -------------------------------------------------------- */
    void decodeftppath();
public:
    /* open from ftp download ------------------------------------------------- */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from ftp download ------------------------------------------------ */
    virtual void StreamClose();
    /* read from ftp download ------------------------------------------------- */
    virtual int StreamRead(unsigned char *buff,int n);
/* Components */
public:
    int state;                          /* state (0:close,1:download,2:complete,3:error) */
    int proto;                          /* protocol (0:ftp,1:http) */
    int error;                          /* error code (0:no error,1-10:wget error,
                                         * 11:no temp dir,12:uncompact error) */
    string addr;                        /* download address (1024) */
    string file;                        /* download file path (1024) */
    string user;                        /* user for ftp (256) */
    string passwd;                      /* password for ftp (256) */
    string local;                       /* local file path (1024) */
    int topts[4];                       /* time options {poff,tint,toff,tretry} (s) */
    gtime_c tnext;                      /* next retry time (gpst) */
    thread_t thread;                    /* download thread */
};
/* stream :: memory buffer class ------------------------------------------------------------------ */
class membuf_c : public stream_c {
/* Constructor */
public:
    membuf_c();
    ~membuf_c();
/* Implementation Functions */
public:
    /* open from memory buffer ------------------------------------------------ */
    virtual int StreamOpen(const char *strpath,int Strtype,int Strmode);
    /* close from memory buffer ----------------------------------------------- */
    virtual void StreamClose();
    /* read from memory buffer ------------------------------------------------ */
    virtual int StreamRead(unsigned char *buff,int n);
    /* write from memory buffer ----------------------------------------------- */
    virtual int StreamWrite(unsigned char *buff,int n);
/* Components */
public:
    int state,wp,rp;                    /* state,write/read pointer */
    int bufsize;                        /* buffer size (bytes) */
    lock_t lock;                        /* lock flag */
    unsigned char *buf;                 /* write buffer */
};

///* stream converter class ------------------------------------------------------------------------- */
//class strconv_t{
///* Constructor */
//    public:
//        strconv_t();
//        ~strconv_t();
///* Components */
//    public:
//        int itype,otype;              /* input and output stream type */
//        int nmsg;                     /* number of output messages */
//        int msgs[32];                 /* output message types */
//        double tint[32];              /* output message intervals (s) */
//        unsigned int tick[32];        /* cycle tick of output message */
//        int ephsat[32];               /* satellites of output ephemeris */
//        int stasel;                   /* station info selection (0:remote,1:local) */
//        rtcm_c *rtcm;                 /* rtcm input data buffer */
//        raw_c *raw;                   /* raw  input data buffer */
//        rtcm_c *out;                  /* rtcm output data buffer */
//};
///* stream server class ---------------------------------------------------------------------------- */
//class strsvr_t{
///* Constructor */
//    public:
//        strsvr_t();
//        ~strsvr_t();
///* Components */
//    public:
//        int state;                    /* server state (0:stop,1:running) */
//        int cycle;                    /* server cycle (ms) */
//        int buffsize;                 /* input/monitor buffer size (bytes) */
//        int nmeacycle;                /* NMEA request cycle (ms) (0:no) */
//        int relayback;                /* relay back of output streams (0:no) */
//        int nstr;                     /* number of streams (1 input + (nstr-1) outputs */
//        int npb;                      /* data length in peek buffer (bytes) */
//        string cmds_periodic[16];     /* periodic commands */
//        double nmeapos[3];            /* NMEA request position (ecef) (m) */
//        unsigned char *buff;          /* input buffers */
//        unsigned char *pbuf;          /* peek buffer */
//        unsigned int tick;            /* start tick */
//        stream_c *stream[16];         /* input/output streams */
//        strconv_t *conv[16];          /* stream converter */
//        thread_t thread;              /* server thread */
//        lock_t lock;                  /* lock flag */
//};

#endif

