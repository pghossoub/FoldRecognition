#ifdef enablemultithread
#define TLS __thread
#else
#define TLS 
#endif

extern TLS int commonAlloc1;
extern TLS int commonAlloc2;
extern TLS int **commonIP;
extern TLS int **commonJP;
