#ifndef PARAMETER_INCLUDED
#define PARAMETER_INCLUDED

#define MEM_THRESHOLD  (2*1024) /* 2048 MB : Algorithm switching threshold */
#ifdef __FUJITSU
#define CPU_CHUNK            64
#define ALIGN_VALUE          64
#else
#define CPU_CHUNK            24
#define ALIGN_VALUE          32
#endif
#define GPU_CHUNK            64
#define BLOCKS          (28*16) /* Only for GPU. */
#define THREADS         (64*16) /* Only for GPU. Must be 2^n */
#define MAX_HOSTNAME_LENGTH 256
#define GEN_GRAPH_ITERS      10
#endif
