#ifndef PARAMETER_INCLUDED
#define PARAMETER_INCLUDED

#define MEM_THRESHOLD  (2*1024) /* 2048 MB : Algorithm switching threshold */
#define CPU_CHUNK            24 /* (multiple of sizeof(uint64_t)*8 for AVX-512). */
#define ALIGN_VALUE          32 /* for posix_memalign() */
#define GPU_CHUNK            64
#define BLOCKS          (28*16) /* Only for GPU. */
#define THREADS         (64*16) /* Only for GPU. Must be 2^n */
#define MAX_HOSTNAME_LENGTH 256
#define GEN_GRAPH_ITERS      10
#endif
