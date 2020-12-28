#ifndef PARAMETER_INCLUDED
#define PARAMETER_INCLUDED

#define MEM_THRESHOLD  (2*1024) /* 2048 MB : Algorithm switching threshold */
#define CHUNK                64 /* (multiple of sizeof(uint64_t)*8 for AVX-512). */
#define BLOCKS          (28*16) /* Only for GPU. */
#define THREADS         (64*16) /* Only for GPU. Must be 2^n */
#define MAX_HOSTNAME_LENGTH 256
#endif
