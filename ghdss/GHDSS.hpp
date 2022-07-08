#include <math.h>

#define ARY_NUM_MIC 12
#define ARY_NUM_SRC 8
#define ARY_NUM_FFT 256
//#define THRESHOLD (1000)
#define ha_iUpdate  (0)     // SS_METHOD
#define ha_fWmyu    (0)     //

typedef struct _HA_Complex_
{
    float re;
    float im;
} HA_Complex;
