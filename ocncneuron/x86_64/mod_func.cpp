#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _tmgexp2syn_reg(void);
extern void _vecevent_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"tmgexp2syn.mod\"");
    fprintf(stderr, " \"vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _tmgexp2syn_reg();
  _vecevent_reg();
}

#if defined(__cplusplus)
}
#endif
