#include <cstdlib>

/* Using this file, you can call in FORTRAN anywhere

   CALL EXAHYPE_ABORT

in order to quit the program. Typically, you would like to
just CALL ABORT but some Fortran versions don't have this
function available. Calling ABORT is better than just EXIT
because it triggers a coredump which allows post mortem
error analysis (i.e. it gives you a stack trace with all
the local variables). */

extern "C" {
void exahype_abort_() {
	std::abort();
}
}
