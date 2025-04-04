#ifndef _STREAM_SRV_H
#define _STREAM_SRV_H

#ifndef _SERVICE_H
#include "service.h"
#endif

int srvInqContents(stream_t *streamptr);
int srvInqTimestep(stream_t *streamptr, int tsID);

int srvInqRecord(stream_t *streamptr, int *varID, int *levelID);
void srvDefRecord(stream_t *streamptr);
void srvCopyRecord(stream_t *streamptr2, stream_t *streamptr1);
void srv_read_record(stream_t *streamptr, int memtype, void *data, size_t *numMissVals);
void srv_write_record(stream_t *streamptr, int memtype, const void *data);

void srvReadVarDP(stream_t *streamptr, int varID, double *data, size_t *numMissVals);
void srvWriteVarDP(stream_t *streamptr, int varID, const double *data);

void srvReadVarSliceDP(stream_t *streamptr, int varID, int levelID, double *data, size_t *numMissVals);
void srvWriteVarSliceDP(stream_t *streamptr, int varID, int levelID, const double *data);

#endif /* _STREAM_SRV_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
