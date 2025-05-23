// Automatically generated by make_fint.c, don't edit!

// clang-format off

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef CDI_H_
#include "cdi.h"
#endif

#ifdef HAVE_CF_INTERFACE

#include <limits.h>
#include <assert.h>

#ifndef __CFORTRAN_LOADED
#  if defined __clang__
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wreserved-id-macro"
#  endif
#  include "cfortran.h"
#  if defined __clang__
#    pragma GCC diagnostic pop
#  endif
#endif
/* These functions are meant to be called from Fortran and don't
 * need an interface declaration in a C header. */
#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

#ifdef CDI_H_

static inline
int SizeType_c2f(SizeType value_SizeType)
{
  assert(value_SizeType < INT_MAX);
  return (int) value_SizeType;
}

#endif

/*  Start of fortran interface for the following routines (make_fint.c)  */


/*  Byte order  */


/*  Error identifier  */


/*  File types  */


/*  Compatibility defines for release 1.8.3 (obsolete defines)  */


/*  Protocols (in filename/URI)  */


/*  Compress types  */


/*  external data types  */


/*  Compatibility defines for release 1.8.3 (obsolete defines)  */


/*  internal data types  */


/*  Chunk types  */


/*  GRID types  */


/*  ZAXIS types  */


/*  SUBTYPE types  */


/*  Data structure defining a key-value search, possibly with multiple
   key-value pairs in combination.

   Currently, only multiple pairs combined by AND are supported.  */


/*  TIME types  */


/*  TSTEP types  */


/*  TAXIS types  */


/*  TUNIT types  */


/*  CALENDAR types  */


/*  number of unsigned char needed to store UUID  */


/*  Structs that are used to return data to the user  */


/*  Opaque types  */


/*  CDI control routines  */

FCALLSCSUB0 (cdiReset, CDIRESET, cdireset)
FCALLSCFUN1 (STRING, cdiStringError, CDISTRINGERROR, cdistringerror, INT)
FCALLSCSUB1 (cdiDebug, CDIDEBUG, cdidebug, INT)
FCALLSCFUN0 (STRING, cdiLibraryVersion, CDILIBRARYVERSION, cdilibraryversion)
FCALLSCSUB0 (cdiPrintVersion, CDIPRINTVERSION, cdiprintversion)
FCALLSCFUN1 (INT, cdiHaveFiletype, CDIHAVEFILETYPE, cdihavefiletype, INT)
FCALLSCSUB1 (cdiDefMissval, CDIDEFMISSVAL, cdidefmissval, DOUBLE)
FCALLSCFUN0 (DOUBLE, cdiInqMissval, CDIINQMISSVAL, cdiinqmissval)
FCALLSCSUB2 (cdiDefGlobal, CDIDEFGLOBAL, cdidefglobal, STRING, INT)
FCALLSCFUN0 (INT, namespaceNew, NAMESPACENEW, namespacenew)
FCALLSCSUB1 (namespaceSetActive, NAMESPACESETACTIVE, namespacesetactive, INT)
FCALLSCFUN0 (INT, namespaceGetActive, NAMESPACEGETACTIVE, namespacegetactive)
FCALLSCSUB1 (namespaceDelete, NAMESPACEDELETE, namespacedelete, INT)

/*  CDI converter routines  */


/*  parameter  */

FCALLSCSUB3 (cdiParamToString, CDIPARAMTOSTRING, cdiparamtostring, INT, PSTRING, INT)
FCALLSCSUB4 (cdiDecodeParam, CDIDECODEPARAM, cdidecodeparam, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeParam, CDIENCODEPARAM, cdiencodeparam, INT, INT, INT)

/*  date format:  YYYYMMDD  */


/*  time format:    hhmmss  */

FCALLSCSUB4 (cdiDecodeDate, CDIDECODEDATE, cdidecodedate, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeDate, CDIENCODEDATE, cdiencodedate, INT, INT, INT)
FCALLSCSUB4 (cdiDecodeTime, CDIDECODETIME, cdidecodetime, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeTime, CDIENCODETIME, cdiencodetime, INT, INT, INT)

/*  STREAM control routines  */

FCALLSCFUN2 (INT, cdiGetFiletype, CDIGETFILETYPE, cdigetfiletype, STRING, PINT)
FCALLSCFUN1 (INT, streamOpenRead, STREAMOPENREAD, streamopenread, STRING)
FCALLSCFUN2 (INT, streamOpenWrite, STREAMOPENWRITE, streamopenwrite, STRING, INT)
FCALLSCFUN1 (INT, streamOpenAppend, STREAMOPENAPPEND, streamopenappend, STRING)
FCALLSCSUB1 (streamClose, STREAMCLOSE, streamclose, INT)
FCALLSCSUB1 (streamSync, STREAMSYNC, streamsync, INT)
FCALLSCSUB2 (streamDefMaxSteps, STREAMDEFMAXSTEPS, streamdefmaxsteps, INT, INT)
FCALLSCSUB2 (streamDefNumWorker, STREAMDEFNUMWORKER, streamdefnumworker, INT, INT)
FCALLSCFUN1 (INT, streamInqNumSteps, STREAMINQNUMSTEPS, streaminqnumsteps, INT)
FCALLSCSUB2 (streamDefVlist, STREAMDEFVLIST, streamdefvlist, INT, INT)
FCALLSCFUN1 (INT, streamInqVlist, STREAMINQVLIST, streaminqvlist, INT)
FCALLSCFUN1 (INT, streamInqFiletype, STREAMINQFILETYPE, streaminqfiletype, INT)
FCALLSCSUB2 (streamDefByteorder, STREAMDEFBYTEORDER, streamdefbyteorder, INT, INT)
FCALLSCFUN1 (INT, streamInqByteorder, STREAMINQBYTEORDER, streaminqbyteorder, INT)
FCALLSCSUB2 (streamDefShuffle, STREAMDEFSHUFFLE, streamdefshuffle, INT, INT)
FCALLSCSUB4 (streamDefFilter, STREAMDEFFILTER, streamdeffilter, INT, INT, INT, PINT)
FCALLSCSUB2 (streamDefCompType, STREAMDEFCOMPTYPE, streamdefcomptype, INT, INT)
FCALLSCFUN1 (INT, streamInqCompType, STREAMINQCOMPTYPE, streaminqcomptype, INT)
FCALLSCSUB2 (streamDefCompLevel, STREAMDEFCOMPLEVEL, streamdefcomplevel, INT, INT)
FCALLSCFUN1 (INT, streamInqCompLevel, STREAMINQCOMPLEVEL, streaminqcomplevel, INT)
FCALLSCFUN2 (INT, streamDefTimestep, STREAMDEFTIMESTEP, streamdeftimestep, INT, INT)
FCALLSCFUN2 (INT, streamInqTimestep, STREAMINQTIMESTEP, streaminqtimestep, INT, INT)
FCALLSCFUN1 (INT, streamInqCurTimestepID, STREAMINQCURTIMESTEPID, streaminqcurtimestepid, INT)
FCALLSCFUN1 (STRING, streamFilename, STREAMFILENAME, streamfilename, INT)
FCALLSCFUN1 (STRING, streamFilesuffix, STREAMFILESUFFIX, streamfilesuffix, INT)
static int streamNvals_fwrap(int streamID)
{
  SizeType v;
  v = streamNvals(streamID);
  return SizeType_c2f(v);
}
FCALLSCFUN1 (INT, streamNvals_fwrap, STREAMNVALS, streamnvals, INT)
FCALLSCFUN1 (INT, streamInqNvars, STREAMINQNVARS, streaminqnvars, INT)

/*  STREAM var I/O routines (random access)  */

static void streamWriteVar_fwrap(int streamID, int varID, const double data[], int numMissVals)
{
  streamWriteVar(streamID, varID, data, (SizeType)numMissVals);
}
FCALLSCSUB4 (streamWriteVar_fwrap, STREAMWRITEVAR, streamwritevar, INT, INT, DOUBLEV, INT)
static void streamWriteVarF_fwrap(int streamID, int varID, const float data[], int numMissVals)
{
  streamWriteVarF(streamID, varID, data, (SizeType)numMissVals);
}
FCALLSCSUB4 (streamWriteVarF_fwrap, STREAMWRITEVARF, streamwritevarf, INT, INT, FLOATV, INT)
static void streamReadVar_fwrap(int streamID, int varID, double data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadVar(streamID, varID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB4 (streamReadVar_fwrap, STREAMREADVAR, streamreadvar, INT, INT, DOUBLEV, PINT)
static void streamReadVarF_fwrap(int streamID, int varID, float data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadVarF(streamID, varID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB4 (streamReadVarF_fwrap, STREAMREADVARF, streamreadvarf, INT, INT, FLOATV, PINT)
static void streamWriteVarSlice_fwrap(int streamID, int varID, int levelID, const double data[], int numMissVals)
{
  streamWriteVarSlice(streamID, varID, levelID, data, (SizeType)numMissVals);
}
FCALLSCSUB5 (streamWriteVarSlice_fwrap, STREAMWRITEVARSLICE, streamwritevarslice, INT, INT, INT, DOUBLEV, INT)
static void streamWriteVarSliceF_fwrap(int streamID, int varID, int levelID, const float data[], int numMissVals)
{
  streamWriteVarSliceF(streamID, varID, levelID, data, (SizeType)numMissVals);
}
FCALLSCSUB5 (streamWriteVarSliceF_fwrap, STREAMWRITEVARSLICEF, streamwritevarslicef, INT, INT, INT, FLOATV, INT)
static void streamReadVarSlice_fwrap(int streamID, int varID, int levelID, double data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadVarSlice(streamID, varID, levelID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB5 (streamReadVarSlice_fwrap, STREAMREADVARSLICE, streamreadvarslice, INT, INT, INT, DOUBLEV, PINT)
static void streamReadVarSliceF_fwrap(int streamID, int varID, int levelID, float data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadVarSliceF(streamID, varID, levelID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB5 (streamReadVarSliceF_fwrap, STREAMREADVARSLICEF, streamreadvarslicef, INT, INT, INT, FLOATV, PINT)
static void streamWriteVarChunk_fwrap(int streamID, int varID, const int rect[][2], const double data[], int numMissVals)
{
  streamWriteVarChunk(streamID, varID, rect, data, (SizeType)numMissVals);
}
FCALLSCSUB5 (streamWriteVarChunk_fwrap, STREAMWRITEVARCHUNK, streamwritevarchunk, INT, INT, INTVV, DOUBLEV, INT)
static void streamWriteVarChunkF_fwrap(int streamID, int varID, const int rect[][2], const float data[], int numMissVals)
{
  streamWriteVarChunkF(streamID, varID, rect, data, (SizeType)numMissVals);
}
FCALLSCSUB5 (streamWriteVarChunkF_fwrap, STREAMWRITEVARCHUNKF, streamwritevarchunkf, INT, INT, INTVV, FLOATV, INT)

/*  STREAM record I/O routines (sequential access)  */

FCALLSCSUB3 (streamDefRecord, STREAMDEFRECORD, streamdefrecord, INT, INT, INT)
FCALLSCSUB3 (streamInqRecord, STREAMINQRECORD, streaminqrecord, INT, PINT, PINT)
static void streamWriteRecord_fwrap(int streamID, const double data[], int numMissVals)
{
  streamWriteRecord(streamID, data, (SizeType)numMissVals);
}
FCALLSCSUB3 (streamWriteRecord_fwrap, STREAMWRITERECORD, streamwriterecord, INT, DOUBLEV, INT)
static void streamWriteRecordF_fwrap(int streamID, const float data[], int numMissVals)
{
  streamWriteRecordF(streamID, data, (SizeType)numMissVals);
}
FCALLSCSUB3 (streamWriteRecordF_fwrap, STREAMWRITERECORDF, streamwriterecordf, INT, FLOATV, INT)
static void streamReadRecord_fwrap(int streamID, double data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadRecord(streamID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB3 (streamReadRecord_fwrap, STREAMREADRECORD, streamreadrecord, INT, DOUBLEV, PINT)
static void streamReadRecordF_fwrap(int streamID, float data[], int *numMissVals)
{
  SizeType numMissVals_SizeType;
  streamReadRecordF(streamID, data, &numMissVals_SizeType);
  assert(numMissVals_SizeType < INT_MAX);
  *numMissVals = numMissVals_SizeType;
}
FCALLSCSUB3 (streamReadRecordF_fwrap, STREAMREADRECORDF, streamreadrecordf, INT, FLOATV, PINT)
FCALLSCSUB2 (streamCopyRecord, STREAMCOPYRECORD, streamcopyrecord, INT, INT)

/*  File driven I/O (may yield better performance than using the streamXXX functions)  */


/*  Creation & Destruction  */


/*  Advancing an iterator  */


/*  Introspecting metadata  */


/*  All outXXX arguments to these functions may be NULL.  */


/*  Reading data  */


/*  TODO[NH]: Add functions to read partial fields.  */


/*  Direct access to grib fields  */


/*  Callthroughs to GRIB-API  */


/*  Convenience functions for accessing GRIB-API keys  */


/*  VLIST routines  */

FCALLSCFUN0 (INT, vlistCreate, VLISTCREATE, vlistcreate)
FCALLSCSUB1 (vlistDestroy, VLISTDESTROY, vlistdestroy, INT)
FCALLSCFUN1 (INT, vlistDuplicate, VLISTDUPLICATE, vlistduplicate, INT)
FCALLSCSUB2 (vlistCopy, VLISTCOPY, vlistcopy, INT, INT)
FCALLSCSUB2 (vlistCopyFlag, VLISTCOPYFLAG, vlistcopyflag, INT, INT)
FCALLSCSUB1 (vlistClearFlag, VLISTCLEARFLAG, vlistclearflag, INT)
FCALLSCSUB2 (vlistCat, VLISTCAT, vlistcat, INT, INT)
FCALLSCSUB2 (vlistMerge, VLISTMERGE, vlistmerge, INT, INT)
FCALLSCSUB1 (vlistPrint, VLISTPRINT, vlistprint, INT)
FCALLSCFUN1 (INT, vlistNumber, VLISTNUMBER, vlistnumber, INT)
FCALLSCFUN1 (INT, vlistNvars, VLISTNVARS, vlistnvars, INT)
FCALLSCFUN1 (INT, vlistNgrids, VLISTNGRIDS, vlistngrids, INT)
FCALLSCFUN1 (INT, vlistNzaxis, VLISTNZAXIS, vlistnzaxis, INT)
FCALLSCFUN1 (INT, vlistNsubtypes, VLISTNSUBTYPES, vlistnsubtypes, INT)
FCALLSCSUB2 (vlistDefNtsteps, VLISTDEFNTSTEPS, vlistdefntsteps, INT, INT)
FCALLSCFUN1 (INT, vlistNtsteps, VLISTNTSTEPS, vlistntsteps, INT)
static int vlistGridsizeMax_fwrap(int vlistID)
{
  SizeType v;
  v = vlistGridsizeMax(vlistID);
  return SizeType_c2f(v);
}
FCALLSCFUN1 (INT, vlistGridsizeMax_fwrap, VLISTGRIDSIZEMAX, vlistgridsizemax, INT)
FCALLSCFUN2 (INT, vlistGrid, VLISTGRID, vlistgrid, INT, INT)
FCALLSCFUN2 (INT, vlistGridIndex, VLISTGRIDINDEX, vlistgridindex, INT, INT)
FCALLSCSUB3 (vlistChangeGridIndex, VLISTCHANGEGRIDINDEX, vlistchangegridindex, INT, INT, INT)
FCALLSCSUB3 (vlistChangeGrid, VLISTCHANGEGRID, vlistchangegrid, INT, INT, INT)
FCALLSCFUN2 (INT, vlistZaxis, VLISTZAXIS, vlistzaxis, INT, INT)
FCALLSCFUN2 (INT, vlistZaxisIndex, VLISTZAXISINDEX, vlistzaxisindex, INT, INT)
FCALLSCSUB3 (vlistChangeZaxisIndex, VLISTCHANGEZAXISINDEX, vlistchangezaxisindex, INT, INT, INT)
FCALLSCSUB3 (vlistChangeZaxis, VLISTCHANGEZAXIS, vlistchangezaxis, INT, INT, INT)
FCALLSCFUN1 (INT, vlistNrecs, VLISTNRECS, vlistnrecs, INT)
FCALLSCFUN2 (INT, vlistSubtype, VLISTSUBTYPE, vlistsubtype, INT, INT)
FCALLSCFUN2 (INT, vlistSubtypeIndex, VLISTSUBTYPEINDEX, vlistsubtypeindex, INT, INT)
FCALLSCSUB2 (vlistDefTaxis, VLISTDEFTAXIS, vlistdeftaxis, INT, INT)
FCALLSCFUN1 (INT, vlistInqTaxis, VLISTINQTAXIS, vlistinqtaxis, INT)
FCALLSCSUB2 (vlistDefTable, VLISTDEFTABLE, vlistdeftable, INT, INT)
FCALLSCFUN1 (INT, vlistInqTable, VLISTINQTABLE, vlistinqtable, INT)
FCALLSCSUB2 (vlistDefInstitut, VLISTDEFINSTITUT, vlistdefinstitut, INT, INT)
FCALLSCFUN1 (INT, vlistInqInstitut, VLISTINQINSTITUT, vlistinqinstitut, INT)
FCALLSCSUB2 (vlistDefModel, VLISTDEFMODEL, vlistdefmodel, INT, INT)
FCALLSCFUN1 (INT, vlistInqModel, VLISTINQMODEL, vlistinqmodel, INT)

/*  VLIST VAR routines  */

FCALLSCFUN5 (INT, vlistDefVarTiles, VLISTDEFVARTILES, vlistdefvartiles, INT, INT, INT, INT, INT)
FCALLSCFUN4 (INT, vlistDefVar, VLISTDEFVAR, vlistdefvar, INT, INT, INT, INT)
FCALLSCSUB3 (vlistChangeVarGrid, VLISTCHANGEVARGRID, vlistchangevargrid, INT, INT, INT)
FCALLSCSUB3 (vlistChangeVarZaxis, VLISTCHANGEVARZAXIS, vlistchangevarzaxis, INT, INT, INT)
FCALLSCSUB5 (vlistInqVar, VLISTINQVAR, vlistinqvar, INT, INT, PINT, PINT, PINT)
FCALLSCFUN2 (INT, vlistInqVarGrid, VLISTINQVARGRID, vlistinqvargrid, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarZaxis, VLISTINQVARZAXIS, vlistinqvarzaxis, INT, INT)

/*  used in MPIOM  */

FCALLSCFUN2 (INT, vlistInqVarID, VLISTINQVARID, vlistinqvarid, INT, INT)
FCALLSCSUB3 (vlistDefVarTimetype, VLISTDEFVARTIMETYPE, vlistdefvartimetype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTimetype, VLISTINQVARTIMETYPE, vlistinqvartimetype, INT, INT)
FCALLSCSUB3 (vlistDefVarTsteptype, VLISTDEFVARTSTEPTYPE, vlistdefvartsteptype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTsteptype, VLISTINQVARTSTEPTYPE, vlistinqvartsteptype, INT, INT)
FCALLSCSUB3 (vlistDefVarCompType, VLISTDEFVARCOMPTYPE, vlistdefvarcomptype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCompType, VLISTINQVARCOMPTYPE, vlistinqvarcomptype, INT, INT)
FCALLSCSUB3 (vlistDefVarCompLevel, VLISTDEFVARCOMPLEVEL, vlistdefvarcomplevel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCompLevel, VLISTINQVARCOMPLEVEL, vlistinqvarcomplevel, INT, INT)
FCALLSCSUB3 (vlistDefVarParam, VLISTDEFVARPARAM, vlistdefvarparam, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarParam, VLISTINQVARPARAM, vlistinqvarparam, INT, INT)
FCALLSCSUB3 (vlistDefVarCode, VLISTDEFVARCODE, vlistdefvarcode, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCode, VLISTINQVARCODE, vlistinqvarcode, INT, INT)
FCALLSCSUB3 (vlistDefVarDatatype, VLISTDEFVARDATATYPE, vlistdefvardatatype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarDatatype, VLISTINQVARDATATYPE, vlistinqvardatatype, INT, INT)
FCALLSCSUB3 (vlistDefVarXYZ, VLISTDEFVARXYZ, vlistdefvarxyz, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarXYZ, VLISTINQVARXYZ, vlistinqvarxyz, INT, INT)
FCALLSCSUB3 (vlistDefVarNSB, VLISTDEFVARNSB, vlistdefvarnsb, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarNSB, VLISTINQVARNSB, vlistinqvarnsb, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarNumber, VLISTINQVARNUMBER, vlistinqvarnumber, INT, INT)
FCALLSCSUB3 (vlistDefVarInstitut, VLISTDEFVARINSTITUT, vlistdefvarinstitut, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarInstitut, VLISTINQVARINSTITUT, vlistinqvarinstitut, INT, INT)
FCALLSCSUB3 (vlistDefVarModel, VLISTDEFVARMODEL, vlistdefvarmodel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarModel, VLISTINQVARMODEL, vlistinqvarmodel, INT, INT)
FCALLSCSUB3 (vlistDefVarTable, VLISTDEFVARTABLE, vlistdefvartable, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTable, VLISTINQVARTABLE, vlistinqvartable, INT, INT)
FCALLSCSUB3 (vlistDefVarName, VLISTDEFVARNAME, vlistdefvarname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarName, VLISTINQVARNAME, vlistinqvarname, INT, INT, PSTRING)
FCALLSCFUN2 (STRING, vlistCopyVarName, VLISTCOPYVARNAME, vlistcopyvarname, INT, INT)
FCALLSCSUB3 (vlistDefVarStdname, VLISTDEFVARSTDNAME, vlistdefvarstdname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarStdname, VLISTINQVARSTDNAME, vlistinqvarstdname, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarLongname, VLISTDEFVARLONGNAME, vlistdefvarlongname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarLongname, VLISTINQVARLONGNAME, vlistinqvarlongname, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarUnits, VLISTDEFVARUNITS, vlistdefvarunits, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarUnits, VLISTINQVARUNITS, vlistinqvarunits, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarMissval, VLISTDEFVARMISSVAL, vlistdefvarmissval, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, vlistInqVarMissval, VLISTINQVARMISSVAL, vlistinqvarmissval, INT, INT)
static int vlistInqVarSize_fwrap(int vlistID, int varID)
{
  SizeType v;
  v = vlistInqVarSize(vlistID, varID);
  return SizeType_c2f(v);
}
FCALLSCFUN2 (INT, vlistInqVarSize_fwrap, VLISTINQVARSIZE, vlistinqvarsize, INT, INT)
FCALLSCSUB4 (vlistDefIndex, VLISTDEFINDEX, vlistdefindex, INT, INT, INT, INT)
FCALLSCFUN3 (INT, vlistInqIndex, VLISTINQINDEX, vlistinqindex, INT, INT, INT)
FCALLSCSUB4 (vlistDefFlag, VLISTDEFFLAG, vlistdefflag, INT, INT, INT, INT)
FCALLSCFUN3 (INT, vlistInqFlag, VLISTINQFLAG, vlistinqflag, INT, INT, INT)
FCALLSCFUN2 (INT, vlistFindVar, VLISTFINDVAR, vlistfindvar, INT, INT)
FCALLSCFUN3 (INT, vlistFindLevel, VLISTFINDLEVEL, vlistfindlevel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistMergedVar, VLISTMERGEDVAR, vlistmergedvar, INT, INT)
FCALLSCFUN3 (INT, vlistMergedLevel, VLISTMERGEDLEVEL, vlistmergedlevel, INT, INT, INT)
FCALLSCSUB0 (cdiClearAdditionalKeys, CDICLEARADDITIONALKEYS, cdiclearadditionalkeys)
FCALLSCSUB1 (cdiDefAdditionalKey, CDIDEFADDITIONALKEY, cdidefadditionalkey, STRING)
FCALLSCSUB4 (vlistDefVarIntKey, VLISTDEFVARINTKEY, vlistdefvarintkey, INT, INT, STRING, INT)
FCALLSCSUB4 (vlistDefVarDblKey, VLISTDEFVARDBLKEY, vlistdefvardblkey, INT, INT, STRING, DOUBLE)
FCALLSCFUN3 (INT, vlistHasVarKey, VLISTHASVARKEY, vlisthasvarkey, INT, INT, STRING)
FCALLSCFUN3 (DOUBLE, vlistInqVarDblKey, VLISTINQVARDBLKEY, vlistinqvardblkey, INT, INT, STRING)
FCALLSCFUN3 (INT, vlistInqVarIntKey, VLISTINQVARINTKEY, vlistinqvarintkey, INT, INT, STRING)

/*  CDI attributes  */

FCALLSCFUN3 (INT, cdiInqNatts, CDIINQNATTS, cdiinqnatts, INT, INT, PINT)
FCALLSCFUN6 (INT, cdiInqAtt, CDIINQATT, cdiinqatt, INT, INT, INT, PSTRING, PINT, PINT)
FCALLSCFUN3 (INT, cdiInqAttLen, CDIINQATTLEN, cdiinqattlen, INT, INT, STRING)
FCALLSCFUN3 (INT, cdiInqAttType, CDIINQATTTYPE, cdiinqatttype, INT, INT, STRING)
FCALLSCFUN3 (INT, cdiDelAtt, CDIDELATT, cdidelatt, INT, INT, STRING)
FCALLSCFUN4 (INT, cdiCopyAtts, CDICOPYATTS, cdicopyatts, INT, INT, INT, INT)
FCALLSCFUN6 (INT, cdiDefAttInt, CDIDEFATTINT, cdidefattint, INT, INT, STRING, INT, INT, INTV)
FCALLSCFUN6 (INT, cdiDefAttFlt, CDIDEFATTFLT, cdidefattflt, INT, INT, STRING, INT, INT, DOUBLEV)
FCALLSCFUN5 (INT, cdiDefAttTxt, CDIDEFATTTXT, cdidefatttxt, INT, INT, STRING, INT, PPSTRING)
FCALLSCFUN5 (INT, cdiInqAttInt, CDIINQATTINT, cdiinqattint, INT, INT, STRING, INT, INTV)
FCALLSCFUN5 (INT, cdiInqAttFlt, CDIINQATTFLT, cdiinqattflt, INT, INT, STRING, INT, DOUBLEV)
FCALLSCFUN5 (INT, cdiInqAttTxt, CDIINQATTTXT, cdiinqatttxt, INT, INT, STRING, INT, PPSTRING)

/*  GRID routines  */

FCALLSCSUB2 (gridName, GRIDNAME, gridname, INT, PSTRING)
FCALLSCFUN1 (STRING, gridNamePtr, GRIDNAMEPTR, gridnameptr, INT)
FCALLSCSUB1 (gridCompress, GRIDCOMPRESS, gridcompress, INT)
FCALLSCSUB2 (gridDefMaskGME, GRIDDEFMASKGME, griddefmaskgme, INT, INTV)
FCALLSCFUN2 (INT, gridInqMaskGME, GRIDINQMASKGME, gridinqmaskgme, INT, INTV)
FCALLSCSUB2 (gridDefMask, GRIDDEFMASK, griddefmask, INT, INTV)
FCALLSCFUN2 (INT, gridInqMask, GRIDINQMASK, gridinqmask, INT, INTV)
static int gridCreate_fwrap(int gridtype, int size)
{
  int v;
  v = gridCreate(gridtype, (SizeType)size);
  return v;
}
FCALLSCFUN2 (INT, gridCreate_fwrap, GRIDCREATE, gridcreate, INT, INT)
FCALLSCSUB1 (gridDestroy, GRIDDESTROY, griddestroy, INT)
FCALLSCFUN1 (INT, gridDuplicate, GRIDDUPLICATE, gridduplicate, INT)
FCALLSCSUB2 (gridDefProj, GRIDDEFPROJ, griddefproj, INT, INT)
FCALLSCFUN1 (INT, gridInqProj, GRIDINQPROJ, gridinqproj, INT)
FCALLSCFUN1 (INT, gridInqProjType, GRIDINQPROJTYPE, gridinqprojtype, INT)
FCALLSCFUN1 (INT, gridInqType, GRIDINQTYPE, gridinqtype, INT)
static int gridInqSize_fwrap(int gridID)
{
  SizeType v;
  v = gridInqSize(gridID);
  return SizeType_c2f(v);
}
FCALLSCFUN1 (INT, gridInqSize_fwrap, GRIDINQSIZE, gridinqsize, INT)
static void gridDefXsize_fwrap(int gridID, int xsize)
{
  gridDefXsize(gridID, (SizeType)xsize);
}
FCALLSCSUB2 (gridDefXsize_fwrap, GRIDDEFXSIZE, griddefxsize, INT, INT)
static int gridInqXsize_fwrap(int gridID)
{
  SizeType v;
  v = gridInqXsize(gridID);
  return SizeType_c2f(v);
}
FCALLSCFUN1 (INT, gridInqXsize_fwrap, GRIDINQXSIZE, gridinqxsize, INT)
static void gridDefYsize_fwrap(int gridID, int ysize)
{
  gridDefYsize(gridID, (SizeType)ysize);
}
FCALLSCSUB2 (gridDefYsize_fwrap, GRIDDEFYSIZE, griddefysize, INT, INT)
static int gridInqYsize_fwrap(int gridID)
{
  SizeType v;
  v = gridInqYsize(gridID);
  return SizeType_c2f(v);
}
FCALLSCFUN1 (INT, gridInqYsize_fwrap, GRIDINQYSIZE, gridinqysize, INT)
FCALLSCSUB2 (gridDefNP, GRIDDEFNP, griddefnp, INT, INT)
FCALLSCFUN1 (INT, gridInqNP, GRIDINQNP, gridinqnp, INT)
FCALLSCSUB2 (gridDefXvals, GRIDDEFXVALS, griddefxvals, INT, DOUBLEV)
static int gridInqXvals_fwrap(int gridID, double xvals[])
{
  SizeType v;
  v = gridInqXvals(gridID, xvals);
  return SizeType_c2f(v);
}
FCALLSCFUN2 (INT, gridInqXvals_fwrap, GRIDINQXVALS, gridinqxvals, INT, DOUBLEV)
static int gridInqXvalsPart_fwrap(int gridID, int start, int size, double xvals[])
{
  SizeType v;
  v = gridInqXvalsPart(gridID, start, (SizeType)size, xvals);
  return SizeType_c2f(v);
}
FCALLSCFUN4 (INT, gridInqXvalsPart_fwrap, GRIDINQXVALSPART, gridinqxvalspart, INT, INT, INT, DOUBLEV)
FCALLSCFUN1 (INT, gridInqXIsc, GRIDINQXISC, gridinqxisc, INT)
FCALLSCSUB2 (gridDefYvals, GRIDDEFYVALS, griddefyvals, INT, DOUBLEV)
static int gridInqYvals_fwrap(int gridID, double yvals[])
{
  SizeType v;
  v = gridInqYvals(gridID, yvals);
  return SizeType_c2f(v);
}
FCALLSCFUN2 (INT, gridInqYvals_fwrap, GRIDINQYVALS, gridinqyvals, INT, DOUBLEV)
static int gridInqYvalsPart_fwrap(int gridID, int start, int size, double yvals[])
{
  SizeType v;
  v = gridInqYvalsPart(gridID, start, (SizeType)size, yvals);
  return SizeType_c2f(v);
}
FCALLSCFUN4 (INT, gridInqYvalsPart_fwrap, GRIDINQYVALSPART, gridinqyvalspart, INT, INT, INT, DOUBLEV)
FCALLSCFUN1 (INT, gridInqYIsc, GRIDINQYISC, gridinqyisc, INT)

/*  CDI var keys  */


/*  String keys  */


/*  Integer keys  */


/*  Floating point keys  */


/*  Byte array keys  */

FCALLSCFUN4 (INT, cdiDefKeyInt, CDIDEFKEYINT, cdidefkeyint, INT, INT, INT, INT)
FCALLSCFUN4 (INT, cdiInqKeyInt, CDIINQKEYINT, cdiinqkeyint, INT, INT, INT, PINT)
FCALLSCFUN4 (INT, cdiDefKeyFloat, CDIDEFKEYFLOAT, cdidefkeyfloat, INT, INT, INT, DOUBLE)

/*  cdiInqKeyFloat Get a float value from a key  */

FCALLSCFUN4 (INT, cdiInqKeyFloat, CDIINQKEYFLOAT, cdiinqkeyfloat, INT, INT, INT, PDOUBLE)
FCALLSCFUN4 (INT, cdiDefKeyString, CDIDEFKEYSTRING, cdidefkeystring, INT, INT, INT, STRING)
FCALLSCFUN5 (INT, cdiInqKeyString, CDIINQKEYSTRING, cdiinqkeystring, INT, INT, INT, PSTRING, PINT)
FCALLSCFUN4 (INT, cdiInqKeyLen, CDIINQKEYLEN, cdiinqkeylen, INT, INT, INT, PINT)
FCALLSCFUN4 (INT, cdiCopyKeys, CDICOPYKEYS, cdicopykeys, INT, INT, INT, INT)
FCALLSCFUN4 (INT, cdiCopyKey, CDICOPYKEY, cdicopykey, INT, INT, INT, INT)
FCALLSCFUN3 (INT, cdiDeleteKey, CDIDELETEKEY, cdideletekey, INT, INT, INT)

/*  GRID routines  */

FCALLSCSUB2 (gridDefXname, GRIDDEFXNAME, griddefxname, INT, STRING)
FCALLSCSUB2 (gridInqXname, GRIDINQXNAME, gridinqxname, INT, PSTRING)
FCALLSCSUB2 (gridDefXlongname, GRIDDEFXLONGNAME, griddefxlongname, INT, STRING)
FCALLSCSUB2 (gridInqXlongname, GRIDINQXLONGNAME, gridinqxlongname, INT, PSTRING)
FCALLSCSUB2 (gridDefXunits, GRIDDEFXUNITS, griddefxunits, INT, STRING)
FCALLSCSUB2 (gridInqXunits, GRIDINQXUNITS, gridinqxunits, INT, PSTRING)
FCALLSCSUB2 (gridDefYname, GRIDDEFYNAME, griddefyname, INT, STRING)
FCALLSCSUB2 (gridInqYname, GRIDINQYNAME, gridinqyname, INT, PSTRING)
FCALLSCSUB2 (gridDefYlongname, GRIDDEFYLONGNAME, griddefylongname, INT, STRING)
FCALLSCSUB2 (gridInqYlongname, GRIDINQYLONGNAME, gridinqylongname, INT, PSTRING)
FCALLSCSUB2 (gridDefYunits, GRIDDEFYUNITS, griddefyunits, INT, STRING)
FCALLSCSUB2 (gridInqYunits, GRIDINQYUNITS, gridinqyunits, INT, PSTRING)
FCALLSCSUB2 (gridDefDatatype, GRIDDEFDATATYPE, griddefdatatype, INT, INT)
FCALLSCFUN1 (INT, gridInqDatatype, GRIDINQDATATYPE, gridinqdatatype, INT)
static double gridInqXval_fwrap(int gridID, int index)
{
  double v;
  v = gridInqXval(gridID, (SizeType)index);
  return v;
}
FCALLSCFUN2 (DOUBLE, gridInqXval_fwrap, GRIDINQXVAL, gridinqxval, INT, INT)
static double gridInqYval_fwrap(int gridID, int index)
{
  double v;
  v = gridInqYval(gridID, (SizeType)index);
  return v;
}
FCALLSCFUN2 (DOUBLE, gridInqYval_fwrap, GRIDINQYVAL, gridinqyval, INT, INT)
FCALLSCFUN1 (DOUBLE, gridInqXinc, GRIDINQXINC, gridinqxinc, INT)
FCALLSCFUN1 (DOUBLE, gridInqYinc, GRIDINQYINC, gridinqyinc, INT)
FCALLSCFUN1 (INT, gridIsCircular, GRIDISCIRCULAR, gridiscircular, INT)
FCALLSCFUN1 (INT, gridInqTrunc, GRIDINQTRUNC, gridinqtrunc, INT)
FCALLSCSUB2 (gridDefTrunc, GRIDDEFTRUNC, griddeftrunc, INT, INT)

/*  Reference of an unstructured grid  */

FCALLSCSUB2 (gridDefNumber, GRIDDEFNUMBER, griddefnumber, INT, INT)
FCALLSCFUN1 (INT, gridInqNumber, GRIDINQNUMBER, gridinqnumber, INT)
FCALLSCSUB2 (gridDefPosition, GRIDDEFPOSITION, griddefposition, INT, INT)
FCALLSCFUN1 (INT, gridInqPosition, GRIDINQPOSITION, gridinqposition, INT)
FCALLSCSUB2 (gridDefReference, GRIDDEFREFERENCE, griddefreference, INT, STRING)
FCALLSCFUN2 (INT, gridInqReference, GRIDINQREFERENCE, gridinqreference, INT, PSTRING)
FCALLSCSUB2 (gridDefUUID, GRIDDEFUUID, griddefuuid, INT, PVOID)
FCALLSCSUB2 (gridInqUUID, GRIDINQUUID, gridinquuid, INT, PVOID)

/*  Rotated Lon/Lat grid  */

FCALLSCSUB4 (gridDefParamRLL, GRIDDEFPARAMRLL, griddefparamrll, INT, DOUBLE, DOUBLE, DOUBLE)
FCALLSCSUB4 (gridInqParamRLL, GRIDINQPARAMRLL, gridinqparamrll, INT, PDOUBLE, PDOUBLE, PDOUBLE)

/*  Hexagonal GME grid  */

FCALLSCSUB5 (gridDefParamGME, GRIDDEFPARAMGME, griddefparamgme, INT, INT, INT, INT, INT)
FCALLSCSUB5 (gridInqParamGME, GRIDINQPARAMGME, gridinqparamgme, INT, PINT, PINT, PINT, PINT)
FCALLSCSUB2 (gridDefArea, GRIDDEFAREA, griddefarea, INT, DOUBLEV)
FCALLSCSUB2 (gridInqArea, GRIDINQAREA, gridinqarea, INT, DOUBLEV)
FCALLSCFUN1 (INT, gridHasArea, GRIDHASAREA, gridhasarea, INT)
FCALLSCSUB2 (gridDefNvertex, GRIDDEFNVERTEX, griddefnvertex, INT, INT)
FCALLSCFUN1 (INT, gridInqNvertex, GRIDINQNVERTEX, gridinqnvertex, INT)
FCALLSCSUB2 (gridDefXbounds, GRIDDEFXBOUNDS, griddefxbounds, INT, DOUBLEV)
static int gridInqXbounds_fwrap(int gridID, double xbounds[])
{
  SizeType v;
  v = gridInqXbounds(gridID, xbounds);
  return SizeType_c2f(v);
}
FCALLSCFUN2 (INT, gridInqXbounds_fwrap, GRIDINQXBOUNDS, gridinqxbounds, INT, DOUBLEV)
static int gridInqXboundsPart_fwrap(int gridID, int start, int size, double xbounds[])
{
  SizeType v;
  v = gridInqXboundsPart(gridID, start, (SizeType)size, xbounds);
  return SizeType_c2f(v);
}
FCALLSCFUN4 (INT, gridInqXboundsPart_fwrap, GRIDINQXBOUNDSPART, gridinqxboundspart, INT, INT, INT, DOUBLEV)
FCALLSCSUB2 (gridDefYbounds, GRIDDEFYBOUNDS, griddefybounds, INT, DOUBLEV)
static int gridInqYbounds_fwrap(int gridID, double ybounds[])
{
  SizeType v;
  v = gridInqYbounds(gridID, ybounds);
  return SizeType_c2f(v);
}
FCALLSCFUN2 (INT, gridInqYbounds_fwrap, GRIDINQYBOUNDS, gridinqybounds, INT, DOUBLEV)
static int gridInqYboundsPart_fwrap(int gridID, int start, int size, double ybounds[])
{
  SizeType v;
  v = gridInqYboundsPart(gridID, start, (SizeType)size, ybounds);
  return SizeType_c2f(v);
}
FCALLSCFUN4 (INT, gridInqYboundsPart_fwrap, GRIDINQYBOUNDSPART, gridinqyboundspart, INT, INT, INT, DOUBLEV)
FCALLSCSUB3 (gridDefReducedPoints, GRIDDEFREDUCEDPOINTS, griddefreducedpoints, INT, INT, INTV)
FCALLSCSUB2 (gridInqReducedPoints, GRIDINQREDUCEDPOINTS, gridinqreducedpoints, INT, INTV)
FCALLSCSUB2 (gridChangeType, GRIDCHANGETYPE, gridchangetype, INT, INT)
FCALLSCSUB2 (gridDefComplexPacking, GRIDDEFCOMPLEXPACKING, griddefcomplexpacking, INT, INT)
FCALLSCFUN1 (INT, gridInqComplexPacking, GRIDINQCOMPLEXPACKING, gridinqcomplexpacking, INT)

/*  ZAXIS routines  */

FCALLSCSUB2 (zaxisName, ZAXISNAME, zaxisname, INT, PSTRING)
FCALLSCFUN1 (STRING, zaxisNamePtr, ZAXISNAMEPTR, zaxisnameptr, INT)
FCALLSCFUN2 (INT, zaxisCreate, ZAXISCREATE, zaxiscreate, INT, INT)
FCALLSCSUB1 (zaxisDestroy, ZAXISDESTROY, zaxisdestroy, INT)
FCALLSCFUN1 (INT, zaxisInqType, ZAXISINQTYPE, zaxisinqtype, INT)
FCALLSCFUN1 (INT, zaxisInqSize, ZAXISINQSIZE, zaxisinqsize, INT)
FCALLSCFUN1 (INT, zaxisDuplicate, ZAXISDUPLICATE, zaxisduplicate, INT)
FCALLSCSUB2 (zaxisDefLevels, ZAXISDEFLEVELS, zaxisdeflevels, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqLevels, ZAXISINQLEVELS, zaxisinqlevels, INT, DOUBLEV)
FCALLSCFUN1 (INT, zaxisInqCLen, ZAXISINQCLEN, zaxisinqclen, INT)
FCALLSCSUB3 (zaxisDefLevel, ZAXISDEFLEVEL, zaxisdeflevel, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, zaxisInqLevel, ZAXISINQLEVEL, zaxisinqlevel, INT, INT)
FCALLSCSUB2 (zaxisDefNlevRef, ZAXISDEFNLEVREF, zaxisdefnlevref, INT, INT)
FCALLSCFUN1 (INT, zaxisInqNlevRef, ZAXISINQNLEVREF, zaxisinqnlevref, INT)
FCALLSCSUB2 (zaxisDefNumber, ZAXISDEFNUMBER, zaxisdefnumber, INT, INT)
FCALLSCFUN1 (INT, zaxisInqNumber, ZAXISINQNUMBER, zaxisinqnumber, INT)
FCALLSCSUB2 (zaxisDefUUID, ZAXISDEFUUID, zaxisdefuuid, INT, PVOID)
FCALLSCSUB2 (zaxisInqUUID, ZAXISINQUUID, zaxisinquuid, INT, PVOID)
FCALLSCSUB2 (zaxisDefName, ZAXISDEFNAME, zaxisdefname, INT, STRING)
FCALLSCSUB2 (zaxisInqName, ZAXISINQNAME, zaxisinqname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefLongname, ZAXISDEFLONGNAME, zaxisdeflongname, INT, STRING)
FCALLSCSUB2 (zaxisInqLongname, ZAXISINQLONGNAME, zaxisinqlongname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefUnits, ZAXISDEFUNITS, zaxisdefunits, INT, STRING)
FCALLSCSUB2 (zaxisInqUnits, ZAXISINQUNITS, zaxisinqunits, INT, PSTRING)
FCALLSCSUB2 (zaxisInqStdname, ZAXISINQSTDNAME, zaxisinqstdname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefDatatype, ZAXISDEFDATATYPE, zaxisdefdatatype, INT, INT)
FCALLSCFUN1 (INT, zaxisInqDatatype, ZAXISINQDATATYPE, zaxisinqdatatype, INT)
FCALLSCSUB2 (zaxisDefPositive, ZAXISDEFPOSITIVE, zaxisdefpositive, INT, INT)
FCALLSCFUN1 (INT, zaxisInqPositive, ZAXISINQPOSITIVE, zaxisinqpositive, INT)
FCALLSCSUB1 (zaxisDefScalar, ZAXISDEFSCALAR, zaxisdefscalar, INT)
FCALLSCFUN1 (INT, zaxisInqScalar, ZAXISINQSCALAR, zaxisinqscalar, INT)
FCALLSCSUB3 (zaxisDefVct, ZAXISDEFVCT, zaxisdefvct, INT, INT, DOUBLEV)
FCALLSCSUB2 (zaxisInqVct, ZAXISINQVCT, zaxisinqvct, INT, DOUBLEV)
FCALLSCFUN1 (INT, zaxisInqVctSize, ZAXISINQVCTSIZE, zaxisinqvctsize, INT)
FCALLSCSUB2 (zaxisDefLbounds, ZAXISDEFLBOUNDS, zaxisdeflbounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqLbounds, ZAXISINQLBOUNDS, zaxisinqlbounds, INT, DOUBLEV)
FCALLSCFUN2 (DOUBLE, zaxisInqLbound, ZAXISINQLBOUND, zaxisinqlbound, INT, INT)
FCALLSCSUB2 (zaxisDefUbounds, ZAXISDEFUBOUNDS, zaxisdefubounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqUbounds, ZAXISINQUBOUNDS, zaxisinqubounds, INT, DOUBLEV)
FCALLSCFUN2 (DOUBLE, zaxisInqUbound, ZAXISINQUBOUND, zaxisinqubound, INT, INT)
FCALLSCSUB2 (zaxisDefWeights, ZAXISDEFWEIGHTS, zaxisdefweights, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqWeights, ZAXISINQWEIGHTS, zaxisinqweights, INT, DOUBLEV)
FCALLSCSUB2 (zaxisChangeType, ZAXISCHANGETYPE, zaxischangetype, INT, INT)

/*  TAXIS routines  */

FCALLSCFUN1 (INT, taxisCreate, TAXISCREATE, taxiscreate, INT)
FCALLSCSUB1 (taxisDestroy, TAXISDESTROY, taxisdestroy, INT)
FCALLSCFUN1 (INT, taxisDuplicate, TAXISDUPLICATE, taxisduplicate, INT)
FCALLSCSUB2 (taxisCopyTimestep, TAXISCOPYTIMESTEP, taxiscopytimestep, INT, INT)
FCALLSCSUB2 (taxisDefType, TAXISDEFTYPE, taxisdeftype, INT, INT)
FCALLSCFUN1 (INT, taxisInqType, TAXISINQTYPE, taxisinqtype, INT)
FCALLSCSUB2 (taxisDefVdate, TAXISDEFVDATE, taxisdefvdate, INT, INT)
FCALLSCSUB2 (taxisDefVtime, TAXISDEFVTIME, taxisdefvtime, INT, INT)
FCALLSCFUN1 (INT, taxisInqVdate, TAXISINQVDATE, taxisinqvdate, INT)
FCALLSCFUN1 (INT, taxisInqVtime, TAXISINQVTIME, taxisinqvtime, INT)
FCALLSCSUB2 (taxisDefRdate, TAXISDEFRDATE, taxisdefrdate, INT, INT)
FCALLSCSUB2 (taxisDefRtime, TAXISDEFRTIME, taxisdefrtime, INT, INT)
FCALLSCFUN1 (INT, taxisInqRdate, TAXISINQRDATE, taxisinqrdate, INT)
FCALLSCFUN1 (INT, taxisInqRtime, TAXISINQRTIME, taxisinqrtime, INT)
FCALLSCFUN1 (INT, taxisHasBounds, TAXISHASBOUNDS, taxishasbounds, INT)
FCALLSCSUB1 (taxisWithBounds, TAXISWITHBOUNDS, taxiswithbounds, INT)
FCALLSCSUB1 (taxisDeleteBounds, TAXISDELETEBOUNDS, taxisdeletebounds, INT)
FCALLSCSUB3 (taxisDefVdateBounds, TAXISDEFVDATEBOUNDS, taxisdefvdatebounds, INT, INT, INT)
FCALLSCSUB3 (taxisDefVtimeBounds, TAXISDEFVTIMEBOUNDS, taxisdefvtimebounds, INT, INT, INT)
FCALLSCSUB3 (taxisInqVdateBounds, TAXISINQVDATEBOUNDS, taxisinqvdatebounds, INT, PINT, PINT)
FCALLSCSUB3 (taxisInqVtimeBounds, TAXISINQVTIMEBOUNDS, taxisinqvtimebounds, INT, PINT, PINT)
FCALLSCSUB2 (taxisDefCalendar, TAXISDEFCALENDAR, taxisdefcalendar, INT, INT)
FCALLSCFUN1 (INT, taxisInqCalendar, TAXISINQCALENDAR, taxisinqcalendar, INT)
FCALLSCSUB2 (taxisDefTunit, TAXISDEFTUNIT, taxisdeftunit, INT, INT)
FCALLSCFUN1 (INT, taxisInqTunit, TAXISINQTUNIT, taxisinqtunit, INT)
FCALLSCSUB2 (taxisDefForecastTunit, TAXISDEFFORECASTTUNIT, taxisdefforecasttunit, INT, INT)
FCALLSCFUN1 (INT, taxisInqForecastTunit, TAXISINQFORECASTTUNIT, taxisinqforecasttunit, INT)
FCALLSCSUB2 (taxisDefForecastPeriod, TAXISDEFFORECASTPERIOD, taxisdefforecastperiod, INT, DOUBLE)
FCALLSCFUN1 (DOUBLE, taxisInqForecastPeriod, TAXISINQFORECASTPERIOD, taxisinqforecastperiod, INT)
FCALLSCSUB2 (taxisDefNumavg, TAXISDEFNUMAVG, taxisdefnumavg, INT, INT)
FCALLSCFUN1 (INT, taxisInqNumavg, TAXISINQNUMAVG, taxisinqnumavg, INT)
FCALLSCFUN1 (STRING, taxisNamePtr, TAXISNAMEPTR, taxisnameptr, INT)
FCALLSCFUN1 (STRING, tunitNamePtr, TUNITNAMEPTR, tunitnameptr, INT)

/*  Institut routines  */

FCALLSCFUN4 (INT, institutDef, INSTITUTDEF, institutdef, INT, INT, STRING, STRING)
FCALLSCFUN4 (INT, institutInq, INSTITUTINQ, institutinq, INT, INT, STRING, STRING)
FCALLSCFUN0 (INT, institutInqNumber, INSTITUTINQNUMBER, institutinqnumber)
FCALLSCFUN1 (INT, institutInqCenter, INSTITUTINQCENTER, institutinqcenter, INT)
FCALLSCFUN1 (INT, institutInqSubcenter, INSTITUTINQSUBCENTER, institutinqsubcenter, INT)
FCALLSCFUN1 (STRING, institutInqNamePtr, INSTITUTINQNAMEPTR, institutinqnameptr, INT)
FCALLSCFUN1 (STRING, institutInqLongnamePtr, INSTITUTINQLONGNAMEPTR, institutinqlongnameptr, INT)

/*  Model routines  */

FCALLSCFUN3 (INT, modelDef, MODELDEF, modeldef, INT, INT, STRING)
FCALLSCFUN3 (INT, modelInq, MODELINQ, modelinq, INT, INT, STRING)
FCALLSCFUN1 (INT, modelInqInstitut, MODELINQINSTITUT, modelinqinstitut, INT)
FCALLSCFUN1 (INT, modelInqGribID, MODELINQGRIBID, modelinqgribid, INT)
FCALLSCFUN1 (STRING, modelInqNamePtr, MODELINQNAMEPTR, modelinqnameptr, INT)

/*  Table routines  */

FCALLSCSUB2 (tableWrite, TABLEWRITE, tablewrite, STRING, INT)
FCALLSCFUN1 (INT, tableRead, TABLEREAD, tableread, STRING)
FCALLSCFUN3 (INT, tableDef, TABLEDEF, tabledef, INT, INT, STRING)
FCALLSCFUN1 (STRING, tableInqNamePtr, TABLEINQNAMEPTR, tableinqnameptr, INT)
FCALLSCFUN3 (INT, tableInq, TABLEINQ, tableinq, INT, INT, STRING)
FCALLSCFUN0 (INT, tableInqNumber, TABLEINQNUMBER, tableinqnumber)
FCALLSCFUN1 (INT, tableInqNum, TABLEINQNUM, tableinqnum, INT)
FCALLSCFUN1 (INT, tableInqModel, TABLEINQMODEL, tableinqmodel, INT)
FCALLSCSUB6 (tableInqEntry, TABLEINQENTRY, tableinqentry, INT, INT, INT, PSTRING, PSTRING, PSTRING)

/*  Subtype routines  */

FCALLSCFUN1 (INT, subtypeCreate, SUBTYPECREATE, subtypecreate, INT)

/*  Gives a textual summary of the variable subtype  */

FCALLSCSUB1 (subtypePrint, SUBTYPEPRINT, subtypeprint, INT)

/*  Compares two subtype data structures  */

FCALLSCFUN2 (INT, subtypeCompare, SUBTYPECOMPARE, subtypecompare, INT, INT)
FCALLSCFUN1 (INT, subtypeInqSize, SUBTYPEINQSIZE, subtypeinqsize, INT)
FCALLSCFUN1 (INT, subtypeInqActiveIndex, SUBTYPEINQACTIVEINDEX, subtypeinqactiveindex, INT)
FCALLSCSUB2 (subtypeDefActiveIndex, SUBTYPEDEFACTIVEINDEX, subtypedefactiveindex, INT, INT)

/*  Generate a "query object" out of a key-value pair  */


/*  Generate an AND-combined "query object" out of two previous query objects  */

FCALLSCFUN3 (INT, subtypeInqTile, SUBTYPEINQTILE, subtypeinqtile, INT, INT, INT)
FCALLSCFUN4 (INT, subtypeInqAttribute, SUBTYPEINQATTRIBUTE, subtypeinqattribute, INT, INT, STRING, PINT)
FCALLSCFUN2 (INT, vlistInqVarSubtype, VLISTINQVARSUBTYPE, vlistinqvarsubtype, INT, INT)
FCALLSCSUB3 (gribapiLibraryVersion, GRIBAPILIBRARYVERSION, gribapilibraryversion, PINT, PINT, PINT)

/*  Compatibility functions for release 1.8.3 (obsolete functions)  */

FCALLSCSUB2 (zaxisDefLtype, ZAXISDEFLTYPE, zaxisdefltype, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTypeOfGeneratingProcess, VLISTINQVARTYPEOFGENERATINGPROCESS, vlistinqvartypeofgeneratingprocess, INT, INT)
FCALLSCSUB3 (vlistDefVarTypeOfGeneratingProcess, VLISTDEFVARTYPEOFGENERATINGPROCESS, vlistdefvartypeofgeneratingprocess, INT, INT, INT)
FCALLSCSUB3 (vlistDefVarProductDefinitionTemplate, VLISTDEFVARPRODUCTDEFINITIONTEMPLATE, vlistdefvarproductdefinitiontemplate, INT, INT, INT)

/*  End of fortran interface  */


#if defined __clang__
#  pragma GCC diagnostic pop
#endif

// clang-format on

#endif
