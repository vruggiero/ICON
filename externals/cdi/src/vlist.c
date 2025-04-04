#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "error.h"
#include "vlist.h"
#include "zaxis.h"
#include "varscan.h"
#include "namespace.h"
#include "resource_handle.h"
#include "vlist_var.h"
#include "cdi_att.h"

#include "resource_unpack.h"
#include "serialize.h"

#ifdef HAVE_LIBGRIB_API
/* list of additional GRIB2 keywords which are read by the open process */
int cdiNAdditionalGRIBKeys = 0;
char *cdiAdditionalGRIBKeys[MAX_OPT_GRIB_ENTRIES];
#endif

static int VLIST_Debug = 0;

static void vlist_initialize(void);

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>

static pthread_once_t _vlist_init_thread = PTHREAD_ONCE_INIT;

#define VLIST_INIT() pthread_once(&_vlist_init_thread, vlist_initialize)

#else

static bool vlistIsInitialized = false;

#define VLIST_INIT() \
  if (!vlistIsInitialized) vlist_initialize()
#endif

static int
vlist_compare(vlist_t *a, vlist_t *b)
{
  int diff = (a->nvars != b->nvars) | (a->ngrids != b->ngrids) | (a->nzaxis != b->nzaxis) | (a->instID != b->instID)
             | (a->modelID != b->modelID) | (a->tableID != b->tableID) | (a->ntsteps != b->ntsteps)
             | (a->atts.nelems != b->atts.nelems);

  int nvars = a->nvars;
  for (int varID = 0; varID < nvars; ++varID) diff |= vlistVarCompare(a, varID, b, varID);

  size_t natts = a->atts.nelems;
  for (size_t attID = 0; attID < natts; ++attID) diff |= cdi_att_compare(&a->atts, &a->atts, (int) attID);

  return diff;
}

static void vlistPrintKernel(vlist_t *vlistptr, FILE *fp);
static void vlist_delete(vlist_t *vlistptr);

static int vlistGetSizeP(void *vlistptr, void *context);
static void vlistPackP(void *vlistptr, void *buff, int size, int *position, void *context);
static int vlistTxCode(void *vlistptr);

#if !defined(__cplusplus)
const
#endif
    resOps vlistOps
    = { (valCompareFunc) vlist_compare,
        (valDestroyFunc) vlist_delete,
        (valPrintFunc) vlistPrintKernel,
        vlistGetSizeP,
        vlistPackP,
        vlistTxCode };

vlist_t *
vlist_to_pointer(int vlistID)
{
  VLIST_INIT();
  return (vlist_t *) reshGetVal(vlistID, &vlistOps);
}

static void
vlist_init_entry(vlist_t *vlistptr)
{
  vlistptr->immutable = 0;
  vlistptr->internal = 0;
  vlistptr->self = CDI_UNDEFID;
  vlistptr->nvars = 0;
  vlistptr->vars = NULL;
  vlistptr->ngrids = 0;
  vlistptr->nzaxis = 0;
  vlistptr->taxisID = CDI_UNDEFID;
  vlistptr->instID = CDI_Default_InstID;
  vlistptr->modelID = CDI_Default_ModelID;
  vlistptr->tableID = CDI_Default_TableID;
  vlistptr->varsAllocated = 0;
  vlistptr->ntsteps = CDI_UNDEFID;
  vlistptr->keys.nalloc = MAX_KEYS;
  vlistptr->keys.nelems = 0;
  for (int i = 0; i < MAX_KEYS; ++i) vlistptr->keys.value[i].length = 0;
  vlistptr->atts.nalloc = MAX_ATTRIBUTES;
  vlistptr->atts.nelems = 0;
  vlistptr->nsubtypes = 0;
  for (int i = 0; i < MAX_SUBTYPES_PS; ++i) vlistptr->subtypeIDs[i] = CDI_UNDEFID;
}

static vlist_t *
vlist_new_entry(cdiResH resH)
{
  vlist_t *vlistptr = (vlist_t *) Malloc(sizeof(vlist_t));
  vlist_init_entry(vlistptr);
  if (resH == CDI_UNDEFID)
    vlistptr->self = reshPut(vlistptr, &vlistOps);
  else
    {
      vlistptr->self = resH;
      reshReplace(resH, vlistptr, &vlistOps);
    }
  return vlistptr;
}

static void
vlist_delete_entry(int vlistID)
{
  reshRemove(vlistID, &vlistOps);

  if (VLIST_Debug) Message("Removed idx %d from vlist list", vlistID);
}

static void
vlist_initialize(void)
{
  char *env = getenv("VLIST_DEBUG");
  if (env) VLIST_Debug = atoi(env);
#ifndef HAVE_LIBPTHREAD
  vlistIsInitialized = true;
#endif
}

static void
vlist_copy(vlist_t *vlistptr2, vlist_t *vlistptr1)
{
  int vlistID2 = vlistptr2->self;
  int vlist2internal = vlistptr2->internal;
  memcpy(vlistptr2, vlistptr1, sizeof(vlist_t));
  vlistptr2->internal = vlist2internal;  // the question who's responsible to destroy the vlist is tied to its containing memory
                                         // region, so we retain this flag
  vlistptr2->immutable = 0;              // this is a copy, so it's mutable, independent of whether the original is mutable or not
  vlistptr2->keys.nelems = 0;
  vlistptr2->atts.nelems = 0;
  vlistptr2->self = vlistID2;
}

void
cdiVlistMakeInternal(int vlistID)
{
  vlist_to_pointer(vlistID)->internal = 1;
}

void
cdiVlistMakeImmutable(int vlistID)
{
  vlist_to_pointer(vlistID)->immutable = 1;
}

/*
@Function  vlistCreate
@Title     Create a variable list

@Prototype int vlistCreate(void)

@Example
Here is an example using @func{vlistCreate} to create a variable list
and add a variable with @func{vlistDefVar}.

@Source
#include "cdi.h"
   ...
int vlistID, varID;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
   ...
streamDefVlist(streamID, vlistID);
   ...
vlistDestroy(vlistID);
   ...
@EndSource
@EndFunction
*/
int
vlistCreate(void)
{
  cdiInitialize();

  VLIST_INIT();

  vlist_t *vlistptr = vlist_new_entry(CDI_UNDEFID);
  if (CDI_Debug) Message("create vlistID = %d", vlistptr->self);
  return vlistptr->self;
}

static void
vlist_delete(vlist_t *vlistptr)
{
  int vlistID = vlistptr->self;
  if (CDI_Debug) Message("call to vlist_delete, vlistID = %d", vlistID);

  cdiDeleteKeys(vlistID, CDI_GLOBAL);
  cdiDeleteAtts(vlistID, CDI_GLOBAL);

  int nvars = vlistptr->nvars;
  var_t *vars = vlistptr->vars;

  for (int varID = 0; varID < nvars; varID++)
    {
      if (vars[varID].levinfo) Free(vars[varID].levinfo);

      if (vlistptr->vars[varID].opt_grib_kvpair)
        {
          for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; i++)
            {
              if (vlistptr->vars[varID].opt_grib_kvpair[i].keyword) Free(vlistptr->vars[varID].opt_grib_kvpair[i].keyword);
            }
          Free(vlistptr->vars[varID].opt_grib_kvpair);
        }
      vlistptr->vars[varID].opt_grib_nentries = 0;
      vlistptr->vars[varID].opt_grib_kvpair_size = 0;
      vlistptr->vars[varID].opt_grib_kvpair = NULL;

      cdiDeleteKeys(vlistID, varID);
      cdiDeleteAtts(vlistID, varID);
    }

  if (vars) Free(vars);

  Free(vlistptr);
}

// destroy a vlist object, should always be called through namespace lookup
void
cdiVlistDestroy_(int vlistID, bool assertInternal)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  static const char warningTextUserByInternal[]
      = "Destroying a vlist object that is owned by the user (vlistID=%d).\n"
        "This is most likely because of a missing vlistDestroy() in the application code.\n"
        "If that's not the case, and you are absolutely certain about it, please report the bug.",
      warningTextInternalByUser[] = "Attempt to destroy an internal vlist object by the user (vlistID=%d).";
  static const char *const wText[2] = { warningTextUserByInternal, warningTextInternalByUser };
  if (vlistptr->internal == assertInternal)
    {
      vlist_delete(vlistptr);
      vlist_delete_entry(vlistID);
    }
  else
    Warning(wText[!assertInternal], vlistID);
}

/*
@Function  vlistDestroy
@Title     Destroy a variable list

@Prototype void vlistDestroy(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.

@EndFunction
*/
void
vlistDestroy(int vlistID)
{
  void (*mycdiVlistDestroy_)(int, bool) = (void (*)(int, bool)) namespaceSwitchGet(NSSWITCH_VLIST_DESTROY_).func;
  mycdiVlistDestroy_(vlistID, false);
}

static void
var_copy_entries(var_t *var2, var_t *var1)
{
  var2->opt_grib_kvpair_size = 0;
  var2->opt_grib_kvpair = NULL;
  var2->opt_grib_nentries = 0;

  resize_opt_grib_entries(var2, var1->opt_grib_nentries);
  var2->opt_grib_nentries = var1->opt_grib_nentries;
  if ((var2->opt_grib_nentries > 0) && CDI_Debug) Message("copy %d optional GRIB keywords", var2->opt_grib_nentries);

  for (int i = 0; i < var1->opt_grib_nentries; i++)
    {
      if (CDI_Debug) Message("copy entry \"%s\" ...", var1->opt_grib_kvpair[i].keyword);
      var2->opt_grib_kvpair[i].keyword = NULL;
      if (var1->opt_grib_kvpair[i].keyword != NULL)
        {
          var2->opt_grib_kvpair[i] = var1->opt_grib_kvpair[i];
          var2->opt_grib_kvpair[i].keyword = strdup(var1->opt_grib_kvpair[i].keyword);
          var2->opt_grib_kvpair[i].update = true;
          if (CDI_Debug) Message("done.");
        }
      else
        {
          if (CDI_Debug) Message("not done.");
        }
    }
}

/*
@Function  vlistCopy
@Title     Copy a variable list

@Prototype void vlistCopy(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
The function @func{vlistCopy} copies all entries from vlistID1 to vlistID2.

@EndFunction
*/
void
vlistCopy(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1);
  vlist_t *vlistptr2 = vlist_to_pointer(vlistID2);
  if (CDI_Debug) Message("call to vlistCopy, vlistIDs %d -> %d", vlistID1, vlistID2);

  var_t *vars1 = vlistptr1->vars;
  var_t *vars2 = vlistptr2->vars;
  vlist_copy(vlistptr2, vlistptr1);

  vlistptr2->keys.nelems = 0;
  cdiCopyKeys(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);
  vlistptr2->atts.nelems = 0;
  cdiCopyAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if (vars1)
    {
      int nvars = vlistptr1->nvars;
      // vlistptr2->varsAllocated = nvars;

      size_t n = (size_t) vlistptr2->varsAllocated;
      vars2 = (var_t *) Realloc(vars2, n * sizeof(var_t));
      memcpy(vars2, vars1, n * sizeof(var_t));
      vlistptr2->vars = vars2;

      for (int varID = 0; varID < nvars; varID++)
        {
          var_copy_entries(&vars2[varID], &vars1[varID]);
          vlistptr2->vars[varID].keys.nelems = 0;
          cdiCopyKeys(vlistID1, varID, vlistID2, varID);

          vlistptr2->vars[varID].atts.nelems = 0;
          cdiCopyAtts(vlistID1, varID, vlistID2, varID);

          if (vars1[varID].levinfo)
            {
              n = (size_t) zaxisInqSize(vars1[varID].zaxisID);
              vars2[varID].levinfo = (levinfo_t *) Malloc(n * sizeof(levinfo_t));
              memcpy(vars2[varID].levinfo, vars1[varID].levinfo, n * sizeof(levinfo_t));
            }
        }
    }
}

/*
@Function  vlistDuplicate
@Title     Duplicate a variable list

@Prototype int vlistDuplicate(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistDuplicate} duplicates the variable list from vlistID1.

@Result
@func{vlistDuplicate} returns an identifier to the duplicated variable list.

@EndFunction
*/
int
vlistDuplicate(int vlistID)
{
  if (CDI_Debug) Message("call to vlistDuplicate");

  int vlistIDnew = vlistCreate();
  vlistCopy(vlistIDnew, vlistID);
  return vlistIDnew;
}

void
vlistClearFlag(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int varID = 0; varID < vlistptr->nvars; varID++)
    {
      vlistptr->vars[varID].flag = false;
      if (vlistptr->vars[varID].levinfo)
        {
          int nlevs = zaxisInqSize(vlistptr->vars[varID].zaxisID);
          for (int levID = 0; levID < nlevs; levID++) vlistptr->vars[varID].levinfo[levID].flag = false;
        }
    }
}

struct vgzSearchState
{
  int resIDValue;
  int zaxistype;
  int nlevels;
  const double *levels;
  const double *lbounds;
  const double *ubounds;
};

static enum cdiApplyRet
vgzZAxisSearch(int id, void *res, void *data)
{
  struct vgzSearchState *state = (struct vgzSearchState *) data;
  (void) res;
  if (zaxis_compare(id, state->zaxistype, state->nlevels, state->levels, state->lbounds, state->ubounds, NULL, NULL, 0, -1)
      == false)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

static int
vlist_generate_zaxis(int vlistID, int zaxistype, int nlevels, const double *levels, const double *lbounds, const double *ubounds,
                     int vctsize, const double *vct, const char **cvals, size_t clen)
{
  int zaxisID = CDI_UNDEFID;
  bool zaxisdefined = false;
  bool zaxisglobdefined = false;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int nzaxis = vlistptr->nzaxis;

  bool hasBounds = (lbounds && ubounds);

  for (int index = 0; index < nzaxis; ++index)
    {
      zaxisID = vlistptr->zaxisIDs[index];

      if (zaxis_compare(zaxisID, zaxistype, nlevels, levels, lbounds, ubounds, NULL, NULL, 0, -1) == false)
        {
          zaxisdefined = true;
          break;
        }
    }

  if (!zaxisdefined)
    {
      struct vgzSearchState query;
      query.zaxistype = zaxistype;
      query.nlevels = nlevels;
      query.levels = levels;
      query.lbounds = lbounds;
      query.ubounds = ubounds;

      if ((zaxisglobdefined = (cdiResHFilterApply(getZaxisOps(), vgzZAxisSearch, &query) == CDI_APPLY_STOP)))
        zaxisID = query.resIDValue;
    }

  if (!zaxisdefined)
    {
      if (!zaxisglobdefined)
        {
          zaxisID = zaxisCreate(zaxistype, nlevels);
          zaxisDefLevels(zaxisID, levels);

          if (zaxistype == ZAXIS_CHAR) zaxisDefCvals(zaxisID, cvals, (int) clen);

          if (hasBounds)
            {
              zaxisDefLbounds(zaxisID, lbounds);
              zaxisDefUbounds(zaxisID, ubounds);
            }

          if (zaxistype == ZAXIS_HYBRID && vctsize > 0) zaxisDefVct(zaxisID, vctsize, vct);
        }

      nzaxis = vlistptr->nzaxis;
      vlistptr->zaxisIDs[nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  return zaxisID;
}

/*
@Function  vlistCopyFlag
@Title     Copy some entries of a variable list

@Prototype void vlistCopyFlag(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
The function @func{vlistCopyFlag} copies all entries with a flag from vlistID1 to vlistID2.

@EndFunction
*/
void
vlistCopyFlag(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1);
  vlist_t *vlistptr2 = vlist_to_pointer(vlistID2);
  var_t *vars1 = vlistptr1->vars;
  var_t *vars2 = vlistptr2->vars;

  vlist_copy(vlistptr2, vlistptr1);

  vlistptr2->keys.nelems = 0;
  cdiCopyKeys(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);
  vlistptr2->atts.nelems = 0;
  cdiCopyAtts(vlistID1, CDI_GLOBAL, vlistID2, CDI_GLOBAL);

  if (vlistptr1->vars)
    {
      vlistptr2->ngrids = 0;
      vlistptr2->nzaxis = 0;

      int nvars = vlistptr1->nvars;
      int nvars2 = 0;
      for (int varID = 0; varID < nvars; varID++) nvars2 += vars1[varID].flag;

      vlistptr2->nvars = nvars2;
      vlistptr2->varsAllocated = nvars2;
      vars2 = (nvars2 > 0) ? (var_t *) Malloc((size_t) nvars2 * sizeof(var_t)) : NULL;

      vlistptr2->vars = vars2;

      int varID2 = 0;
      for (int varID = 0; varID < nvars; varID++)
        if (vars1[varID].flag)
          {
            vlistptr2->vars[varID2].flag = false;
            int zaxisID = vlistptr1->vars[varID].zaxisID;
            int gridID = vlistptr1->vars[varID].gridID;
            int subtypeID = vlistptr1->vars[varID].subtypeID;

            memcpy(&vars2[varID2], &vars1[varID], sizeof(var_t));

            vars1[varID].fvarID = varID2;
            vars2[varID2].fvarID = varID;

            vars2[varID2].mvarID = varID2;

            var_copy_entries(&vars2[varID2], &vars1[varID]);
            vlistptr2->vars[varID2].keys.nelems = 0;
            cdiCopyKeys(vlistID1, varID, vlistID2, varID2);

            vlistptr2->vars[varID2].atts.nelems = 0;
            cdiCopyAtts(vlistID1, varID, vlistID2, varID2);

            int nlevs = zaxisInqSize(vars1[varID].zaxisID);
            int nlevs2 = 0;
            if (vars1[varID].levinfo)
              for (int levID = 0; levID < nlevs; levID++) nlevs2 += vars1[varID].levinfo[levID].flag;

            vars2[varID2].levinfo = (levinfo_t *) Malloc((size_t) nlevs2 * sizeof(levinfo_t));

            if (nlevs != nlevs2)
              {
                int nvct = 0;
                double *levels = NULL;
                double *lbounds = NULL, *ubounds = NULL;
                const double *vct = NULL;

                if (!vars1[varID].levinfo) cdiVlistCreateVarLevInfo(vlistptr1, varID);

                zaxisID = vars1[varID].zaxisID;
                int zaxisType = zaxisInqType(zaxisID);

                int levID2 = 0;
                for (int levID = 0; levID < nlevs; levID++)
                  if (vars1[varID].levinfo[levID].flag)
                    {
                      vars1[varID].levinfo[levID].flevelID = levID2;
                      vars1[varID].levinfo[levID].mlevelID = levID2;
                    }

                if (zaxisInqLevels(zaxisID, NULL))
                  {
                    levels = (double *) Malloc((size_t) nlevs2 * sizeof(double));

                    levID2 = 0;
                    for (int levID = 0; levID < nlevs; ++levID)
                      if (vars1[varID].levinfo[levID].flag) levels[levID2++] = zaxisInqLevel(zaxisID, levID);
                  }

                if (zaxisType == ZAXIS_HYBRID)
                  {
                    nvct = zaxisInqVctSize(zaxisID);
                    vct = zaxisInqVctPtr(zaxisID);
                  }

                size_t clen2 = 0;
                char **cvals2 = NULL;
#ifndef USE_MPI
                if (zaxisType == ZAXIS_CHAR)
                  {
                    char **cvals1 = zaxisInqCValsPtr(zaxisID);
                    size_t clen1 = (size_t) zaxisInqCLen(zaxisID);
                    for (int levID = 0; levID < nlevs; ++levID)
                      if (vars1[varID].levinfo[levID].flag)
                        {
                          size_t testlen = clen1;
                          while (cvals1[levID][testlen] == ' ') testlen--;
                          if (clen2 < testlen) clen2 = testlen;
                        }
                    cvals2 = (char **) Malloc((size_t) nlevs2 * sizeof(char *));
                    levID2 = 0;

                    for (int levID = 0; levID < nlevs; ++levID)
                      if (vars1[varID].levinfo[levID].flag)
                        {
                          cvals2[levID2] = (char *) Malloc((size_t) (clen2) * sizeof(char));
                          memcpy(cvals2[levID2], cvals1[levID], clen2 * sizeof(char));
                          levID2++;
                        }
                  }
#endif

                if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
                  {
                    lbounds = (double *) Malloc(2 * (size_t) nlevs2 * sizeof(double));
                    ubounds = lbounds + nlevs2;

                    double *lbounds1 = (double *) Malloc(2 * (size_t) nlevs * sizeof(double)), *ubounds1 = lbounds1 + nlevs;

                    zaxisInqLbounds(zaxisID, lbounds1);
                    zaxisInqUbounds(zaxisID, ubounds1);

                    levID2 = 0;
                    for (int levID = 0; levID < nlevs; ++levID)
                      if (vars1[varID].levinfo[levID].flag)
                        {
                          lbounds[levID2] = lbounds1[levID];
                          ubounds[levID2] = ubounds1[levID];
                          levID2++;
                        }

                    Free(lbounds1);
                  }

                int zaxisID2 = vlist_generate_zaxis(vlistID2, zaxisType, nlevs2, levels, lbounds, ubounds, nvct, vct,
                                                    (const char **) cvals2, clen2);
                if (levels) Free(levels);
                if (lbounds) Free(lbounds);
                if (cvals2)
                  {
                    for (int levID = 0; levID < nlevs2; ++levID) Free(cvals2[levID]);
                    Free(cvals2);
                  }

                char ctemp[CDI_MAX_NAME];
                int length = CDI_MAX_NAME;
                cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, ctemp, &length);
                cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_NAME, ctemp);
                length = CDI_MAX_NAME;
                cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, ctemp, &length);
                cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_LONGNAME, ctemp);
                length = CDI_MAX_NAME;
                cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, ctemp, &length);
                cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_UNITS, ctemp);

                zaxisDefDatatype(zaxisID2, zaxisInqDatatype(zaxisID));
                zaxisDefPositive(zaxisID2, zaxisInqPositive(zaxisID));

                if (zaxisType == ZAXIS_CHAR)
                  {
                    char dimname[CDI_MAX_NAME + 3];
                    length = sizeof(dimname);
                    cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_DIMNAME, dimname, &length);
                    if (dimname[0] == 0)
                      {
                        memcpy(dimname, "area_type", 10);
                        dimname[10] = 0;
                      }
                    cdiDefKeyString(zaxisID2, CDI_GLOBAL, CDI_KEY_DIMNAME, dimname);
                  }

                if (zaxisType == ZAXIS_GENERIC) cdiCopyKey(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, zaxisID2);

                cdiCopyAtts(zaxisID, CDI_GLOBAL, zaxisID2, CDI_GLOBAL);

                zaxisID = zaxisID2;
                vars2[varID2].zaxisID = zaxisID2;
              }

            for (int levID = 0; levID < nlevs2; levID++)
              {
                vars2[varID2].levinfo[levID].flag = false;
                vars2[varID2].levinfo[levID].index = -1;
              }

            int levID2 = 0;
            for (int levID = 0; levID < nlevs; levID++)
              if (vars1[varID].levinfo[levID].flag)
                {
                  vars2[varID2].levinfo[levID2].flevelID = levID;
                  vars2[varID2].levinfo[levID2].mlevelID = levID2;
                  levID2++;
                }

            vlistAdd2GridIDs(vlistptr2, gridID);
            vlistAdd2ZaxisIDs(vlistptr2, zaxisID);
            vlistAdd2SubtypeIDs(vlistptr2, subtypeID);

            varID2++;
          }
    }
}

/*
@Function  vlistCat
@Title     Concatenate two variable lists

@Prototype void vlistCat(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
Concatenate the variable list vlistID1 at the end of vlistID2.

@EndFunction
*/
void
vlistCat(int vlistID2, int vlistID1)
{
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1);
  vlist_t *vlistptr2 = vlist_to_pointer(vlistID2);
  var_t *vars1 = vlistptr1->vars;
  var_t *vars2 = vlistptr2->vars;
  int nvars1 = vlistptr1->nvars;
  int nvars2 = vlistptr2->nvars;
  int nvars = nvars1 + nvars2;
  vlistptr2->nvars = nvars;

  if (nvars > vlistptr2->varsAllocated)
    {
      vlistptr2->varsAllocated = nvars;
      vars2 = (var_t *) Realloc(vars2, (size_t) nvars * sizeof(var_t));
      vlistptr2->vars = vars2;
    }
  memcpy(vars2 + nvars2, vars1, (size_t) nvars1 * sizeof(var_t));

  for (int varID = 0; varID < nvars1; varID++)
    {
      int varID2 = varID + nvars2;
      vars1[varID].fvarID = varID2;
      vars2[varID2].fvarID = varID;

      vars1[varID].mvarID = varID2;
      vars2[varID2].mvarID = varID;

      if (vars1[varID].param < 0)
        {
          int pnum, pcat, pdis;
          cdiDecodeParam(vars1[varID].param, &pnum, &pcat, &pdis);
          pnum = -(varID2 + 1);
          vars2[varID2].param = cdiEncodeParam(pnum, pcat, pdis);
        }

      var_copy_entries(&vars2[varID2], &vars1[varID]);
      vars2[varID2].keys.nelems = 0;
      cdiCopyKeys(vlistID1, varID, vlistID2, varID2);

      if (vars1[varID].levinfo)
        {
          size_t nlevs = (size_t) zaxisInqSize(vars1[varID].zaxisID);
          vars2[varID2].levinfo = (levinfo_t *) Malloc(nlevs * sizeof(levinfo_t));
          memcpy(vars2[varID2].levinfo, vars1[varID].levinfo, nlevs * sizeof(levinfo_t));
        }

      vars2[varID2].atts.nelems = 0;
      cdiCopyAtts(vlistID1, varID, vlistID2, varID2);

      vlistAdd2GridIDs(vlistptr2, vars1[varID].gridID);
      vlistAdd2ZaxisIDs(vlistptr2, vars1[varID].zaxisID);
      vlistAdd2SubtypeIDs(vlistptr2, vars1[varID].subtypeID);
    }
}

/*
@Function  vlistMerge
@Title     Merge two variable lists

@Prototype void vlistMerge(int vlistID2, int vlistID1)
@Parameter
    @Item  vlistID2  Target variable list ID.
    @Item  vlistID1  Source variable list ID.

@Description
Merge the variable list vlistID1 to the variable list vlistID2.

@EndFunction
*/
void
vlistMerge(int vlistID2, int vlistID1)
{
  int varID = 0;
  vlist_t *vlistptr1 = vlist_to_pointer(vlistID1);
  vlist_t *vlistptr2 = vlist_to_pointer(vlistID2);
  var_t *vars1 = vlistptr1->vars;
  var_t *vars2 = vlistptr2->vars;
  int nvars1 = vlistptr1->nvars;
  int nvars2 = vlistptr2->nvars;

  if (nvars1 == nvars2)
    {
      char name1[CDI_MAX_NAME], name2[CDI_MAX_NAME];
      for (varID = 0; varID < nvars2; varID++)
        {
          size_t ngp1 = gridInqSize(vars1[varID].gridID);
          size_t ngp2 = gridInqSize(vars2[varID].gridID);
          if (ngp1 != ngp2) break;

          int length = CDI_MAX_NAME;
          (void) cdiInqKeyString(vlistID1, varID, CDI_KEY_NAME, name1, &length);
          length = CDI_MAX_NAME;
          (void) cdiInqKeyString(vlistID2, varID, CDI_KEY_NAME, name2, &length);

          if (*name1 && *name2)
            {
              if (!str_is_equal(name1, name2)) break;
            }
          else
            {
              if (vars1[varID].param != vars2[varID].param) break;
            }
        }
    }

  if (varID == nvars2) /* same variables in vlistID1 and vlistID2 */
    {
      for (varID = 0; varID < nvars2; varID++)
        {
          vars1[varID].fvarID = varID;
          vars2[varID].fvarID = varID;

          vars1[varID].mvarID = varID;
          vars2[varID].mvarID = varID;

          int nlevs1 = zaxisInqSize(vars1[varID].zaxisID);
          int nlevs2 = zaxisInqSize(vars2[varID].zaxisID);

          int nlevs = nlevs1 + nlevs2;

          /*
          fprintf(stderr, "var %d %d %d %d %d\n", varID, nlevs1, nlevs2, nlevs, sizeof(levinfo_t));
          */
          if (vars1[varID].levinfo)
            {
              vars2[varID].levinfo = (levinfo_t *) Realloc(vars2[varID].levinfo, (size_t) nlevs * sizeof(levinfo_t));

              memcpy(vars2[varID].levinfo + nlevs2, vars1[varID].levinfo, (size_t) nlevs1 * sizeof(levinfo_t));
            }
          else
            cdiVlistCreateVarLevInfo(vlistptr1, varID);

          for (int levID = 0; levID < nlevs1; levID++) vars1[varID].levinfo[levID].mlevelID = nlevs2 + levID;
        }

      bool *lvar = (bool *) Calloc((size_t) nvars2, sizeof(bool));

      for (varID = 0; varID < nvars2; varID++)
        {
          if (lvar[varID] == true) continue;

          int zaxisID1 = vars1[varID].zaxisID;
          int zaxisID2 = vars2[varID].zaxisID;
          // nlevs1 = zaxisInqSize(vars1[varID].zaxisID);
          // nlevs2 = zaxisInqSize(vars2[varID].zaxisID);
          int nlevs1 = zaxisInqSize(zaxisID1);
          int nlevs2 = zaxisInqSize(zaxisID2);
          // fprintf(stderr, "zaxis %d %d %d %d\n", zaxisID1, zaxisID2, nlevs1, nlevs2);

          int nlevs = nlevs1 + nlevs2;

          int zaxisID = zaxisDuplicate(zaxisID2);
          zaxisResize(zaxisID, nlevs);

          if (zaxisInqLevels(zaxisID1, NULL))
            {
              double *levels = (double *) Malloc((size_t) nlevs1 * sizeof(double));

              zaxisInqLevels(zaxisID1, levels);
              /*
                for (int levID = 0; levID < nlevs1; levID++)
                  fprintf(stderr, "%d %d %d %d %d %g\n", varID, levID, nlevs1, nlevs2, vars2[varID].nlevs, levels[levID]);
              */
              for (int levID = 0; levID < nlevs1; levID++) zaxisDefLevel(zaxisID, nlevs2 + levID, levels[levID]);

              Free(levels);
            }

          for (int index = 0; index < vlistptr2->nzaxis; index++)
            if (vlistptr2->zaxisIDs[index] == zaxisID2) vlistptr2->zaxisIDs[index] = zaxisID;

          for (int varID2 = 0; varID2 < nvars2; varID2++)
            if (lvar[varID2] == false && vars2[varID2].zaxisID == zaxisID2)
              {
                vars2[varID2].zaxisID = zaxisID;
                lvar[varID2] = true;
              }
        }

      Free(lvar);
    }
  else
    {
      vlistCat(vlistID2, vlistID1);
    }
}

/*
@Function  vlistNvars
@Title     Number of variables in a variable list

@Prototype int vlistNvars(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNvars} returns the number of variables in the variable list vlistID.

@Result
@func{vlistNvars} returns the number of variables in a variable list.

@EndFunction
*/
int
vlistNvars(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return vlistptr->nvars;
}

int
vlistNrecs(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int nrecs = 0;
  for (int varID = 0; varID < vlistptr->nvars; varID++) nrecs += zaxisInqSize(vlistptr->vars[varID].zaxisID);

  return nrecs;
}

int
vlistNumber(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int datatype = vlistptr->vars[0].datatype;
  int number = (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) ? CDI_COMP : CDI_REAL;

  for (int varID = 1; varID < vlistptr->nvars; varID++)
    {
      datatype = vlistptr->vars[varID].datatype;
      int number2 = (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) ? CDI_COMP : CDI_REAL;
      if (number2 != number)
        {
          number = CDI_BOTH;
          break;
        }
    }

  return number;
}

/*
@Function  vlistNgrids
@Title     Number of grids in a variable list

@Prototype int vlistNgrids(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNgrids} returns the number of grids in the variable list vlistID.

@Result
@func{vlistNgrids} returns the number of grids in a variable list.

@EndFunction
*/
int
vlistNgrids(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return vlistptr->ngrids;
}

/*
@Function  vlistNzaxis
@Title     Number of zaxis in a variable list

@Prototype int vlistNzaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistNzaxis} returns the number of zaxis in the variable list vlistID.

@Result
@func{vlistNzaxis} returns the number of zaxis in a variable list.

@EndFunction
*/
int
vlistNzaxis(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return vlistptr->nzaxis;
}

int
vlistNsubtypes(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return vlistptr->nsubtypes;
}

void
vlistDefNtsteps(int vlistID, int nts)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->ntsteps != nts)
    {
      vlistptr->ntsteps = nts;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

// This function is used in CDO!
int
vlistNtsteps(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (int) vlistptr->ntsteps;
}

static void
vlistPrintKernel(vlist_t *vlistptr, FILE *fp)
{
  int vlistID = vlistptr->self;
  fprintf(fp, "#\n# vlistID %d\n#\n", vlistID);

  int nvars = vlistptr->nvars;

  fprintf(fp,
          "nvars    : %d\n"
          "ngrids   : %d\n"
          "nzaxis   : %d\n"
          "nsubtypes: %d\n"
          "taxisID  : %d\n"
          "instID   : %d\n"
          "modelID  : %d\n"
          "tableID  : %d\n",
          nvars, vlistptr->ngrids, vlistptr->nzaxis, vlistptr->nsubtypes, vlistptr->taxisID, vlistptr->instID, vlistptr->modelID,
          vlistptr->tableID);

  if (nvars > 0)
    {
      fprintf(fp, " varID param    gridID zaxisID stypeID tsteptype flag name     longname         units\n");
      for (int varID = 0; varID < nvars; varID++)
        {
          int param = vlistptr->vars[varID].param;
          int gridID = vlistptr->vars[varID].gridID;
          int zaxisID = vlistptr->vars[varID].zaxisID;
          int subtypeID = vlistptr->vars[varID].subtypeID;
          int tsteptype = vlistptr->vars[varID].tsteptype;
          char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
          int length = CDI_MAX_NAME;
          (void) cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length);
          length = CDI_MAX_NAME;
          (void) cdiInqKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname, &length);
          length = CDI_MAX_NAME;
          (void) cdiInqKeyString(vlistID, varID, CDI_KEY_UNITS, units, &length);
          int flag = vlistptr->vars[varID].flag;

          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          fprintf(fp, "%6d %-8s %6d  %6d  %6d  %6d  %5d %-8s %s [%s]\n", varID, paramstr, gridID, zaxisID, subtypeID, tsteptype,
                  flag, name, longname, units);
        }

      fputs("\n"
            " varID  levID fvarID flevID mvarID mlevID  index  dtype  flag  level\n",
            fp);
      for (int varID = 0; varID < nvars; varID++)
        {
          int zaxisID = vlistptr->vars[varID].zaxisID;
          int nlevs = zaxisInqSize(zaxisID);
          int fvarID = vlistptr->vars[varID].fvarID;
          int mvarID = vlistptr->vars[varID].mvarID;
          int dtype = vlistptr->vars[varID].datatype;
          for (int levID = 0; levID < nlevs; levID++)
            {
              levinfo_t li;
              if (vlistptr->vars[varID].levinfo)
                li = vlistptr->vars[varID].levinfo[levID];
              else
                li = DEFAULT_LEVINFO(levID);
              int flevID = li.flevelID;
              int mlevID = li.mlevelID;
              int index = li.index;
              int flag = li.flag;

              double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levID) : levID + 1;

              fprintf(fp, "%6d %6d %6d %6d %6d %6d %6d %6d %5d  %.9g\n", varID, levID, fvarID, flevID, mvarID, mlevID, index, dtype,
                      flag, level);
            }
        }

      fputs("\n"
            " varID  size\n",
            fp);
      for (int varID = 0; varID < nvars; varID++)
        fprintf(fp, "%3d %8zu\n", varID,
                (size_t) zaxisInqSize(vlistptr->vars[varID].zaxisID) * gridInqSize(vlistptr->vars[varID].gridID));
    }
}

void
vlistPrint(int vlistID)
{
  if (vlistID == CDI_UNDEFID) return;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  vlistPrintKernel(vlistptr, stdout);
}

/*
@Function  vlistDefTaxis
@Title     Define the time axis

@Prototype void vlistDefTaxis(int vlistID, int taxisID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}.

@Description
The function @func{vlistDefTaxis} defines the time axis of a variable list.

@EndFunction
*/
void
vlistDefTaxis(int vlistID, int taxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->taxisID != taxisID)
    {
      // FIXME: This code seems to leak a taxis_t object if `vlistptr->taxisID` was valid before the call to vlistDefTaxis.
      vlistptr->taxisID = taxisID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistInqTaxis
@Title     Get the time axis

@Prototype int vlistInqTaxis(int vlistID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.

@Description
The function @func{vlistInqTaxis} returns the time axis of a variable list.

@Result
@func{vlistInqTaxis} returns an identifier to the time axis.

@EndFunction
*/
int
vlistInqTaxis(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return vlistptr->taxisID;
}

void
vlistDefTable(int vlistID, int tableID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->tableID != tableID)
    {
      vlistptr->tableID = tableID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqTable(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return vlistptr->tableID;
}

void
vlistDefInstitut(int vlistID, int instID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->instID != instID)
    {
      vlistptr->instID = instID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqInstitut(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int instID = vlistptr->instID;

  if (instID == CDI_UNDEFID)
    {
      instID = vlistInqVarInstitut(vlistID, 0);

      for (int varID = 1; varID < vlistptr->nvars; varID++)
        if (instID != vlistInqVarInstitut(vlistID, varID))
          {
            instID = CDI_UNDEFID;
            break;
          }
      vlistDefInstitut(vlistID, instID);
    }

  return instID;
}

void
vlistDefModel(int vlistID, int modelID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->modelID != modelID)
    {
      vlistptr->modelID = modelID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqModel(int vlistID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int modelID = vlistptr->modelID;

  if (modelID == CDI_UNDEFID)
    {
      modelID = vlistInqVarModel(vlistID, 0);

      for (int varID = 1; varID < vlistptr->nvars; varID++)
        if (modelID != vlistInqVarModel(vlistID, varID))
          {
            modelID = CDI_UNDEFID;
            break;
          }

      vlistDefModel(vlistID, modelID);
    }

  return modelID;
}

SizeType
vlistGridsizeMax(int vlistID)
{
  SizeType gridsizemax = 0;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int index = 0; index < vlistptr->ngrids; index++)
    {
      int gridID = vlistptr->gridIDs[index];
      SizeType gridsize = gridInqSize(gridID);
      if (gridsize > gridsizemax) gridsizemax = gridsize;
    }

  return gridsizemax;
}

int
vlistGrid(int vlistID, int index)
{
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (index < vlistptr->ngrids && index >= 0) gridID = vlistptr->gridIDs[index];

  return gridID;
}

int
vlistGridIndex(int vlistID, int gridID)
{
  int index;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (index = 0; index < vlistptr->ngrids; index++)
    if (gridID == vlistptr->gridIDs[index]) break;

  if (index == vlistptr->ngrids) index = -1;

  return index;
}

void
vlistChangeGridIndex(int vlistID, int index, int gridID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int gridIDold = vlistptr->gridIDs[index];
  if (gridIDold != gridID)
    {
      vlistptr->gridIDs[index] = gridID;

      int nvars = vlistptr->nvars;
      for (int varID = 0; varID < nvars; varID++)
        if (vlistptr->vars[varID].gridID == gridIDold)
          {
            vlistptr->vars[varID].gridID = gridID;
            int chunkSize = 0;
            cdiInqKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE, &chunkSize);
            if (chunkSize > 0) cdiDeleteKey(vlistID, varID, CDI_KEY_CHUNKSIZE);
          }
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

void
vlistChangeGrid(int vlistID, int gridID1, int gridID2)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (gridID1 != gridID2)
    {
      int ngrids = vlistptr->ngrids;
      for (int index = 0; index < ngrids; index++)
        {
          if (vlistptr->gridIDs[index] == gridID1)
            {
              vlistptr->gridIDs[index] = gridID2;
              break;
            }
        }
      int nvars = vlistptr->nvars;
      for (int varID = 0; varID < nvars; varID++)
        if (vlistptr->vars[varID].gridID == gridID1) vlistptr->vars[varID].gridID = gridID2;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistZaxis(int vlistID, int index)
{
  int zaxisID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (index < vlistptr->nzaxis && index >= 0) zaxisID = vlistptr->zaxisIDs[index];

  return zaxisID;
}

int
vlistZaxisIndex(int vlistID, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int index;
  for (index = 0; index < vlistptr->nzaxis; index++)
    if (zaxisID == vlistptr->zaxisIDs[index]) break;

  if (index == vlistptr->nzaxis) index = -1;

  return index;
}

void
vlistChangeZaxisIndex(int vlistID, int index, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int zaxisIDold = vlistptr->zaxisIDs[index];
  if (zaxisIDold != zaxisID)
    {
      vlistptr->zaxisIDs[index] = zaxisID;

      int nlevs = zaxisInqSize(zaxisID), nlevsOld = zaxisInqSize(zaxisIDold);
      int nvars = vlistptr->nvars;
      for (int varID = 0; varID < nvars; varID++)
        if (vlistptr->vars[varID].zaxisID == zaxisIDold)
          {
            vlistptr->vars[varID].zaxisID = zaxisID;
            if (vlistptr->vars[varID].levinfo && nlevs != nlevsOld)
              {
                vlistptr->vars[varID].levinfo
                    = (levinfo_t *) Realloc(vlistptr->vars[varID].levinfo, (size_t) nlevs * sizeof(levinfo_t));

                for (int levID = 0; levID < nlevs; levID++) vlistptr->vars[varID].levinfo[levID] = DEFAULT_LEVINFO(levID);
              }
          }
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

void
vlistChangeZaxis(int vlistID, int zaxisID1, int zaxisID2)
{
  int nlevs1 = zaxisInqSize(zaxisID1), nlevs2 = zaxisInqSize(zaxisID2);
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int nzaxis = vlistptr->nzaxis;
  for (int index = 0; index < nzaxis; index++)
    {
      if (vlistptr->zaxisIDs[index] == zaxisID1)
        {
          vlistptr->zaxisIDs[index] = zaxisID2;
          break;
        }
    }

  int nvars = vlistptr->nvars;
  for (int varID = 0; varID < nvars; varID++)
    if (vlistptr->vars[varID].zaxisID == zaxisID1)
      {
        vlistptr->vars[varID].zaxisID = zaxisID2;

        if (vlistptr->vars[varID].levinfo && nlevs2 != nlevs1)
          {
            vlistptr->vars[varID].levinfo
                = (levinfo_t *) Realloc(vlistptr->vars[varID].levinfo, (size_t) nlevs2 * sizeof(levinfo_t));

            for (int levID = 0; levID < nlevs2; levID++) vlistptr->vars[varID].levinfo[levID] = DEFAULT_LEVINFO(levID);
          }
      }
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

int
vlistSubtype(int vlistID, int index)
{
  int subtypeID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (index < vlistptr->nsubtypes && index >= 0) subtypeID = vlistptr->subtypeIDs[index];

  return subtypeID;
}

int
vlistSubtypeIndex(int vlistID, int subtypeID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int index;
  for (index = vlistptr->nsubtypes; index--;)
    if (subtypeID == vlistptr->subtypeIDs[index]) break;

  return index;
}

int
vlistHasTime(int vlistID)
{
  bool hastime = false;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (!(CDI_Reduce_Dim && vlistptr->ntsteps == 1))
    {
      size_t nvars = vlistptr->nvars > 0 ? (size_t) vlistptr->nvars : (size_t) 0;
      var_t *restrict vars = vlistptr->vars;
      for (size_t varID = 0; varID < nvars; varID++)
        if (vars[varID].timetype != TIME_CONSTANT)
          {
            hastime = true;
            break;
          }
    }

  return (int) hastime;
}

enum
{
  VLIST_PACK_INT_SELF,
  VLIST_PACK_INT_NVARS,
  VLIST_PACK_INT_TAXISID,
  VLIST_PACK_INT_TABLEID,
  VLIST_PACK_INT_INSTID,
  VLIST_PACK_INT_MODELID,
  vlistNints,
};

static int
vlistTxCode(void *vlistptr)
{
  (void) vlistptr;
  return VLIST;
}

static int
vlistGetSizeP(void *vlistptr, void *context)
{
  vlist_t *p = (vlist_t *) vlistptr;
  int txsize = serializeGetSize(vlistNints, CDI_DATATYPE_INT, context);
  txsize += serializeGetSize(1, CDI_DATATYPE_LONG, context);
  txsize += cdiAttsGetSize(p, CDI_GLOBAL, context);
  for (int varID = 0; varID < p->nvars; varID++) txsize += vlistVarGetPackSize(p, varID, context);
  return txsize;
}

static void
vlistPackP(void *vlistptr, void *buf, int size, int *position, void *context)
{
  int tempbuf[vlistNints];
  vlist_t *p = (vlist_t *) vlistptr;
  tempbuf[VLIST_PACK_INT_SELF] = p->self;
  tempbuf[VLIST_PACK_INT_NVARS] = p->nvars;
  tempbuf[VLIST_PACK_INT_TAXISID] = p->taxisID;
  tempbuf[VLIST_PACK_INT_TABLEID] = p->tableID;
  tempbuf[VLIST_PACK_INT_INSTID] = p->instID;
  tempbuf[VLIST_PACK_INT_MODELID] = p->modelID;
  serializePack(tempbuf, vlistNints, CDI_DATATYPE_INT, buf, size, position, context);
  serializePack(&p->ntsteps, 1, CDI_DATATYPE_LONG, buf, size, position, context);

  cdiAttsPack(p, CDI_GLOBAL, buf, size, position, context);
  for (int varID = 0; varID < p->nvars; varID++)
    {
      vlistVarPack(p, varID, (char *) buf, size, position, context);
    }
}

int
vlistUnpack(char *buf, int size, int *position, int originNamespace, void *context, int force_id)
{
#define adaptKey(key) (namespaceAdaptKey((key), originNamespace))
  int tempbuf[vlistNints];
  serializeUnpack(buf, size, position, tempbuf, vlistNints, CDI_DATATYPE_INT, context);
  int nvars = tempbuf[VLIST_PACK_INT_NVARS];
  int targetID = force_id ? adaptKey(tempbuf[VLIST_PACK_INT_SELF]) : CDI_UNDEFID;
  vlist_t *p = vlist_new_entry(targetID);
  xassert(!force_id || p->self == targetID);
  if (!force_id) targetID = p->self;
  cdiVlistMakeInternal(p->self);
  p->taxisID = adaptKey(tempbuf[VLIST_PACK_INT_TAXISID]);
  p->tableID = tempbuf[VLIST_PACK_INT_TABLEID];
  p->instID = adaptKey(tempbuf[VLIST_PACK_INT_INSTID]);
  p->modelID = adaptKey(tempbuf[VLIST_PACK_INT_MODELID]);
  serializeUnpack(buf, size, position, &p->ntsteps, 1, CDI_DATATYPE_LONG, context);
  cdiAttsUnpack(targetID, CDI_GLOBAL, buf, size, position, context);
  for (int varID = 0; varID < nvars; varID++) vlistVarUnpack(targetID, buf, size, position, originNamespace, context);
  reshSetStatus(targetID, &vlistOps, reshGetStatus(targetID, &vlistOps) & ~RESH_SYNC_BIT);
#undef adaptKey
  return targetID;
}

void
vlist_check_contents(int vlistID)
{
  int nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; index++)
    {
      int zaxisID = vlistZaxis(vlistID, index);
      if (zaxisInqType(zaxisID) == ZAXIS_GENERIC) cdiCheckZaxis(zaxisID);
    }
}

/* Resizes and initializes opt_grib_kvpair data structure. */
void
resize_opt_grib_entries(var_t *var, int nentries)
{
  if (var->opt_grib_kvpair_size >= nentries)
    {
      return; /* nothing to do; array is still large enough */
    }
  else
    {
      if (CDI_Debug) Message("resize data structure, %d -> %d", var->opt_grib_kvpair_size, nentries);

      int new_size = (2 * var->opt_grib_kvpair_size) > nentries ? (2 * var->opt_grib_kvpair_size) : nentries;
      opt_key_val_pair_t *tmp = (opt_key_val_pair_t *) Malloc((size_t) new_size * sizeof(opt_key_val_pair_t));
      for (int i = 0; i < var->opt_grib_kvpair_size; ++i)
        {
          tmp[i] = var->opt_grib_kvpair[i];
        }
      for (int i = var->opt_grib_kvpair_size; i < new_size; ++i)
        {
          tmp[i].int_val = 0;
          tmp[i].dbl_val = 0;
          tmp[i].update = false;
          tmp[i].keyword = NULL;
        }  // for
      var->opt_grib_kvpair_size = new_size;
      Free(var->opt_grib_kvpair);
      var->opt_grib_kvpair = tmp;
    }
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
