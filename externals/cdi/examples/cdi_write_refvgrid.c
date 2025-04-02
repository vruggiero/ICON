#include <stdio.h>
#include "cdi.h"

#define nlon 12  // Number of longitudes
#define nlat 6   // Number of latitudes
#define nlev 5   // Number of levels
#define nts 3    // Number of time steps

int
main(void)
{
  int gridID, zaxisID1, zaxisID2, taxisID;
  int vlistID, varID1, varID2, streamID, tsID;
  size_t numMissVals = 0;
  double levs[nlev] = { 1, 2, 3, 4, 5 };
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Create a grid reference
  gridID = gridCreate(GRID_UNSTRUCTURED, nlon * nlat);
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, 123);
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, 3);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, "http://www.x.y/gridfile.nc");
  cdiDefKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, "1234569887654321", CDI_UUID_SIZE);

  // Create a surface level Z-axis
  zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a general vertical height Z-axis
  // zaxisID2 = zaxisCreate(ZAXIS_HEIGHT, nlev);
  zaxisID2 = zaxisCreate(ZAXIS_REFERENCE, nlev);
  zaxisDefLevels(zaxisID2, levs);

  cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_NLEV, nlev);
  cdiDefKeyInt(zaxisID2, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, 71);
  // zaxisDefReference(zaxisID2, "http://www.x.y/vgridfile.nc");
  cdiDefKeyBytes(zaxisID2, CDI_GLOBAL, CDI_KEY_UUID, "8765432112345678", CDI_UUID_SIZE);

  // Create a variable list
  vlistID = vlistCreate();

  // Define the variables
  varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TIME_VARYING);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARYING);

  // Define the variable names
  vlistDefVarName(vlistID, varID1, "varname1");
  vlistDefVarName(vlistID, varID2, "varname2");

  // Create a Time axis
  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset in netCDF format
  streamID = streamOpenWrite("example.grb2", CDI_FILETYPE_GRB2);
  // streamID = streamOpenWrite("example.nc", CDI_FILETYPE_NC);
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return (1);
    }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  // Loop over the number of time steps
  for (int tsID = 0; tsID < nts; tsID++)
    {
      // Set the verification date to 1985-01-01 + tsID
      taxisDefVdate(taxisID, 19850101 + tsID);
      // Set the verification time to 12:00:00
      taxisDefVtime(taxisID, 120000);
      // Define the time step
      streamDefTimestep(streamID, tsID);

      // Init var1 and var2
      for (size_t i = 0; i < nlon * nlat; i++) var1[i] = 1.1;
      for (size_t i = 0; i < nlon * nlat * nlev; i++) var2[i] = 2.2;

      // Write var1 and var2
      streamWriteVar(streamID, varID1, var1, numMissVals);
      streamWriteVar(streamID, varID2, var2, numMissVals);
    }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID1);
  zaxisDestroy(zaxisID2);
  gridDestroy(gridID);

  return 0;
}
