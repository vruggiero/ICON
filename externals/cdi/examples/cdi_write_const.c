#include <stdio.h>
#include "cdi.h"

#define nlon 12  // Number of longitudes
#define nlat 6   // Number of latitudes
#define nlev 5   // Number of levels

int
main(void)
{
  int gridID, zaxisID1, zaxisID2;
  int vlistID, varID1, varID2, streamID;
  size_t i, numMissVals = 0;
  double lons[nlon] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
  double lats[nlat] = { -75, -45, -15, 15, 45, 75 };
  double levs[nlev] = { 101300, 92500, 85000, 50000, 20000 };
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Create a regular lon/lat grid
  gridID = gridCreate(GRID_LONLAT, nlon * nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);

  // Create a surface level Z-axis
  zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a pressure level Z-axis
  zaxisID2 = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID2, levs);

  // Create a variable list
  vlistID = vlistCreate();

  // Define the variables
  varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TIME_CONSTANT);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_CONSTANT);

  // Define the variable names
  vlistDefVarName(vlistID, varID1, "varname1");
  vlistDefVarName(vlistID, varID2, "varname2");

  // Create a Time axis
  int taxisID = taxisCreate(TAXIS_RELATIVE);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset in netCDF format
  streamID = streamOpenWrite("example_const.nc", CDI_FILETYPE_NC);
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return 1;
    }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  // Set the verification date to 1985-01-01
  taxisDefVdate(taxisID, 19850101);
  // Set the verification time to 12:00:00
  taxisDefVtime(taxisID, 120000);
  // Define the time step
  streamDefTimestep(streamID, 0);

  // Init var1 and var2
  for (i = 0; i < nlon * nlat; i++) var1[i] = 1.1;
  for (i = 0; i < nlon * nlat * nlev; i++) var2[i] = 2.2;

  // Write var1 and var2
  streamWriteVar(streamID, varID1, var1, numMissVals);
  streamWriteVar(streamID, varID2, var2, numMissVals);

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  vlistDestroy(vlistID);
  zaxisDestroy(zaxisID1);
  zaxisDestroy(zaxisID2);
  gridDestroy(gridID);

  return 0;
}
