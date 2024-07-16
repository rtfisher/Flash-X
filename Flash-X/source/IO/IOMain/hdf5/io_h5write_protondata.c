#include <stdio.h>
#include <hdf5.h>
#include <stdlib.h>
#include "mangle_names.h"
#include "constants.h"
#include "hdf5_flash.h"
#include "Simulation.h"

#ifdef FLASH_IO_ASYNC_HDF5
extern hid_t io_es_id;
#endif

void FTOC(io_h5write_protondata)(hid_t *pFileID,
                                 int *rank,
			         int *localPoints,
			         int *globalPoints,
			         int *startPoint,
			         int *error,
			         double t[],
			         double x[],
			         double y[], 
			         double z[])
{
  /* Open the 'ProtonData' dataset inside the HDF5 plot file.
     It should already exist */
#ifdef FLASH_IO_ASYNC_HDF5
  hid_t dataset = H5Dopen_async(*pFileID, "ProtonData",H5P_DEFAULT, io_es_id);
#else
  hid_t dataset = H5Dopen(*pFileID, "ProtonData",H5P_DEFAULT);
#endif

  if (dataset < 0)
      {
       printf (" Error in opening HDF5 dataset: ProtonData");
       *error = dataset;
       return;
      }

  /* Find out how big the dataset already is. The result will be placed
     in 2D array dims. Array dims will contain the current dimensions.
     The maximum dimensions are not needed -> NULL is passed to indicate
     no interest in these. */

  hid_t dataspace = H5Dget_space (dataset);
  if (dataspace < 0)
      {
       printf (" Error in getting HDF5 dataspace for: ProtonData");
       *error = dataspace;
       return;
      }

  hsize_t dims[2];
  *error = H5Sget_simple_extent_dims (dataspace, dims, NULL);
  if (*error < 0)
      {
       printf (" Error in getting HDF5 dimensions for dataspace: ProtonData");
       return;
      }

  /* Return the old dataspace handle. We must get a new one after we increase
     the dimensions of the dataset 'ProtonData'. */

  H5Sclose (dataspace);

  /* Extend the dimensions of the 'ProtonData' dataset by adding
     all the new 'globalPoints' entries. Get a new handle on the
     extended dataspace.*/

  hsize_t newdims[2] = {dims[0] + *globalPoints, 4};
  *error = H5Dextend (dataset, newdims);
  if (*error < 0)
      {
       printf (" Error in extending HDF5 dimensions for dataspace: ProtonData");
       return;
      }

  dataspace = H5Dget_space (dataset);
  if (dataspace < 0)
      {
       printf (" Error in getting HDF5 extended dataspace for: ProtonData");
       *error = dataspace;
       return;
      }

  /* Set offsets, hyperslab size and maximum dimensions. */

  hsize_t offset  [2] = {dims[0] + *startPoint, 0};
  hsize_t slabsize[2] = {*localPoints, 4};
  hsize_t maxsize [2] = {H5S_UNLIMITED, 4};

  /* printf("%i START_POS        = %i\n", *rank, *startPoint); */
  /* printf("%i OLD DATASET SIZE = %i\t%i\n", *rank, dims[0], dims[1]); */
  /* printf("%i NEW DATASET SIZE = %i\t%i\n", *rank, newdims[0], newdims[1]); */
  /* printf("%i HYPERSLAB OFFSET = %i\t%i\n", *rank, offset[0], offset[1]); */
  /* printf("%i HYPERSLAB SIZE   = %i\t%i\n\n", *rank, slabsize[0], slabsize[1]); */

  /* Select a hyperslab region in the current extended HDF5 dataspace and
     create a handle to the memory space. */

  *error = H5Sselect_hyperslab (dataspace,
                                H5S_SELECT_SET,
                                offset,          /* starting point of the hyperslab         */
                                NULL,            /* implicit stride -> 1 value              */
			        slabsize,        /* size of the 2D hyperslab                */
                                NULL);           /* indicates 1x1 blocks -> single elements */
  if (*error < 0)
      {
       printf (" Error in selecting HDF5 hyperslab for: ProtonData");
       return;
      }

  hid_t memspace = H5Screate_simple (2, slabsize, maxsize);
  if (memspace < 0)
      {
       printf (" Error in creating memory space in HDF5 hyperslab for: ProtonData");
       *error = memspace;
       return;
      }

  /* Create the local slab and write the proton data to it. */

  int i;
  double *slab = (double*) malloc( (*localPoints)*4 * sizeof(double));

  for (i = 0; i < *localPoints; i++)
       {
        slab [i*4 + 0] = t [i];
        slab [i*4 + 1] = x [i];
        slab [i*4 + 2] = y [i];
        slab [i*4 + 3] = z [i];
       }

  /* Prepare for data transfer to HDF5 file and write the hyperslab to it. */

  hid_t datatransfer = H5Pcreate (H5P_DATASET_XFER);

#ifdef H5_HAVE_PARALLEL
    if (HDF5_MODE == COLLECTIVE)
        {
         *error = H5Pset_dxpl_mpio (datatransfer, H5FD_MPIO_COLLECTIVE);
         if (*error < 0)
	     {
              printf (" Error in setting MPI collective data transfer mode for: HDF5 ProtonData");
              return;
             }
        }
    else
        {
         *error = H5Pset_dxpl_mpio (datatransfer, H5FD_MPIO_INDEPENDENT);
         if (*error < 0)
	     {
              printf (" Error in setting MPI independent data transfer mode for: HDF5 ProtonData");
              return;
             }
        }
#endif

#ifdef FLASH_IO_ASYNC_HDF5
  *error = H5Dwrite_async (dataset,
                     H5T_NATIVE_DOUBLE,
                     memspace,
                     dataspace,
                     datatransfer,
                     slab, io_es_id);
#else
  *error = H5Dwrite (dataset,
                     H5T_NATIVE_DOUBLE,
                     memspace,
                     dataspace,
  		     datatransfer,
                     slab);
#endif


  if (*error < 0)
      {
       printf (" Error in writing local slab to HDF5 dataset -> ProtonData");
       return;
      }

  /* Smooth and clean ending. */

  free (slab);

  H5Sclose (dataspace);
  H5Sclose (memspace);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif 
  H5Pclose (datatransfer);
}
