#include "mangle_names.h"
#include <string.h>
#include <assert.h>
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Simulation.h"
#include "constants.h"
#include "io_flash.h"
#include "io_h5_attribute.h"

#ifdef FLASH_HDF5_ZFP
#include "H5Zzfp_lib.h"
#include "H5Zzfp_props.h"
#endif

#ifdef FLASH_HDF5_SZ
#include "H5Z_SZ.h"
#endif

#ifdef FLASH_HDF5_SZ3
#include "hdf5_sz3/include/H5Z_SZ3.hpp"
#endif

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

int Driver_abortC(char* message);

/* 
   This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from an unk variable
   
   
   The dimensions of the unknowns array (nxb, nyb, nzb, maxblocks)
   are passed through as arguments.  The dataspace (what gets written
   to disk) and the memory space (how the unknowns array is accessed in
   local memory) are defined based on these passed values. 
   We also pass in the maximum and minimum values of the dataset 

*/


void FTOC(io_h5write_unknowns)(int* myPE,
			       hid_t* file_identifier,
			       int* nxb,             /* num of zones in x dir*/
			       int* nyb,             /* num of zones in y */
			       int* nzb,             /* num of zones in z */
			       double* unknowns,     /* [NZB][NYB][NXB][var] */
			       double* varMin, 
			       double* varMax,
			       char record_label[5], /* add char-null termination */
			       int* local_blocks,  
			       int* total_blocks,
			       int* global_offset,
			       int* compression_type,
			       double* compression_val,
			       double* compression_val2,
			       int* dowrite)
{
  int rank;
  hid_t dataspace, dataset, memspace, dxfer_template, dataset_plist;
  herr_t status, ierr, herr;
  hsize_t dimens_4d[4], dimens_5d[5], chunk_dims[4], start_4d[4], stride_4d[4], count_4d[4];
  size_t cd_nelmts;
  unsigned int* cd_values = NULL;
  char record_label_new[5];
  const char minLabel[] = "minimum";
  const char maxLabel[] = "maximum";
  const int attType = IO_FLASH_DOUBLE;
  const int dims = 1;
  const int diskSize[] = {1};
  int parallelIO;
#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif

#define TRY_AND_CHECK_IGNORE(call, message) if ((call) < 0) \
        {fprintf(stderr,"Error (ignored) from %s in io_h5write_unknowns.\n",message);}

  /* 
     the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination 
  */
  strncpy(record_label_new, record_label,4);
  *(record_label_new + 4) = '\0';

  /* set the dimensions of the dataset */
  rank = 4;
  dimens_4d[0] = *total_blocks;
  dimens_4d[1] = *nzb;
  dimens_4d[2] = *nyb;
  dimens_4d[3] = *nxb;

  /* if (dimens_4d[0] == 1) { */
  /*   dimens_4d[0] = 2; */
  /* } */

  /* printf("dspace: %lu %lu %lu %lu\n",  dimens_4d[0], dimens_4d[1], dimens_4d[2], dimens_4d[3]); */

  dataspace = H5Screate_simple(rank, dimens_4d, NULL);
  if(dataspace < 0) Driver_abortC("Error: negative return from dataspace H5Screate_simple");

  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);
  if(dataset_plist < 0) Driver_abortC("Error: negative return from dataset_plist H5Pcreate");
  H5Pset_fill_time(dataset_plist, H5D_FILL_TIME_NEVER);

  if ((*myPE == MASTER_PE) && (*global_offset != 0)) {
    
    /* open a parallel hdf5 dataset */
#ifdef FLASH_IO_ASYNC_HDF5
    dataset = H5Dopen_async(*file_identifier, record_label_new,H5P_DEFAULT, io_es_id); 
#else
    dataset = H5Dopen(*file_identifier, record_label_new,H5P_DEFAULT); 
#endif    
    if(dataset < 0) {
      Driver_abortC("Error: H5Dopen io_h5write_unk\n");
    }
  } else {
    
#ifdef FLASH_HDF5_ZFP
    if (H5Z_zfp_initialize() < 0)
      Driver_abortC("Error: H5Z_zfp_initialize\n");
#endif

#if (defined FLASH_HDF5_ZFP) || (defined FLASH_HDF5_SZ) || (defined FLASH_HDF5_SZ3) || (defined CHUNK)
    chunk_dims[0] = 1;
    chunk_dims[1] = *nzb;
    chunk_dims[2] = *nyb;
    chunk_dims[3] = *nxb;
    H5Pset_chunk(dataset_plist, rank, chunk_dims);
    H5Pset_alloc_time(dataset_plist, H5D_ALLOC_TIME_LATE);
    /* printf("chunk size: %lu %lu %lu %lu\n", chunk_dims[0], chunk_dims[1], chunk_dims[2], chunk_dims[3]); */
#endif    

    switch(*compression_type) {
#ifdef FLASH_HDF5_SZ
      case 1:
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, ABS, *compression_val, 0, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;
            
      case 2:
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, REL, 0, *compression_val, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 3:
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, PW_REL, 0, 0, *compression_val, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 4:
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, PSNR, 0, 0, 0, *compression_val);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;
#endif    

#ifdef FLASH_HDF5_SZ3
      // EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL
      case 1:
        // ABS
        if(*myPE == MASTER_PE) 
            printf("Using SZ ABS %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 0, *compression_val, 0, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;
            
      case 2:
        // REL
        if(*myPE == MASTER_PE) 
            printf("Using SZ REL %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 1, 0, *compression_val, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 3:
        // Not supported with SZ3
        if(*myPE == MASTER_PE) 
            printf("Not supported by SZ3 %f %f\n",  *compression_val, *compression_val2);
        break;

      case 4:
        // NORM
        if(*myPE == MASTER_PE) 
            printf("Using SZ NORM %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 3, 0, 0, 0, *compression_val);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 5:
        // PSNR
        if(*myPE == MASTER_PE) 
            printf("Using SZ PSNR %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 2, 0, 0, *compression_val, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 6:
        // ABS AND REL
        if(*myPE == MASTER_PE) 
            printf("Using SZ ABS_AND_REL %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 4, *compression_val, *compression_val2, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;

      case 7:
        // ABS OR REL
        if(*myPE == MASTER_PE) 
            printf("Using SZ ABS_OR_REL %f %f\n",  *compression_val, *compression_val2);
        SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 5, *compression_val, *compression_val2, 0, 0);
        H5Pset_filter(dataset_plist, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        break;
#endif    
            
#ifdef FLASH_HDF5_ZFP
      case 101:
        H5Pset_zfp_accuracy(dataset_plist, *compression_val);
        break;

      case 102:
        H5Pset_zfp_precision(dataset_plist, *compression_val);
        break;

      case 103:
        H5Pset_zfp_rate(dataset_plist, *compression_val);
        break;

      case 104:
        H5Pset_zfp_reversible(dataset_plist);
        break;
#endif

      default:
        break;
    }

#ifdef FLASH_IO_ASYNC_HDF5
    dataset = H5Dcreate_async(*file_identifier, record_label_new,
                   H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dataset_plist, H5P_DEFAULT, io_es_id); 
#else
    dataset = H5Dcreate(*file_identifier, record_label_new,
                  H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dataset_plist, H5P_DEFAULT); 
#endif    

    if(dataset < 0) {
      Driver_abortC("Error: H5Dcreate io_h5write_unk\n");
    }    
  }

  if (*local_blocks > 0) {
    start_4d[0] = (hsize_t) (*global_offset);
    start_4d[1] = 0;
    start_4d[2] = 0;
    start_4d[3] = 0;

    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;

    count_4d[0] = (hsize_t) (*local_blocks);
    count_4d[1] = *nzb;
    count_4d[2] = *nyb;
    count_4d[3] = *nxb;

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
			       stride_4d, count_4d, NULL);

    if(ierr < 0){
      printf("%s\n", "Error: unable to select hyperslab for unknowns dataspace");
      Driver_abortC("Error: unable to select hyperslab for unknowns dataspace");
    }

    /* rank = 5; */
    /* dimens_5d[0] = *local_blocks; */
    /* dimens_5d[1] = *nzb; */
    /* dimens_5d[2] = *nyb; */
    /* dimens_5d[3] = *nxb; */
    /* dimens_5d[4] = 1; */

    rank = 4;
    dimens_5d[0] = *local_blocks;
    dimens_5d[1] = *nzb;
    dimens_5d[2] = *nyb;
    dimens_5d[3] = *nxb;
    dimens_5d[4] = 1;

    } else {
    /* MPI rank owns zero blocks - create a memory space of
       size 1 so that H5Screate_simple does not abort */
    /* rank = 5; */
    /* dimens_5d[0] = 1; */
    /* dimens_5d[1] = 1; */
    /* dimens_5d[2] = 1; */
    /* dimens_5d[3] = 1; */
    /* dimens_5d[4] = 1; */

    rank = 4;
    dimens_5d[0] = 1;
    dimens_5d[1] = 1;
    dimens_5d[2] = 1;
    dimens_5d[3] = 1;
  }

  memspace = H5Screate_simple(rank, dimens_5d, NULL);
  if(memspace < 0) Driver_abortC("Error: negative return from memspace H5Screate_simple");

  /* obtain a copy of the file transfer property list */ 
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);
  if(dxfer_template < 0) Driver_abortC("Error: negative return from dxfer_template H5Pcreate");
#ifdef IO_HDF5_PARALLEL
  if (HDF5_MODE == COLLECTIVE) {
    herr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
    assert(herr >= 0);
  }
#endif

  /* default is 'independent' mode */
  /*ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);*/
  /*ierr = H5Pset_preserve(dxfer_template, 0u);*/


  /*!!DEV -kda I've removed the max and min parts for now */

  /* write the data */
  if(!*dowrite || *local_blocks <= 0) {
    ierr = H5Sselect_none(dataspace);
    if(ierr < 0) Driver_abortC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
    ierr = H5Sselect_none(memspace);
    if(ierr < 0) Driver_abortC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
  }
#ifdef FLASH_IO_ASYNC_HDF5
  status = H5Dwrite_async(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                dxfer_template, unknowns,io_es_id);
#else
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                dxfer_template, unknowns);
#endif  

#ifdef DEBUG_IO
  printf("UNKNOWNS: wrote unknowns, status = %d\n", (int) status);
#endif
  if (status < 0){
      Driver_abortC("Error: Unable to write unknowns\n");
  }

  
  H5Pclose(dxfer_template);



/*  TRY_AND_CHECK_IGNORE( H5Fflush(*file_identifier,H5F_SCOPE_GLOBAL) , "H5Fflush" ) */

  H5Pclose(dataset_plist);

  H5Sclose(memspace); 
  H5Sclose(dataspace);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif  

  /* write the attributes -- minimum and maximum values of the data */
  if (parallelIO ||
      (!parallelIO && (*myPE == MASTER_PE) && (*global_offset == 0))) {
    /* The logic ensures that creating and writing attributes work
       in both parallel and serial FLASH I/O modes.  In parallel all MPI
       processes create and write the attribute.  In serial the master
       MPI process creates and writes the attribute on the first
       entry to io_h5write_unknowns only.  Note: The io_h5_attribute_XXX
       functions open the dataset "record_label_new" before creating /
       writing attributes, so "record_label_new" should be closed before
       the call. */

    io_h5_attribute_create(*myPE, *file_identifier, attType, dims,
			   diskSize, record_label_new, minLabel);
    io_h5_attribute_write(*myPE, *file_identifier, attType,
			  record_label_new, minLabel, varMin);

    io_h5_attribute_create(*myPE, *file_identifier, attType, dims,
			   diskSize, record_label_new, maxLabel);
    io_h5_attribute_write(*myPE, *file_identifier, attType,
			  record_label_new, maxLabel, varMax);
  }
}
