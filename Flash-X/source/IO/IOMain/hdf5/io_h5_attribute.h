#ifndef IO_H5_ATTRIBUTE_H
#define IO_H5_ATTRIBUTE_H

#include "constants.h"
#include "Simulation.h"
#include "hdf5_flash.h"
#include "io_flash.h"
#include "io_h5_type.h"
#include <assert.h>
#include <string.h>

void io_h5_attribute_create(const int myPE,
			    const hid_t fileID,
			    const int diskType,
			    const int dims,
			    const int diskSize[],
			    const char datasetName[],
			    const char attDatasetName[]);

void io_h5_attribute_write(const int myPE,
			   const hid_t fileID,
			   const int memType,
			   const char datasetName[],
			   const char attDatasetName[],
			   const void * const pData);

#endif
