/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_mhd
#define H_data_mhd

#include "define.h"
#include "struct.h"
#include "File_process.h"

enum MHD_ElementType{ 
	NotDefined,
	MET_FLOAT,
	MET_DOUBLE,
	MET_INT,
	MET_SHORT
};


typedef struct MHD_header MHD_header;
struct MHD_header{

  // Header informations
  int NDims;
  int ElementNumberOfChannels;
  int DimSize[3];
  VAR_DATA ElementSpacing[3];
  VAR_DATA Offset[3];
  enum MHD_ElementType ElementType;
  int ElementByteOrderMSB;
  char ElementDataFile[100];

};

void export_MHD_image(char *file_name, int GridSize[3], VAR_DATA VoxelLength[3], VAR_DATA Offset[3], VAR_SCORING *data);
int Parse_MHD_header(char *file_name, MHD_header *header);
VAR_DATA *import_MHD_image(char *file_name, int *GridSize, VAR_DATA *VoxelLength, VAR_DATA *Origin);

#endif
