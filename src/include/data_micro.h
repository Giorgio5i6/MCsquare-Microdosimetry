/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_micro
#define H_data_micro

#include "define.h"
#include "struct.h"

typedef struct BioModel_parameters BioModel_parameters;
struct BioModel_parameters
{
	DATA_Micro LUT_micro_P; //LUT for P
    	DATA_Micro LUT_micro_He; //LUT for He
};


int Read_Micro_LUT(char *lut_name, DATA_Micro *data_micro, DATA_config *config);
//void Free_BioModel_Parameters(BioModel_parameters *bio);

#endif
