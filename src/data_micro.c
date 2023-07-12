/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/




#include "include/data_micro.h"

int Read_Micro_LUT(char* lut_name, DATA_Micro *data_micro, DATA_config *config){

  FILE *file = NULL;
  VAR_DATA tmp1, tmp2, tmp3;
  int i = 0, nbr_lines = 0;

  char FileName[200]; //Alphas

  strcpy(FileName, config->Micro_LUT_Folder);

  strcat(FileName, lut_name);


 file = fopen(FileName, "r");
 
 if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", FileName);
    return 1;
  }



#if VAR_DATA_PRECISION==1
  while(fscanf(file, "%f %f %f", &tmp1, &tmp2, &tmp3) > 0){
#else
  while(fscanf(file, "%lf %lf %lf", &tmp1, &tmp2, &tmp3) > 0){
#endif

    nbr_lines++;					// compte le nombre de ligne pour pouvoir allouer la mémoire dynamiquement
  }
  //printf("\n\n Nbr lignes: %d\n\n", nbr_lines);

  rewind(file);	// restart fle from the beginning


  data_micro->Nbr_Data = nbr_lines;
  data_micro->Energy_List = (VAR_DATA*) malloc(nbr_lines * sizeof(VAR_DATA));
  data_micro->Alpha = (VAR_DATA*) malloc(nbr_lines * sizeof(VAR_DATA));
  data_micro->Sqrt_Beta = (VAR_DATA*) malloc(nbr_lines * sizeof(VAR_DATA));

#if VAR_DATA_PRECISION==1
  while(fscanf(file, "%f %f %f", &(data_micro->Energy_List[i]), &(data_micro->Alpha[i]), &(data_micro->Sqrt_Beta[i])) > 0){
#else
  while(fscanf(file, "%lf %lf %lf", &(data_micro->Energy_List[i]), &(data_micro->Alpha[i]), &(data_micro->Sqrt_Beta[i])) > 0){
#endif
	i++;
  }

  fclose(file);



  return 0;
}

/*
void Free_BioModel_Parameters(BioModel_parameters *bio){

	if(bio->LUT_micro_P != NULL) free(bio->LUT_micro_P);
	if(bio->LUT_micro_He != NULL) free(bio->LUT_micro_He);

}
*/
