/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_scenario_simulation
#define H_compute_scenario_simulation

#include "define.h"
#include "struct.h"
#include "compute_random.h"
#include "compute_simulation.h"
#include "compute_simulation_beamlet.h"
#include "data_beam_model.h"
#include "compute_beam_model.h"
#include "compute_4D.h"
#include "compute_treatment_uncertainties.h"
#include "data_micro.h" 

enum Scenario_selection_shape{ 
	HyperSphere,
	HyperCylinder,
	HyperCube
} ;

void Scenarios_selection_all(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path, enum Scenario_selection_shape Selection_shape);
void Scenarios_selection_reduced(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path, enum Scenario_selection_shape Selection_shape);
void Scenarios_selection_random(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path);
void Scenario_simulation(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields);

#endif
