# aggregating-urban-heat-demand
The script included here aims at aggregating buildings in urban space, so that heat demand characteristics can be made publicly available without breaking data protection requirements. 



## Script inputs
- **input_layer_name**
- **gml_id_field**: used to uniquely identify each building
- **gfk_class_field**: used in the estimation of the units per building. See below
- **build_num_floors_field**: number of floors, also used for the estimation of units per building
- **urb_block_field**: the unique id of urban block in which each building resides, used in forming the id of the newly formed groups
- **plot_num_field**: used for reweighting. Giving a unique string would negate any reweigting
- **res_fl_area_field**: total residential floor area of the building
- **nonres_fl_area_field**: total non-residential floor area of the building
- **heww_unsan_field**: energy for heating and domestic hot water in a non-renovated building state
- **heww_san1_field**: energy for heating and domestic hot water in a renovated building state
- **min_group_size**: minimum number of units per building group
- **max_dist_factor**: the distance factor, refer to the paper.
- **cluster_output_file**: path to save the output. Per default in a sqlite database. 
- **building_output_file**: makes a copy of the input buildings with an added group field "cluster". Per default in a sqlite database.

## Script Outputs
