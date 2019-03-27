'''
#
#   Copyright 2018 Ivan Dochev	ivan.dochev@hcu-hamburg.de
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
'''


##input_layer_name=vector

##gml_id_field=string UUID
##gfk_klasse_field=string GFK_Klasse
##build_num_floors_field=string AOG
##urb_block_field=string baublockbe
##plot_num_field=string flurstueck
##wohnfl_field=string wohnfl
##ngw_ngf_field=string nwg_fl
##heww_unsan_field=string HEWW_total_UNSAN
##heww_san1_field=string HEWW_total_SAN1_RW

##min_cl_size=number 5
##max_dist_factor=number 5

##cluster_output_file=output vector
##ids_table_file=output file


import math 
import timeit
import operator
import itertools
from itertools import groupby

import numpy as np
from numpy import unravel_index

from qgis.core import *
from PyQt4.QtCore import QVariant

import scipy as sc
from scipy import ndimage
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse.csgraph import connected_components

import shapely.geometry as geometry
import shapely.wkb
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.ops import cascaded_union

import csv

#FUNCTIONS
def MinWE(floors, constr_type_field, function_class, geom):
    
    if geom.area() < 40:
        return 0
    #wohnen and wohnen_gemischt. In ArcGIS Tool uses wohnen_gemischt_0.85 etc..
    if 'wohnen' in function_class:
        if floors <= 3:
            return 1
        elif floors > 3 and floors <= 5:
            return 3
        elif floors > 5:
            return floors
    else:
        return 1
    
def shply_centroid(geom):
    #clean with buffer
    as_wkb = geom.buffer(0,0).asWkb()
    shply_geom = shapely.wkb.loads(as_wkb)
    centroid = QgsGeometry.fromWkt(str(shply_geom.representative_point()))
    return centroid

def max_dist_func(geom1, geom2, floors1, floors2):
    '''Function to make a flexible max dist based on the areas of the buildings
    So if two buildings are big allow a larger distance as max dist, if they 
    are rather small, then a shorter distance. The implemented way of achieving
    this is to find the diameters of the circles with the same areas as the 
    footprint area of the geometries and have the max_dist be a function of 
    their sum or equal to it...'''
    
    diameter1 = 2*math.sqrt(geom1.area() * floors1/3.14)
    diameter2 = 2*math.sqrt(geom2.area() * floors2/3.14)
    summed = diameter1 + diameter2
    
    return summed * max_dist_factor
    
def split_graph(build_list, atts, plot_num_field):
    
    #Input with attributes as array converted from list
    att = np.array(atts)
    
    #Find distances between polygon geoms and create complete graph
    weighted_dists = []
    cartesian_dists = []
    for build_pair in itertools.combinations(build_list, 2):
        distance = build_pair[0].geometry().distance(build_pair[1].geometry())

        #check distance and change '0' (when the polygons are bordering) to 0.5,
        #0 is taken as 'unconnected'
        if distance == 0:
            weighted_dists.append(0.5)
            cartesian_dists.append(0.5)
        elif build_pair[0][plot_num_field] == build_pair[1][plot_num_field]:
            #+0.5 so that they are not closer than connected buildings
            weighted_dists.append(0.05*distance + 0.5) 
            cartesian_dists.append(distance)
        else:
            weighted_dists.append(distance)
            cartesian_dists.append(distance)
            
    weighted_dists = np.array(weighted_dists)
    cartesian_dists = np.array(cartesian_dists)
        
    X = squareform(weighted_dists)
    cartesian_X = squareform(cartesian_dists)
    
    X2 = csr_matrix(X)
    
    #Minimum Spanning Tree
    mstgraph = minimum_spanning_tree(X2)
    
    #Sort edges, argsort is used on arrays not csr_matrix so .toarray()
    mstgraph = mstgraph.toarray()
    
    '''The graph is weighted and the MST is ready, however I need to get the 
    original lengths (in a way I need a suboptimal MST to begin with) so that I
    can then decide spatially if a building remains in a cluster or is too far,
    for that I use the Cartesian matrix cartesian_X, which was prepared in 
    parallel with the weighted one and take the MST and make it, well 
    'non-minimal', by introducing the real Cartesian distances of buildings 
    within the same plot. The distances between connected buildings remain 
    a fixed 0.5.'''
    
    np.copyto(mstgraph,cartesian_X,'same_kind',mstgraph>0)
    
    #Get an array of sorted (min to max) flattened indexes 
    sort_ids = np.argsort(mstgraph, axis=None)
    
    #Get them as list in max to min order (b[::-1]) and as tuples (zip)
    idlist = zip(*np.unravel_index(sort_ids[::-1], mstgraph.shape))
   
    #Start splitting
    allowed_too_small = 0
    outlabel = []
    for i in idlist:
        
        max_dist = max_dist_func(
                                build_list[i[0]].geometry(), 
                                build_list[i[1]].geometry(),
								build_list[i[0]][build_num_floors_field], 
                                build_list[i[1]][build_num_floors_field]
                                )
        
        if mstgraph[i] == 0:
            break

        else:
            #Keep Track of original weight
            origWeight = mstgraph[i]
            
            #Keep track of what is the situation before this edge was removed
            #find connected Components
            components, labels = connected_components(mstgraph, directed=False)
            #labels is an array that is ordered
            #unique are the unique values in labels used for aggregation see below
            unique = np.unique(labels)
            #aggregate with sc (scipy) based on the connected components (labels)
            aggreg = sc.ndimage.sum(att, labels, index=unique)

            too_small_count_before_this_edge = np.where(aggreg < min_cl_size)[0].size
            
            
            #Remove this edge
            mstgraph[i] = 0
    
			#Now with this edge removed
            #find connected Components
            components, labels = connected_components(mstgraph, directed=False)
            #labels is an array that is ordered
            #unique are the unique values in labels used for aggregation see below
            unique = np.unique(labels)
            #aggregate with sc (scipy) based on the connected components (labels)
            aggreg = sc.ndimage.sum(att, labels, index=unique)

            too_small_count = np.where(aggreg < min_cl_size)[0].size
            too_small_count_diff = too_small_count - too_small_count_before_this_edge
            
            '''So if the amount of tooSmall segments after the split with the 
            current edge is larger than what is allowed (default 0 allowed) but
            the edge is too long - then split anyways and 'remember' that you 
            allow now a number of tooSmall. If the amount of too Small segments
            after the split with the current edge is larger than what is 
            allowed but the edge does not 'force' a split (the length is 
            smaller than MaxDist, then return the edge - 
            mstgraph[i] = origWeight.
            Else proceed with split, so in all cases when 
            too_small_count <= allowed_too_small'''
            

            if too_small_count > allowed_too_small and origWeight >= max_dist:
                outlabel = labels
                allowed_too_small += too_small_count_diff
                #print 'split due to MaxDist + SPLIT'
                
            elif too_small_count > allowed_too_small and origWeight < max_dist:
                mstgraph[i] = origWeight
                #print 'Edge restored'
                #Redo components with restored edge
                components, labels = connected_components(mstgraph, directed=False)
                outlabel = labels
            
            else:
                #Note the output, for each iteration this changes if the new
                #labels from the last split do no break the condition above
                outlabel = labels
                #print 'normal case SPLIT'
            
            #print '--'
            
    #get the lengths of the longest edge of every cluster in order to 
    #give this to the buffer length when creating the geometry 
    unique = np.unique(outlabel)
    if np.size(unique): #check since indimage fails if unique is empty
        max_cluster_lengths = sc.ndimage.maximum(mstgraph, outlabel, index=unique)
    else:
        max_cluster_lengths = [0]
        
    
    return (list(outlabel), max_cluster_lengths)

def aggregation(
                layer, build_num_floors_field, 
                constr_type_field, 
                heww_unsan_field, 
                plot_num_field
                ):

    iter = layer.getFeatures()

    feat_list = [[feat, feat[urb_block_field], feat.geometry().area()]
                    for feat in iter 
                    if feat[heww_unsan_field]
                    and feat[heww_unsan_field] > 0] 
                    #Take only the ones that have a demand, NO min area set
    
    #SORT feat_list
    feat_list = sorted(feat_list, key=operator.itemgetter(1,2))

    #GROUP BASED ON BAUBLOCK
    grouped_list = []
    for bbgroup, group in groupby(feat_list, lambda x: x[1]):
        grouped_list.append(list(group))

    #ITERATE over urban blocks
    out_rows = []
    out_feats = []
    fields = []
    
    for group in grouped_list:
        
        urb_block = group[0][1]
        #print urb_block
        feats = [f[0] for f in group]
        attrs =[MinWE(f[0][build_num_floors_field], 
                      f[0][constr_type_field], 
                      f[0][gfk_klasse_field], 
                      f[0].geometry()) for f in group]
        
        #output from split graph is a tuple
        cluster_ids, max_lengths = split_graph(feats, attrs, plot_num_field) 

        #If single building in baublock - so no edges
        if not cluster_ids:
          cluster_ids = [0]

     
        #prepare list for geometry creation
        clusters = []
        for feat, cl_num in zip(feats, cluster_ids):
            clusters.append([feat, cl_num])
        
        #CREATE GEOMETRY - clusters are a list like [feat, cluster_id]
        clusters = sorted(clusters, key=operator.itemgetter(1))
        grouped = []
        for k, cluster in groupby(clusters, lambda x: x[1]):
		
            feats = [cl[0] for cl in cluster]
            grouped.append([k, feats])
    
        for cluster in grouped:
		
		    #get cluster_id and cluster lengths
            cluster_id = cluster[0]
            cluster_name = '%s_%s' % (urb_block, cluster_id)
            
            #CHECK COUNTS
            we = sum([MinWE(f[build_num_floors_field], 
                     f[constr_type_field], 
                     f[gfk_klasse_field], 
                     f.geometry()) for f in cluster[1]])
            
            
            if we < min_cl_size:
                for feat in cluster[1]:
                    out_rows.append((feat[gml_id_field], 'Anonymized', we))
                continue
            else:
                for feat in cluster[1]:
                    out_rows.append((feat[gml_id_field], cluster_name, we))

            #SUM areas and demands
            floor_area = 0
            heww_unsan = 0
            heww_san1 = 0
            res_fl_area = 0
            
            for f in cluster[1]:
                heww_unsan += f[heww_unsan_field]
                heww_san1 += f[heww_san1_field]
                floor_area += (f[wohnfl_field] + f[ngw_ngf_field])
                res_fl_area += f[wohnfl_field]

            unsan_m2 = heww_unsan/floor_area
            san1_m2 = heww_san1/floor_area

            #get cluster max lengths
            cl_max_length = int(max_lengths[cluster_id])
            geom_dist = cl_max_length/2.0 + 1
            
            #geometries
            new_geoms = [shapely.wkb.loads(feat.geometry().buffer(0,0).asWkb())
                                    for feat in cluster[1]]
            geoms_merged = cascaded_union(new_geoms).wkt
            qgs_geom = QgsGeometry.fromWkt(geoms_merged)

            #if not all buildings in the cluster are connected
            if qgs_geom.isMultipart():
                #concave hull via 2 buffers with mitre join
                hull = qgs_geom.buffer(geom_dist, 8,2,2,2.5).buffer(geom_dist*-1,8,2,2,2.5)
                hull = hull.combine(qgs_geom)
                #try with double distance 
                if hull.isMultipart():
                    hull = qgs_geom.buffer(geom_dist*2, 8,2,2,2.5).buffer(geom_dist*-2,8,2,2,2.5)
                    hull = hull.combine(qgs_geom)
                #try with bevel join and double distance
                if hull.isMultipart():
                    hull = qgs_geom.buffer(geom_dist*2,8).buffer(geom_dist*-2,8,2,3,1)
                    hull = hull.combine(qgs_geom)
                #remove inner rings with shapely, if still multipart, then convexHull
                try:
                    hull = QgsGeometry.fromWkt(Polygon(shapely.wkb.loads(hull.asWkb()).exterior).wkt)
                except: 
                    #print 'cluster as convex hull'
                    hull = qgs_geom.convexHull()

            #if no need for any hull, all buildings with connected walls
            else:
                hull = QgsGeometry.fromWkt(Polygon(shapely.wkb.loads(qgs_geom.asWkb()).exterior).wkt)
          
    
            #Fields
            fields = QgsFields()
            fields.append(QgsField('cluster', QVariant.String)),
            fields.append(QgsField('unit_count', QVariant.Int)),
            fields.append(QgsField('floor_area', QVariant.Int))
            fields.append(QgsField('res_fl_area', QVariant.Int))
            fields.append(QgsField('heww_unsan', QVariant.Int))
            fields.append(QgsField('heww_san1', QVariant.Int))
            fields.append(QgsField('unsan_m2', QVariant.Int))
            fields.append(QgsField('san1_m2', QVariant.Int))


            #Features
            feat = QgsFeature(fields)
            feat.setGeometry(hull)
            feat.setAttribute('cluster', cluster_name)
            feat.setAttribute('unit_count', we)
            feat.setAttribute('floor_area', floor_area)
            feat.setAttribute('res_fl_area', res_fl_area)
            feat.setAttribute('heww_unsan', heww_unsan)
            feat.setAttribute('heww_san1', heww_san1)
            feat.setAttribute('unsan_m2', unsan_m2)

            feat.setAttribute('san1_m2', san1_m2)

            
            out_feats.append(feat)
        
    return out_rows, out_feats, fields

#UPDATE FEATURES
tic=timeit.default_timer() # Timer begins
input_layer = processing.getObject(input_layer_name)


out_rows, out_feats, fields = aggregation(
                                          input_layer, 
                                          build_num_floors_field, 
                                          constr_type_field, 
                                          heww_unsan_field, 
                                          plot_num_field
                                          )

#ADD TABLE OUTPUT FROM out_rows
out_rows_header = [('gml_id', 'cluster_id', 'unit_count')] + out_rows
with open(ids_table_file,'wb') as f:
    w = csv.writer(f, delimiter=';')
    w.writerows(out_rows_header)


#NEW GEOMS
cluster_layer = QgsVectorLayer('Polygon', 'clusters', 'memory')
cluster_layer.dataProvider().addAttributes(fields) #[2] is QgsFields
cluster_layer.updateFields()
cluster_layer.dataProvider().addFeatures(out_feats)

error = QgsVectorFileWriter.writeAsVectorFormat(
                                                cluster_layer, 
                                                cluster_output_file, 
                                                'utf-8', 
                                                None, 
                                                'ESRI Shapefile'
                                                )


if error == QgsVectorFileWriter.NoError:
    #print 'clusters success!'
    pass

toc=timeit.default_timer()
#print 'Clusters generated and id table exported, Time elapsed: ' + str(toc - tic)