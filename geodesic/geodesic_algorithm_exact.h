//Copyright (C) 2008 Danil Kirsanov, MIT License
#pragma once
#ifndef GEODESIC_ALGORITHM_EXACT_20071231
#define GEODESIC_ALGORITHM_EXACT_20071231

#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
#include <vector>
#include <cmath>
#include <assert.h>
#include <set>

namespace geodesic{

class GeodesicAlgorithmExact : public GeodesicAlgorithmBase
{
public:
	GeodesicAlgorithmExact(geodesic::Mesh* mesh):
	  	GeodesicAlgorithmBase(mesh),
		m_memory_allocator((unsigned int)mesh->edges().size(), (unsigned int)mesh->edges().size()),
		m_edge_interval_lists(mesh->edges().size())
	{
		m_type = EXACT;

		for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
		{
			m_edge_interval_lists[i].initialize(&mesh->edges()[i]);
		}
		m_collisionTwoPathEdges.clear();
	};	

	~GeodesicAlgorithmExact(){};

	void propagate(std::vector<SurfacePoint>& sources,
   				   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
				   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

	void trace_back(SurfacePoint& destination,		//trace back piecewise-linear path
					std::vector<SurfacePoint>& path);
	void trace_back_edge(edge_pointer destination,		//trace back piecewise-linear edge path
						std::vector<edge_pointer>& edge_path, int level = 0);
	void trace_back_interval(	interval_pointer destination,		//trace back piecewise-linear edge path
								std::vector<edge_pointer>& edge_path, int level = 0);

	unsigned best_source(SurfacePoint& point,			//quickly find what source this point belongs to and what is the distance to this source
		double& best_source_distance); 

	void print_statistics();

	const std::vector<edge_pointer>&  CollisionEdges(){return m_collisionTwoPathEdges;}
private:
	typedef std::set<interval_pointer, Interval> IntervalQueue;
	Interval::DirectionType get_edge_source_direction(edge_pointer e);
	void update_list_and_queue(list_pointer list,
							   IntervalWithStop* candidates,	//up to two candidates
							   unsigned num_candidates);

	unsigned compute_propagated_parameters(double pseudo_x, 
											double pseudo_y, 
											double d,		//parameters of the interval
											double start, 
											double end,		//start/end of the interval
											double alpha,	//corner angle
											double L,		//length of the new edge
											bool first_interval,		//if it is the first interval on the edge
											bool last_interval,
											bool turn_left,
											bool turn_right,
											IntervalWithStop* candidates);		//if it is the last interval on the edge

	void construct_propagated_intervals(bool invert, 
									  edge_pointer edge, 
									  face_pointer face,		//constructs iNew from the rest of the data
									  IntervalWithStop* candidates,
									  unsigned& num_candidates,
									  interval_pointer source_interval);

	double compute_positive_intersection(double start,
										 double pseudo_x,
										 double pseudo_y,
										 double sin_alpha,
										 double cos_alpha);		//used in construct_propagated_intervals

	unsigned intersect_intervals(interval_pointer zero, 
								    IntervalWithStop* one);			//intersecting two intervals with up to three intervals in the end

	interval_pointer best_first_interval(SurfacePoint& point, 
										double& best_total_distance, 
										double& best_interval_position,
										unsigned& best_source_index);
	interval_pointer best_min_interval_ofEdge(edge_pointer edge, 
											  double& best_total_distance);
	interval_pointer best_max_interval_ofEdge(edge_pointer edge, 
											  double& best_total_distance);
	bool check_stop_conditions(unsigned& index);
	void compute_collision_edges();
	void compute_collision_edges2();
	void clear()
	{
		m_memory_allocator.clear();
		m_queue.clear();
		for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
		{
			m_edge_interval_lists[i].clear();
		}
		m_propagation_distance_stopped = GEODESIC_INF;
	};

	list_pointer interval_list(edge_pointer e)
	{
		return &m_edge_interval_lists[e->id()];
	};

	void set_sources(std::vector<SurfacePoint>& sources)
	{
		m_sources.initialize(sources);
	}

	void initialize_propagation_data();		

	void list_edges_visible_from_source(MeshElementBase* p,
										std::vector<edge_pointer>& storage); //used in initialization

	long visible_from_source(SurfacePoint& point);	//used in backtracing

	void best_point_on_the_edge_set(SurfacePoint& point, 
									std::vector<edge_pointer> const& storage,
									interval_pointer& best_interval,
									double& best_total_distance,
									double& best_interval_position);

	void possible_traceback_edges(SurfacePoint& point, 
								  std::vector<edge_pointer>& storage);

	bool erase_from_queue(interval_pointer p);

	IntervalQueue m_queue;	//interval queue

	MemoryAllocator<Interval> m_memory_allocator;			//quickly allocate and deallocate intervals 
	std::vector<IntervalList> m_edge_interval_lists;		//every edge has its interval data 

	enum MapType {OLD, NEW};		//used for interval intersection
	MapType map[5];		
	double start[6];
	interval_pointer i_new[5];

	unsigned m_queue_max_size;			//used for statistics
	unsigned m_iterations;			//used for statistics

	SortedSources m_sources;

	std::vector<edge_pointer> m_collisionTwoPathEdges;
};



}		//geodesic

#endif //GEODESIC_ALGORITHM_EXACT_20071231
