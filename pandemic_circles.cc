// --------------------------------------------------------------------------
// Copyright (C) 2021 by Sebastian Stark
//
// This file is part of PandemicCircles
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <graph_lib/graph.h>
#include <iostream>
#include <math.h>
#include <bits/stdc++.h>

using namespace GraphLib;
using namespace std;

int main()
{

	// input data : longitude of center, latitude of center, radius of circle, mean radius of earth
	const double lat = 51.05871;
	const double lon = 13.77054;
	const double r = 15.0;
	const double r_earth = 6371.0;
	assert( r < M_PI * r_earth );

	// define an azimuthal equidistant projection
	const auto azimuthal_projection = [&](const coordinate_t& p) -> coordinate_t
	{
		coordinate_t p_r, res;
		const double lat_r = lat * M_PI / 180.0;
		const double lon_r = lon * M_PI / 180.0;
		p_r.first = p.first * M_PI / 180.0;
		p_r.second = p.second * M_PI / 180.0;
		const double c = acos( sin(lat_r) * sin(p_r.first) + cos(lat_r) * cos(p_r.first) * cos(p_r.second - lon_r) );
		const double k = r_earth * c / sin(c);
		res.first = k * cos(p_r.first) * sin(p_r.second - lon_r);
		res.second = k * ( cos(lat_r) * sin(p_r.first) - sin(lat_r) * cos(p_r.first) * cos(p_r.second - lon_r) );
		return res;
	};

	// standard (Euclidean) metric
	const auto standard_metric = [&](const coordinate_t& p_1, const coordinate_t& p_2) -> double
	{
		const double delta_lat = (p_1.first - p_2.first) * M_PI / 180.0;
		const double delta_lon = (p_1.second - p_2.second) * M_PI / 180.0;
		const double lat_1_r = p_1.first * M_PI / 180.0;
		const double lat_2_r = p_1.second * M_PI / 180.0;
		const double a = pow(sin(0.5 * delta_lat), 2.0) + pow(sin(0.5 * delta_lon), 2.0) * cos(lat_1_r) * cos(lat_2_r);
		const double c = 2.0 * asin(sqrt(a));
		return r_earth * c;
	};

	// distorted metric
	const auto distorted_metric = [&](const coordinate_t& p_1, const coordinate_t& p_2) -> double
	{
		auto p_1_ap = azimuthal_projection(p_1);
		auto p_2_ap = azimuthal_projection(p_2);
		const double dx = p_1_ap.first - p_2_ap.first;
		const double dy = p_1_ap.second - p_2_ap.second;
		const double r_1 = sqrt( p_1_ap.first * p_1_ap.first + p_1_ap.second * p_1_ap.second );
		const double r_2 = sqrt( p_2_ap.first * p_2_ap.first + p_2_ap.second * p_2_ap.second );
		const double r_avg = 0.5 * (r_1 + r_2);
		// this is only approximate as analytical integration is only possible for special cases and, even then, the formulae are quite complex
		return sqrt( (dx * dx + dy * dy) / r_avg / r_avg ) / pow(r_avg / r, 10.0);
	};


	// define the function filtering the ways to be included based on the tags of the way
	const map<string, set<string>> exclude_tags = { {"bicycle", {"no", "dismount"}}, {"area", {"yes"}}, {"access", {"no", "private"}}};
	const map<string, set<string>> include_tags = { {"highway", {"residential", "trunk", "primary", "secondary", "tertiary", "unclassified", "trunk_link", "primary_link", "secondary_link", "tertiary_link", "living_street", "road"}} };
	const auto include_way = [&](const osmium::Way& way) -> bool
	{
		bool include = false;
		for(const osmium::Tag& t : way.tags())
		{
			const auto t_it_ex = exclude_tags.find(string(t.key()));
			if( (t_it_ex != exclude_tags.end()) && (t_it_ex->second.find(string(t.value())) != t_it_ex->second.end()) )
				return false;
			const auto t_it_in = include_tags.find(string(t.key()));
			if( (t_it_in != include_tags.end()) && (t_it_in->second.find(string(t.value())) != t_it_in->second.end()) )
				include = true;
		}
		return include;
	};

	// define the function filtering the nodes to be included
	const auto include_node = [&](const osmium::NodeRef& node) -> bool
	{
		const double delta_lat = (node.lat() - lat) * M_PI / 180.0;
		const double delta_lon = (node.lon() - lon) * M_PI / 180.0;
		const double lat_1_r = lat * M_PI / 180.0;
		const double lat_2_r = node.lat() * M_PI / 180.0;
		const double a = pow(sin(0.5 * delta_lat), 2.0) + pow(sin(0.5 * delta_lon), 2.0) * cos(lat_1_r) * cos(lat_2_r);
		const double c = 2.0 * asin(sqrt(a));
		const double d = r_earth * c;
		if( (d < r) && (d > 0.001)) // safety distance around origin to later avoid potential issues with singular metric
			return true;
		else
			return false;
	};

	// set up graph (start with distorted metric since the exterior circle is to be computed first)
	Graph graph(distorted_metric);

	// read the osm.pbf file
	graph.read_graph_from_osm("dresden.osm.pbf", include_way, include_node, true);

	// cut graph along (0, y>0) before computation of exterior circle
	priority_queue< pair<double, pair<const Edge*, Direction>>, vector <pair<double, pair<const Edge*, Direction>>> > cuts;
	for(const auto& edge : graph)
	{
		const auto& coordinates = edge.get_coordinates();
		unsigned int cut_intersections = 0;
		for(unsigned int m = 1; m < coordinates.size(); ++m)
		{
			const auto p_1 = azimuthal_projection(coordinates[m-1]);
			const auto p_2 = azimuthal_projection(coordinates[m]);
			// note: shift any path point sitting closer than 1e-12 to the cutting edge slightly towards positive x-direction in order to get well-defined cutting procedure
			double p_1_x = fabs(p_1.first) < 1e-12 ? 1e-12 : p_1.first;
			double p_1_y = p_1.second;
			double p_2_x = fabs(p_2.first) < 1e-12 ? 1e-12 : p_2.first;
			double p_2_y = p_2.second;

			// there can only be an intersection if the sign of the x-coordinate of the two points is different
			if(p_1_x * p_2_x < 0.0)
			{
				// cannot run into division by zero here since p_1_x * p_2_x < 0 and |p_1_x| > 1e-12 and |p_2_x| > 1e-12
				const double t = -p_1_x / (p_2_x - p_1_x);
				assert( ( "This should not have happened and seems to be a bug", (t > 0.0) && (t < 1.0) ) );
				if(p_1_y + t * (p_2_y - p_1_y) > 0.0)
					++cut_intersections;
			}
		}
		// if the edge makes an odd number of intersections with the cutting line, block it
		if( (cut_intersections % 2 != 0) && (edge.get_direction() != Direction::none))
		{
			const auto p_1 = azimuthal_projection(coordinates.front());
			const auto p_2 = azimuthal_projection(coordinates.back());
			const double y_cut = 0.5 * (p_1.second + p_2.second);
			cuts.push(make_pair(y_cut, make_pair(&edge, edge.get_direction())));
			edge.set_direction(Direction::none);
		}
	}
	// keep a copy of cuts in order to undo them later
	auto cuts_copy = cuts;

	// calculate circular path (start with the cut with the largest distance from the center; if this does not work, try the next one, etc.)
	list<pair<const Edge*, Direction>> circular_path;
	while(!cuts.empty())
	{
		const auto& cut = cuts.top();
		const auto& edge = *(cut.second.first);
		const auto& direction = cut.second.second;
		Node node_start, node_end;
		if( (direction == Direction::forward) || (direction == Direction::both))
		{
			node_start = edge.get_node_2();
			node_end = edge.get_node_1();
		}
		else if(direction == Direction::backward)
		{
			node_start = edge.get_node_1();
			node_end = edge.get_node_2();
		}
		else
			continue;
		circular_path = graph.compute_shortest_path(node_start, node_end);
		if(circular_path.size() > 0)	// successful
		{
			// close the circle in break the loop
			if(edge.get_node_1() == node_end)
				circular_path.push_back(make_pair(&edge, Direction::forward));
			else
				circular_path.push_back(make_pair(&edge, Direction::backward));
			break;
		}
		cuts.pop();
	}
	if(circular_path.size() == 0)
		assert(("Did not find a circular path", false));

	// undo blocking of edges
	while(!cuts_copy.empty())
	{
		const auto& cut = cuts_copy.top();
		const auto& edge = *(cut.second.first);
		const auto& direction = cut.second.second;
		edge.set_direction(direction);
		cuts_copy.pop();
	}

	// use standard metric for calculation of route from center to exterior circle and back
	graph.set_metric(standard_metric);

	// calculate way from center to exterior circle
	double min_dist = numeric_limits<double>::max();
	Node min_dist_node;
	list<pair<const Edge*, Direction>>::iterator start_edge;
	for(auto edge_it = circular_path.begin(); edge_it != circular_path.end(); ++edge_it)
	{
		if(edge_it->second == Direction::forward)
		{
			const auto& p = edge_it->first->get_coordinates().front();
			if(standard_metric(p, make_pair(lat,lon)) < min_dist)
			{
				min_dist = standard_metric(p, make_pair(lat,lon));
				min_dist_node = edge_it->first->get_node_1();
				start_edge = edge_it;
			}
		}
		else
		{
			const auto& p = edge_it->first->get_coordinates().back();
			if(standard_metric(p, make_pair(lat,lon)) < min_dist)
			{
				min_dist = standard_metric(p, make_pair(lat,lon));
				min_dist_node = edge_it->first->get_node_2();
				start_edge = edge_it;
			}
		}
	}
	list<pair<const Edge*, Direction>> way_out, way_in;
	set<Node> tried_start_nodes;
	// Todo: strictly, it would be necessary to check here that the start node is not already on the exterior circle
	while(tried_start_nodes.size() < graph.get_n_nodes())
	{
		const Node start_node = graph.get_closest_node(make_pair(lat, lon), tried_start_nodes);
		tried_start_nodes.insert(start_node);
		way_out = graph.compute_shortest_path(start_node, min_dist_node);
		way_in = graph.compute_shortest_path(min_dist_node, start_node);
		if( (way_out.size() != 0) && (way_in.size() != 0) )
			break;
	}

	// finally assemble everything into a single path and calculate distance
	list<pair<const Edge*, Direction>> path;
	double distance = 0.0;
	const auto path_adder = [&](list<pair<const Edge*, Direction>>::iterator begin_it, list<pair<const Edge*, Direction>>::iterator end_it) -> void
	{
		for(auto edge_it = begin_it; edge_it != end_it; ++edge_it)
		{
			path.push_back(*edge_it);
			distance += edge_it->first->get_length(standard_metric);
		}
	};
	path_adder(way_out.begin(), way_out.end());
	path_adder(start_edge, circular_path.end());
	path_adder(circular_path.begin(), start_edge);
	path_adder(way_in.begin(), way_in.end());
	cout << "Calculated loop with total distance of " << distance << " km" << endl;

	// write to svg file
	for(const auto& edge : path)
		edge.first->set_user_flag();
	graph.write_svg("map.svg", azimuthal_projection, true);




/*	graph.add_edge(1, 0, {{1.0, 2.0},{0.0, 1.0}}, Direction::both);
	graph.add_edge(0, 7, {{0.0, 1.0},{1.0, 0.0}}, Direction::both);
	graph.add_edge(1, 7, {{1.0, 2.0},{1.0, 0.0}}, Direction::both);
	graph.add_edge(1, 2, {{1.0, 2.0},{2.0, 2.0}}, Direction::both);
	graph.add_edge(7, 6, {{1.0, 0.0},{2.0, 0.0}}, Direction::both);
	graph.add_edge(8, 7, {{2.0, 1.0},{1.0, 0.0}}, Direction::forward);
	graph.add_edge(6, 8, {{2.0, 0.0},{2.0, 1.0}}, Direction::both);
	graph.add_edge(2, 8, {{2.0, 2.0},{2.0, 1.0}}, Direction::both);
	graph.add_edge(6, 5, {{2.0, 0.0},{3.0, 0.0}}, Direction::backward);
	graph.add_edge(2, 5, {{2.0, 2.0},{3.0, 0.0}}, Direction::forward);
	graph.add_edge(2, 3, {{2.0, 2.0},{3.0, 2.0}}, Direction::backward);
	graph.add_edge(5, 3, {{3.0, 0.0},{3.0, 2.0}}, Direction::both);
	graph.add_edge(3, 4, {{3.0, 2.0},{4.0, 1.0}}, Direction::both);
	graph.add_edge(5, 4, {{3.0, 0.0},{4.0, 1.0}}, Direction::backward);*/



/*	auto shortest_path = graph.compute_shortest_path(src, dest);
	for(const auto& e : shortest_path)
		e.first->set_user_flag();*/



	return 0;
}
