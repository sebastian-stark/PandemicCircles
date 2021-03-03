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

/**
 * A small program computing routes approximating a circle given a center point and a radius. The program is based on the library GraphLib, which allows
 * to convert OSM data into a graph.
 *
 * The computation first starts with computing a closed "exterior circle". For this computation, it is tried to
 * put all way points close to the maximum radius without ever exceeding the maximum distance (=radius of the circle) from the center.
 * For the computation of the circle, the map is cut along a single radial line starting at the center point by "blocking" edges crossing that radial line.
 * Then, a pair of corresponding nodes on the cut faces (i.e., nodes belonging to a blocked edge) are used as start and end point for the
 * computation of the circle. Since the cut face is completely blocked, the resulting route must generally enclose the center point. However,
 * this will not necessarily give a good approximation to a circle if the standard metric is used and the circle is computed with a "shortest path algorithm" (here: Dijkstra's algorithm).
 * Therefore, a "distorted metric" is used, which gives edges far away from the center point small weight, while edges close to the center get large weights assigned.
 * Furthermore, the start/end node pair of the circle is chosen such that it is as far away from the center point as possible.
 * With this approach, for most cases, the Dijkstra algorithm should give a reasonable approximation to a circle with maximum distance from the center point.
 *
 * Subsequent to the computation of the exterior circle, a "way out" from a node close to the center point to a point on the exterior circle is calculated
 * as well as a "way in" from the same point on the exterior circle back to the center point. The start node for the way out (corresponds to the end node for
 * the way in) is chosen to be as close as possible to the center point with regard to the standard metric. The node on the exterior circle (end node for way out and
 * start node for way in) is chosen to be as close to the center point as well in order to minimize the distance covered to reach the exterior circle and come back.
 * For the computations of the ways in and out, Dijkstra's algorithm is used together with the standard metric.
 *
 * Finally, the way out, the exterior circle and the way in are joined to a single route.
 *
 * The algorithm has no guarantee to succeed. In general, a sufficiently dense network of ways will be required. Also, the success of the algorithm
 * depends on the choice of the radial cut line (which is currently fixed to be the line given by points x=0, y>0 in an azimuthal equidistant projection,
 * although this may easily be changed). The implementation will try several start/end nodes for the computation of the exterior circle (starting from
 * the pair furthest away from the center point) until an exterior circle is successfully computed. Also, for the computation of the ways out and in different
 * nodes close to the center points will be tried (starting from the node closest to the center point) until the ways in and out are successfully computed.
 * Todo: Even with the described measures, it is not hard to construct examples, where the algorithm will not find any route despite the presence of a reasonable route.
 * This issue could be largely eliminated by implementing smarter algorithms for choosing the start node for the route, the cut line for the computation of the
 * exterior circle, and the node where the ways in and out hit the exterior circle (or one could use the following manual procedure: Firstly, choose the cut line, secondly compute
 * the exterior circle, thirdly manually choose the start and end point for the ways in and out, finally compute the ways in and out and join the route). Alternatively, other
 * algorithms than the one used here may be better for the computation of the exterior circle (one could e.g. devise a "contour following" algorithm).
 */
int main()
{

/* input data : longitude of center, latitude of center, radius of circle, mean radius of earth */

	const double lat = 51.05871;				// latitude of center of circle (degrees)
	const double lon = 13.77054;				// longitutde of center of circle (degrees)
	const double r = 15.0;						// radius of circle (km)
	const double r_earth = 6371.0;				// radius of earth	(km)
	const string osm_file = "dresden.osm.pbf";	// the input OSM file for the way and node data
	assert( r < M_PI * r_earth );

/* filters for extraction of way and node data from OSM */

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

/* projection functions and metric functions to be used */

	// Azimuthal equidistant projection
	// Projects input point p consisting of (latitude, longitude) onto a plane (x, y)
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

	// Standard metric
	// Computes shorted distance between two points p_1 and p_2 given in terms of (latitude, longitude)
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

	// Distorted metric
	// Points p_1 and p_2 given in terms of (latitude, longitude) are first projected onto the plane (x, y) using
	// an azimuthal equidistant projection. Then, the infinitesimal length element is defined to be dl = (dx^2 + dy^2)^(1/2) / (x^2 + y^2)^(m/2).
	// If m>>1, This gives ways further away from the origin (x=0, y=0) less weight than those close to the origin, so that ways away from the origin are preferred
	// when the total weight of a route is to be minimized.
	// Here, m=10 is chosen.
	// The implementation is approximate in that the integration over straight length segments is not performed exactly, but rather the approximation
	// delta_l = 2 (delta_x^2 + delta_y^2)^(1/2) / [(x_0^2 + y_0^2)^(1/2) + (x_1^2 + y_1^2)^(1/2)]^m is used where (x_0, y_0) and (x_1, y_1) are the two
	// projected points, and delta_x = x_0 - x_1, y = y_0 - y_1. This will work only well if the points are sufficiently close to each other.
	const auto distorted_metric = [&](const coordinate_t& p_1, const coordinate_t& p_2) -> double
	{
		auto p_1_ap = azimuthal_projection(p_1);
		auto p_2_ap = azimuthal_projection(p_2);
		const double dx = p_1_ap.first - p_2_ap.first;
		const double dy = p_1_ap.second - p_2_ap.second;
		const double r_1 = sqrt( p_1_ap.first * p_1_ap.first + p_1_ap.second * p_1_ap.second );
		const double r_2 = sqrt( p_2_ap.first * p_2_ap.first + p_2_ap.second * p_2_ap.second );
		const double r_avg = 0.5 * (r_1 + r_2);
		return sqrt( (dx * dx + dy * dy) / r_avg / r_avg ) / pow(r_avg / r, 10.0);
	};

/* Set up graph and read file */

	// set up graph (start with distorted metric since the exterior circle is to be computed first)
	Graph graph(distorted_metric);

	// read the osm.pbf file (eliminate "dead ends" as this reduces the possibility for error of the algorithm computing the exterior circle)
	// Todo: for the calculation of the ways in and out, it is likely better to take "dead ends" into account
	graph.read_graph_from_osm(osm_file, include_way, include_node, true);

/* Compute exterior circle */

	// cut graph along (0, y>0) before computation of exterior circle
	// Todo: allow for different cutting lines by allowing the azimuthal projection to be rotated
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
	// keep a copy of cuts (or, rather "edge blocks") in order to undo them later
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

/* Compute ways in and out */

	// use standard metric for calculation of route from center to exterior circle and back
	graph.set_metric(standard_metric);

	// get the node on the exterior circle being closest to the center point,
	// also store an iterator to the corresponding edge (which makes it easier to join the different route portions later)
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

	// compute the ways in and out trying different start nodes
	list<pair<const Edge*, Direction>> way_out, way_in;
	set<Node> tried_start_nodes;
	// Todo: Strictly, it would be necessary to check here that the start node is not already on the exterior circle
	// Todo: There is the possibility that the ways in and out share edges with the exterior circle, which means
	//       that one is using some edges three times. These shared edges could be eliminated from the ways in and out
	while(tried_start_nodes.size() < graph.get_n_nodes())
	{
		const Node start_node = graph.get_closest_node(make_pair(lat, lon), tried_start_nodes);
		tried_start_nodes.insert(start_node);
		way_out = graph.compute_shortest_path(start_node, min_dist_node);
		way_in = graph.compute_shortest_path(min_dist_node, start_node);
		if( (way_out.size() != 0) && (way_in.size() != 0) )
			break;
	}

/* Join the way out, the exterior circle and the way in into a single route; also compute distance of the route */

	list<pair<const Edge*, Direction>> route;
	double distance = 0.0;
	const auto path_adder = [&](list<pair<const Edge*, Direction>>::iterator begin_it, list<pair<const Edge*, Direction>>::iterator end_it) -> void
	{
		for(auto edge_it = begin_it; edge_it != end_it; ++edge_it)
		{
			route.push_back(*edge_it);
			distance += edge_it->first->get_length(standard_metric);
		}
	};
	path_adder(way_out.begin(), way_out.end());
	path_adder(start_edge, circular_path.end());
	path_adder(circular_path.begin(), start_edge);
	path_adder(way_in.begin(), way_in.end());
	cout << "Calculated loop with total distance of " << distance << " km" << endl;

/* Write to svg file */

	// mark computed path
	for(const auto& edge : route)
		edge.first->set_user_flag();
	// write
	graph.write_svg("map.svg", azimuthal_projection, true);

	return 0;
}
