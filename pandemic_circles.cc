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
 * A small program computing routes approximating a circle given a center point and a radius. For further information, see README.md
 */
int main(int argc, char *argv[])
{
/* filters for ways to be used */

	const map<string, set<string>> include_tags = { {"highway", {"track", "residential", "trunk", "primary", "secondary", "tertiary", "unclassified", "trunk_link", "primary_link", "secondary_link", "tertiary_link", "living_street", "road"}} };
	const map<string, set<string>> exclude_tags = { {"bicycle", {"no", "dismount"}},
													{"area", {"yes"}},
													{"access", {"no", "private"}}};

/* input data : longitude of center, latitude of center, radius of circle, mean radius of earth */

	vector<string> arguments;
	for(unsigned int arg = 0; arg < (unsigned int)argc; ++arg)
		arguments.push_back(argv[arg]);
	double lat = std::numeric_limits<double>::max();	// latitude of center of circle (degrees)
	double lon = std::numeric_limits<double>::max();	// longitutde of center of circle (degrees)
	double r = std::numeric_limits<double>::max();		// radius of circle (km)
	string osm_file = "";								// the input OSM file for the way and node data
	string svg_output_file = "";						// the file into which svg output is written (if not provided, no svg is written)
	string gpx_output_file = "";						// the file into which gpx output is written (if not provided, no gpx is written)
	for(const auto& argument : arguments)
	{
		if(argument.find("lat=") == 0)
			lat = stod(argument.substr(4));
		else if(argument.find("lon=") == 0)
			lon = stod(argument.substr(4));
		else if(argument.find("r=") == 0)
			r = stod(argument.substr(2));
		else if(argument.find("file=") == 0)
			osm_file = argument.substr(5);
		else if(argument.find("svg_output_file=") == 0)
			svg_output_file = argument.substr(16);
		else if(argument.find("gpx_output_file=") == 0)
			gpx_output_file = argument.substr(16);
	}
	if(lat == std::numeric_limits<double>::max())
	{
		cout << "Could not determine latitude of center point. Did you provide lat=... as command line argument?" << endl;
		return 0;
	}
	if(lon == std::numeric_limits<double>::max())
	{
		cout << "Could not determine longitude of center point. Did you provide lon=... as command line argument?" << endl;
		return 0;
	}
	if(r == std::numeric_limits<double>::max())
	{
		cout << "Could not determine radius of circle. Did you provide r=... as command line argument?" << endl;
		return 0;
	}
	if(osm_file == "")
	{
		cout << "Could not determine OSM file name. Did you provide file=... as command line argument?" << endl;
		return 0;
	}
	const double r_earth = 6371.0;				// mean radius of earth	(km)
	if(r >= M_PI * r_earth)
	{
		cout << "The radius may not be larger than " << M_PI * r_earth << " km." << endl;
		return 0;
	}

/* filters for extraction of way and node data from OSM */

	// define the function filtering the ways to be included based on the tags of the way
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
	cout << "Read " << graph.get_n_nodes() << " nodes and " << graph.get_n_edges() << " edges from file " << osm_file << endl;
	if(graph.get_n_nodes() == 0)
	{
		cout << "The graph does not have any edges and nodes. Is the region corresponding to the center point and the radius included in your *.osm.pbf file?" << endl;
		return 0;
	}

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
	{
		cout << "Did not find a circular path" << endl;
		return 0;
	}
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
	const auto path_adder = [&](list<pair<const Edge*, Direction>>::iterator begin_it, list<pair<const Edge*, Direction>>::iterator end_it) -> void
	{
		for(auto edge_it = begin_it; edge_it != end_it; ++edge_it)
			route.push_back(*edge_it);
	};
	path_adder(way_out.begin(), way_out.end());
	path_adder(start_edge, circular_path.end());
	path_adder(circular_path.begin(), start_edge);
	path_adder(way_in.begin(), way_in.end());
	cout << "Successfully generated route" << endl;

/* Write to files */

	// mark computed path
	for(const auto& edge : route)
		edge.first->set_user_flag();
	// write svg if requested
	if(svg_output_file != "")
	{
		graph.write_svg(svg_output_file, azimuthal_projection, true);
		cout << "svg file written" << endl;
	}
	// write gpx if requested
	if(gpx_output_file != "")
	{
		graph.write_gpx(gpx_output_file, route, gpx_output_file);
		cout << "gpx file written" << endl;
	}

	return 0;
}
