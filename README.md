# PandemicCircles

A small program computing routes approximating a circle given a center point and a radius. The program is based on the library GraphLib, which allows
to convert OSM data into a graph.

The computation first starts with computing a closed "exterior circle". For this computation, it is tried to
put all way points close to the maximum radius without ever exceeding the maximum distance (=radius of the circle) from the center.
For the computation of the circle, the map is cut along a single radial line starting at the center point by "blocking" edges crossing that radial line.
Then, a pair of corresponding nodes on the cut faces (i.e., nodes belonging to a blocked edge) are used as start and end point for the
computation of the circle. Since the cut face is completely blocked, the resulting route must generally enclose the center point. However,
this will not necessarily give a good approximation to a circle if the standard metric is used and the circle is computed with a "shortest path algorithm" (here: Dijkstra's algorithm).
Therefore, a "distorted metric" is used, which gives edges far away from the center point small weight, while edges close to the center get large weights assigned.
Furthermore, the start/end node pair of the circle is chosen such that it is as far away from the center point as possible.
With this approach, for most cases, the Dijkstra algorithm should give a reasonable approximation to a circle with maximum distance from the center point.

Subsequent to the computation of the exterior circle, a "way out" from a node close to the center point to a point on the exterior circle is calculated
as well as a "way in" from the same point on the exterior circle back to the center point. The start node for the way out (corresponds to the end node for
the way in) is chosen to be as close as possible to the center point with regard to the standard metric. The node on the exterior circle (end node for way out and
start node for way in) is chosen to be as close as possible to the center point as well in order to minimize the distance covered to reach the exterior circle and come back.
For the computations of the ways in and out, Dijkstra's algorithm is used together with the standard metric.

Finally, the way out, the exterior circle and the way in are joined to a single route.

The algorithm has no guarantee to succeed. In general, a sufficiently dense network of ways will be required. Also, the success of the algorithm
depends on the choice of the radial cut line (which is currently fixed to be the line given by points x=0, y>0 in an azimuthal equidistant projection,
although this may easily be changed). The implementation will try several start/end nodes for the computation of the exterior circle (starting from
the pair furthest away from the center point) until an exterior circle is successfully computed. Also, for the computation of the ways out and in different
nodes close to the center points will be tried (starting from the node closest to the center point) until the ways in and out are successfully computed.
Todo: Even with the described measures, it is not hard to construct examples, where the algorithm will not find any route despite the presence of a reasonable route.
This issue could be largely eliminated by implementing smarter algorithms for choosing the start node for the route, the cut line for the computation of the
exterior circle, and the node where the ways in and out hit the exterior circle (or one could use the following manual procedure: Firstly, choose the cut line, secondly compute
the exterior circle, thirdly manually choose the start and end point for the ways in and out, finally compute the ways in and out and join the route). Alternatively, other
algorithms than the one used here may be better for the computation of the exterior circle (one could e.g. devise a "contour following" algorithm).

To compile the program:

1. install the library GraphLib (together with its dependencies libosmium and ZLIB), see https://github.com/sebastian-stark/GraphLib.
2. place library source files into some folder /path/to/folder/PandemicCircles (you can use e.g. git clone https://github.com/sebastian-stark/PandemicCircles.git from /path/to/folder for this)
3. cd /path/to/folder/PandemicCircles
4. optionally: if you want to modify the ways to be included into the considerations, modify the include_tags and exclude_tags variables in /path/to/folder/PandemicCircles/pandemic_circles.cc appropriately. If a way is tagged with any of the tag/key combinations contained in include_tags, it is included into the considerations, unless the way also involves a tag/key combination contained in exclude_tags. Attention: With the predefined tags, not only paved roads, but also tracks are included. So, the resulting routes are typically not suitable for road bikes. If you want something rideable with a road bike, remove "track" from the include_tags.
5. cmake -DGRAPH_LIB_DIR=~/path/to/graph_lib .
6. make

To run the program:

1. Generate an \*.osm.pbf file with the relevant map information. This file should only cover the region really needed for the generation of the route. Use e.g. Osmium tool (command "osmium extract") to clip a larger \*.osm.pbf file. Put the \*.osm.pbf file into /path/to/folder/PandemicCircles.
2. cd /path/to/folder/PandemicCircles
3. Use command "./pandemic_circles lat=00.00000 lon=00.00000 r=00.00000 file=your_file.osm.pbf svg_output_file=your_svg_file.svg gpx_output_file=your_gpx_file.gpx" to run the program, with appropriate values for the latitude of the center of the circle in degrees (lat), the longitude of the center of the circle in degrees (lon), the radius of the circle in km ( r ), the file name of your \*.osm.pbf file (file), optionally an svg file into which the map and the route are written (svg_output_file), and optionally a gpx file into which the route is written (gpx_output_file).
4. If successful and the respective file names are provided, this will write an svg file /path/to/folder/PandemicCircles/your_svg_file.svg (however, if there are too many edges in the graph, the svg file may become too large for being displayed) and a gpx file /path/to/folder/PandemicCircles/your_gpx_file.gpx


