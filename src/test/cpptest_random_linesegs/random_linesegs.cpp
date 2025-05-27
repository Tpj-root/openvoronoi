
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "voronoidiagram.hpp"
#include "version.hpp"
#include "common/point.hpp"
#include "utility/vd2svg.hpp"

#include <boost/random.hpp>
#include <boost/timer/timer.hpp>
//#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

typedef std::pair<ovd::Point,ovd::Point> segment;

// return true if segment s1 intersects with segment s2
bool intersects( segment& s1, segment& s2 ) {
    ovd::Point p1 = s1.first;
    ovd::Point p2 = s1.second;
    ovd::Point p = p1;
    ovd::Point r = p2-p1;
    ovd::Point q1 = s2.first;
    ovd::Point q2 = s2.second;
    ovd::Point q = q1;
    ovd::Point s = q2-q1;
    // t = (q-p) cross (s) / (r cross s)
    // u = (q-p) cross (r) / (r cross s)
    if ( r.cross(s) == 0 ) { //parallel lines
        if ( (q-p).cross(r) == 0 ) //collinear
            return true;
        else
            return false; // parallel lines that never intersect
    }
    double t = (q-p).cross(s) / (r.cross(s));
    double u = (q-p).cross(r) / (r.cross(s));
    if ( (0.0<=t) && (t<=1.0) && (0.0<=u) && (u<=1.0) )
        return true;
    return false;
}

// test if s intersects with any of the segments in segs
// return true if there is an intersection, otherwise false
bool segment_intersects(std::vector<segment>& segs, segment& s) {    
    BOOST_FOREACH(segment& seg, segs) {
        if (intersects(seg,s))
            return true;
    }
    return false; // no intersections found
}

// create a random segment
segment random_segment(double far, double r1,double r2,double r3,double r4) {
    double pradius = (1.0/sqrt(2))*far;
    segment out;
    double x1 = -pradius+2*pradius*r1;
    double y1 = -pradius+2*pradius*r2;
    double x2 = -pradius+2*pradius*r3;
    double y2 = -pradius+2*pradius*r4;
    out.first = ovd::Point(x1,y1);
    out.second = ovd::Point(x2,y2);
    return out;
}

// create Nmax random segments with a far-radius circle
std::vector<segment> random_segments(double far,int Nmax) {
    boost::mt19937 rng(42);
    boost::uniform_01<boost::mt19937> rnd(rng);
    
    std::vector<segment> segs; 
    for (int n=0;n<Nmax;n++) { 
        segment seg;
        do {
            seg = random_segment(far, rnd(),rnd(),rnd(),rnd()); // try a new segment until one is found
        } while (segment_intersects(segs, seg)); // that does not intersect with any previous segment.
        // NOTE: for high Nmax this can be very slow
        segs.push_back(seg);
    }
    return segs;
}

// random points
int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("This program calculates the voronoi diagram for n random non-intersecting line-segments\nAllowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<int>(), "set number of line-segments")
        ("d",  "run in debug-mode")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
    unsigned int nmax = 100;

    if (vm.count("n")) 
        nmax = vm["n"].as<int>();
    else {
        std::cout << desc << "\n";
        return 1;
    }
    
    
    

    std::cout << "Number of random line segments: " << nmax << "\n";
    int bins = (int)sqrt(nmax);
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,10*bins);
    
    if (vm.count("d")) {
        std::cout << "running in debug mode!\n";
        vd->debug_on();
    }
    
    std::cout << "OpenVoronoi version: " << ovd::version() << "\n";
    
    boost::mt19937 rng(42);
    boost::uniform_01<boost::mt19937> rnd(rng);
    std::cout << "Waiting for " << nmax << " random line segments..." << std::flush;
    
    boost::timer::cpu_timer tmr;

    // Generate segments and time it
    std::vector<segment> segs = random_segments(1, nmax);
    auto elapsed = tmr.elapsed();
    double t_points = elapsed.wall / 1e9;  // timing for segment generation
    std::cout << "done in " << t_points << " seconds\n" << std::flush;
    std::cout << "number of segs: " << segs.size() << "\n" << std::flush;

    typedef std::pair<int,int> IdSeg;
    typedef std::vector<IdSeg> IdSegments;
    IdSegments segment_ids;

    // Insert point sites timing
    tmr.start();
    BOOST_FOREACH(segment s, segs) {
        IdSeg id;
        id.first = vd->insert_point_site(s.first);
        id.second = vd->insert_point_site(s.second);
        segment_ids.push_back(id);
    }
    elapsed = tmr.elapsed();
    t_points = elapsed.wall / 1e9;  // reusing t_points for point insertion time

    std::cout << "all point-sites inserted!\n" << std::flush;

    // Insert line sites timing
    tmr.start();
    BOOST_FOREACH(IdSeg id, segment_ids) {
        vd->insert_line_site(id.first, id.second);
    }
    elapsed = tmr.elapsed();
    double t_lines = elapsed.wall / 1e9;

    std::cout << "all line-sites inserted!\n" << std::flush;

    std::cout << "Points: " << t_points << " seconds \n";
    std::cout << "Lines: " << t_lines << " seconds \n";

    double norm = 2 * nmax * log(2 * (double)nmax) / log(2.0);
    std::cout << "Points: " << 1e6 * t_points / norm << " us * n*log2(n)\n";
    std::cout << "Lines: " << 1e6 * t_lines / norm << " us * n*log2(n)\n";

    std::cout << vd->print();
    vd2svg("random_segments.svg", vd);

    delete vd;
    return 0;
}

