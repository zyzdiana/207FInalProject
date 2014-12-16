/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <memory>
#include <X11/Xlib.h>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
  double c;        //< damping constant
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K=100;               //< Edge K
  double initial_len;     //< Edge initial length
};
// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint Function object adjust all nodes violating certain constraint
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    //if (n. position () != Point (0 ,0 ,0) && n. position () != Point (1 ,0 ,0)){ 
    if (1){      
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().velocity * dt;
    }
  }

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  constraint(g, t);
  return t + dt;
}
//functor combine two forces
template<typename F1, typename F2>
struct MetaForce{
  F1 f1;
  F2 f2;
  MetaForce(F1 f1_, F2 f2_):f1(f1_), f2(f2_){}
  Point operator()(Node n, double t){
    return f1(n,t)+f2(n,t);
  }
};
template<typename F1, typename F2, typename F3, typename F4>
struct FourForce{
  F1& f1;
  F2& f2;
  F3& f3;
  F4& f4;
  FourForce(F1& f1_, F2& f2_, F3& f3_, F4& f4_):f1(f1_), f2(f2_), f3(f3_), f4(f4_){}
  Point operator()(Node n, double t){
    return f1(n,t)+f2(n,t)+f3(n,t)+f4(n,t);
  }
};
template<typename F1, typename F2>
MetaForce<F1,F2> make_combined_force(F1 f1, F2 f2){
  return MetaForce<F1,F2>(f1,f2);
}
template<typename F1, typename F2, typename F3>
MetaForce<MetaForce<F1,F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3){
  return MetaForce<MetaForce<F1,F2>,F3>( MetaForce<F1,F2>(f1,f2), f3 );
}
template<typename F1, typename F2, typename F3, typename F4>
FourForce<F1,F2,F3,F4> make_combined_force(F1& f1, F2& f2, F3& f3, F4& f4){
  return FourForce<F1,F2,F3,F4>( f1,f2,f3,f4 );
}

//Define Gravity Force
struct GravityForce {
  Point operator()(Node n, double t) {
    Point g(0,0,-n.value().mass*grav);
    return g;
    (void) t;
  }
};
//Define Mass-spring Force
struct MassSpringForce {
  Point operator()(Node n, double t) {
    Point Spring(0,0,0);
    for(auto it=n.edge_begin();it!=n.edge_end();++it){
      auto e=*it;
      Spring=Spring-e.value().K*(e.node1().position()-e.node2().position())/e.length()*(e.length()-e.value().initial_len);
    }
    return Spring;
    (void) t;
  }
};
//Define Damping Force
struct DampingForce {
  Point operator()(Node n, double t) {
    return -n.value().velocity*n.value().c;
    (void) t;
  }
};

//functor combine constraints
template<typename C1, typename C2>
struct MetaConstraint{
  C1& c1;
  C2& c2;
  MetaConstraint(C1& c1_, C2& c2_):c1(c1_), c2(c2_){}
  void operator()(GraphType& g, double t){
    c1(g,t);
    c2(g,t);
  }
};

template<typename C1, typename C2>
MetaConstraint<C1,C2> make_combined_constraint(C1& c1, C2& c2){
  return MetaConstraint<C1, C2>(c1, c2);
}
template<typename C1, typename C2, typename C3>
MetaConstraint<C1,MetaConstraint<C2,C3>> make_combined_constraint(C1& c1, C2& c2, C3& c3){
  return MetaConstraint<C1, MetaConstraint<C2,C3>>(c1, MetaConstraint<C2, C3>(c2, c3));
}




//Ground constraint
struct GroundConstraint {
  void operator()(GraphType& g, double t) {
    for(auto it=g.node_begin();it!=g.node_end();++it){
      if((*it).position().z < -10.75){
        (*it).position().z = -10.75;
        (*it).value().velocity.z=0;
      }
    }
    (void) t;
  }
};
//Sphere constraint
struct SphereConstraint{
  double x, y, z;
  SphereConstraint(double x_=0.5, double y_=0.5, double z_=-0.5):x(x_), y(y_), z(z_){}
  void operator()(GraphType& g, double t) {
    Point center=Point(x,y,z);
    for(auto it=g.node_begin();it!=g.node_end();++it){
      if(norm((*it).position()-center)<0.15){
        Point temp=(*it).position()-center;
        (*it).position()=0.15/norm(temp)*temp+center;
        (*it).value().velocity-= dot((*it).value().velocity,temp/norm(temp))*temp/norm(temp);
      }
    }
    (void) t;
  }
};
//Constraint delete nodes if it is inside sphere
struct SphereDeleteConstraint{
  void operator()(GraphType& g, double t) {
    Point center=Point(0.5,0.5,-0.5);
    for(auto it=g.node_begin();it!=g.node_end();){
      if(norm((*it).position()-center)<0.15){
        g.remove_node(it);
      }
      else{
        ++it;
      }
    }
    (void) t;
  }
};



//Define Wind Force
struct WindForce {
  double level;
  double loc;
  WindForce(double l, double lo=-0.2): level(l), loc(lo){};
  Point operator()(Node n, double t) {
    if( n.position().z<loc+0.15 || n.position().z>loc-0.15 ) return Point(0,level,0);
    return Point(0,0,0);
    (void) t;
  }
};

//customize my listener
struct Listener_Wind: public CS207::SDLViewer::Listener{
  WindForce& wind;
  double pre_level;
  Listener_Wind(WindForce& w): wind(w), pre_level(wind.level){}
  void handle(SDL_Event e){
    switch (e.type) {
      case SDL_MOUSEBUTTONDOWN: {
        // Left mouse button is down
        if (e.button.button == SDL_BUTTON_LEFT) {
          std::cout<<"Stop Wind"<<std::endl;
          wind.level=0;
        }
        // Right mouse button is down
        if (e.button.button == SDL_BUTTON_RIGHT) {
          std::cout<<"Resume Wind"<<std::endl;
          wind.level=pre_level;
        }
      } break;

      case SDL_KEYDOWN: {
        // Keyboard 'arrow right' to increase wind
        if (e.key.keysym.sym == SDLK_RIGHT){
          wind.level+=0.001;
          pre_level=wind.level;
        }
        // Keyboard 'arrow left' to decrease wind
        if (e.key.keysym.sym == SDLK_LEFT){
          wind.level-=0.001;
          pre_level=wind.level;
        }
        // Keyboard 'arrow up' to increase location of wind
        if (e.key.keysym.sym == SDLK_UP){
          wind.loc+=0.001;
        }
        // Keyboard 'arrow down' to decrease location of wind
        if (e.key.keysym.sym == SDLK_DOWN){
          wind.loc-=0.001;
        }
      } break;
    }
  }
};

struct Listener_Constraint: public CS207::SDLViewer::Listener{
  SphereConstraint& con;
  Listener_Constraint(SphereConstraint& con_): con(con_){}
  void handle(SDL_Event e){
    switch (e.type) {
      case SDL_KEYDOWN: {
        // Keyboard 'w' to increase z
        if (e.key.keysym.sym == SDLK_w){
          con.z+=0.05;
        }
        // Keyboard 's' to decrease z
        if (e.key.keysym.sym == SDLK_s){
          con.z-=0.05;
        }
        // Keyboard 'd' to increase y
        if (e.key.keysym.sym == SDLK_d){
          con.y+=0.05;
        }
        // Keyboard 'a' to decrease y
        if (e.key.keysym.sym == SDLK_a){
          con.y-=0.05;
        }
      } break;
    }
  }
};






int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);

      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);

      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // Initialization
  for(auto it=graph.node_begin();it!=graph.node_end();++it){
    (*it).value().mass=float(1)/graph.num_nodes();
    (*it).value().c=float(1)/graph.num_nodes();
    (*it).value().velocity=Point(0);
  }
  for(auto it=graph.edge_begin();it!=graph.edge_end();++it){
      (*it).value().initial_len=(*it).length();
  }
  // Customize Forces
  WindForce wind(0.01);
  GravityForce gravity;
  MassSpringForce ms;
  DampingForce damping;
  
  // Customize Constraint
  GroundConstraint ground;
  SphereConstraint sphere;

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  // Initialize Wind Listener
  std::shared_ptr<CS207::SDLViewer::Listener> ptr_wind(new Listener_Wind(wind));
  std::shared_ptr<CS207::SDLViewer::Listener> ptr_con(new Listener_Constraint(sphere));
  viewer.add_listener(ptr_wind);
  viewer.add_listener(ptr_con);
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 5.001;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cerr << "t = " << t << std::endl;
    symp_euler_step(graph, t, dt, make_combined_force(gravity, ms, damping, wind), 
                   make_combined_constraint(ground, sphere));
    viewer.clear();
    node_map.clear();
    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
    //viewer.set_label(t);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }
  //std::cout<<graph.node(1).position()<<std::endl;
  return 0;
}


