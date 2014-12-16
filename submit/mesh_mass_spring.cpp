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
#include "Mesh.hpp"
#include "tet_mesh.hpp"


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
  double K=700;               //< Edge K
  double initial_len;     //< Edge initial length
};

/** Custom structure of data to store with Triangles */
struct TriData
{
Point n; //the outward surface normal vector
};

/** Custom structure of data to store with Tetrahedral */
struct TetrahedralData {
  double initialVolume;   //< Initial Tetrahedral Volume
};


// Define my tetrahedral mesh type
typedef T_Mesh<NodeData, EdgeData, TetrahedralData> T_MeshType;
typedef typename T_MeshType::Node T_Node;
typedef typename T_MeshType::Edge T_Edge;

// Define my mesh type
typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;

// Define my Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type G_Node;
typedef typename GraphType::edge_type G_Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint Function object adjust all nodes violating certain constraint
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports
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
    if (n. position () != Point (0 ,0 ,0) && n. position () != Point (1 ,0 ,0) && 
      n. position () != Point (0,1,0) && n.position()!=Point(1,1,0)){ 
    //if (1){      
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

//structure combines two forces
template<typename F1, typename F2>
struct MetaForce{
  F1& f1;
  F2& f2;
  MetaForce(F1& f1_, F2& f2_):f1(f1_), f2(f2_){}
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t)+f2(n,t);
  }
};

//structure combines three forces
template<typename F1, typename F2, typename F3>
struct ThrForce{
  F1& f1;
  F2& f2;
  F3& f3;
  ThrForce(F1& f1_, F2& f2_, F3& f3_):f1(f1_), f2(f2_), f3(f3_){}
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t)+f2(n,t)+f3(n,t);
  }
};

//structure combines four forces
template<typename F1, typename F2, typename F3, typename F4>
struct FourForce{
  F1& f1;
  F2& f2;
  F3& f3;
  F4& f4;
  FourForce(F1& f1_, F2& f2_, F3& f3_, F4& f4_):f1(f1_), f2(f2_), f3(f3_), f4(f4_){}
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t)+f2(n,t)+f3(n,t)+f4(n,t);
  }
};

//functor combines two forces
template<typename F1, typename F2>
MetaForce<F1,F2> make_combined_force(F1& f1, F2& f2){
  return MetaForce<F1,F2>(f1,f2);
}

//functor combines three forces
template<typename F1, typename F2, typename F3>
ThrForce<F1, F2, F3> make_combined_force(F1& f1, F2& f2, F3& f3){
  return ThrForce<F1, F2, F3>(f1, f2, f3 );
}

//functor combines four forces
template<typename F1, typename F2, typename F3, typename F4>
FourForce<F1,F2,F3,F4> make_combined_force(F1& f1, F2& f2, F3& f3, F4& f4){
  return FourForce<F1,F2,F3,F4>( f1,f2,f3,f4 );
}

//Define Gravity Force
struct GravityForce {
  template<typename NODE>
  Point operator()(NODE n, double) {
    Point g(0,0,-n.value().mass*grav);
    return g;
  }
};

//Define Mass-spring Force
struct MassSpringForce {
  template<typename NODE>
  Point operator()(NODE n, double) {
    Point Spring(0,0,0);
    for(auto it=n.edge_begin();it!=n.edge_end();++it){
      auto e=*it;
      Spring=Spring-e.value().K*(e.node1().position()-e.node2().position())/e.length()*(e.length()-e.value().initial_len);
    }
    return Spring;
  }
};

//Define Damping Force
struct DampingForce {
  template<typename NODE>
  Point operator()(NODE n, double) {
    return -n.value().velocity*n.value().c;
  }
};

//Define Wind Force
struct WindForce {
  double level;
  WindForce(double l): level(l){};
  template<typename NODE>
  Point operator()(NODE, double) {
    return Point(0,0, level);
  }
};


/** Dashpot Force Functor that returns the Dashpot Force
 */
struct DashpotForce {
  double K_;
  double C_;
  /** DashpotForce Constructor.
   * @param[in] K Spring constant in N/m
   * @param[in] K Damping constant in in N*s/m
   */
  DashpotForce(double K = 100, double C=100) : K_(K), C_(C) {}

  /** Calculates Dashpot Force
   * @param[in] n Valid node.
   * @param[in] t Valid time.
   * @return Point object that represents the dashpot force.
   */
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    Point dashpotForce(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      NODE incidentNode = (*it).node2();
      if( (*it).node1().index() == n.index() ) {
        incidentNode = (*it).node2();
      } else {
        incidentNode = (*it).node1();
      }

      double distance = norm(n.position() - incidentNode.position());
      //Spring component
      double spring_comp = K_ * (distance - (*it).value().initial_len);
      //Damping component
      Point c_comp = C_* (n.value().velocity - incidentNode.value().velocity) *
                 (n.position() - incidentNode.position()) / distance;
      Point unitVector = (n.position() - incidentNode.position()) / distance;

      dashpotForce += (-1.0) * (spring_comp + c_comp) * unitVector;
    }
    return dashpotForce;
  }
};

/** Volume Penalty Force Functor that returns the Volume Penalty Force
 */
struct VolumePenaltyForce {
  T_MeshType* m_;
  double K_;
  /** VolumePenaltyForce Constructor.
   * @param[in] g Gravity in m/s^2.
   */
  VolumePenaltyForce(T_MeshType* m, double K) : m_(m), K_(K) {}

  /** Calculates VolmePenalty Force
   * @param[in] n Valid node.
   * @param[in] t Valid time.
   * @return Point object that represents the volume penalty force.
   */
  Point operator()(T_Node n, double t) {
    (void) t;     // silence compiler warnings
    Point volumePenaltyForce(0,0,0);
    auto AdjTetrahedral = n.nodeAdjTetrahedral();

    // Loop through all of the node's adj tetrahedrals
    for(unsigned k = 0; k < AdjTetrahedral.size(); ++k) {
      T_MeshType::Tetrahedral tet = AdjTetrahedral[k];
      // Mass weighted centroid
      Point baryCenter(0,0,0);
      double totalMass = 0;

      // Loop through each of the tet's node to calculate baryCenter
      for(unsigned i = 0; i < NUM_TET_ADJ_TET; ++i) {
        baryCenter += tet.node(i).position() * tet.node(i).value().mass;
        totalMass += tet.node(i).value().mass;
      }
      baryCenter = baryCenter/totalMass;

      Point unitVector = (n.position() - baryCenter)/norm(n.position() - baryCenter);
      volumePenaltyForce += -K_ * (tet.volume() - tet.value().initialVolume) * unitVector;

    }
    return volumePenaltyForce;
  }
};

//Structure combine constraints
template<typename C1, typename C2>
struct MetaConstraint{
  C1& c1;
  C2& c2;
  MetaConstraint(C1& c1_, C2& c2_):c1(c1_), c2(c2_){}
  template<typename MESH>
  void operator()(MESH& g, double t){
    c1(g,t);
    c2(g,t);
  }
};

// Functors combines constraints
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
  template<typename G>
  void operator()(G& g, double) {
    for(auto it=g.node_begin();it!=g.node_end();++it){
      if((*it).position().z < -0.75){
        (*it).position().z = -0.75;
        (*it).value().velocity.z=0;
      }
    }
  }
};

//Sphere constraint
struct SphereConstraint{
  double x, y, z, r;
  SphereConstraint(double x_=0.5, double y_=0.5, double z_=-0.5, double r_=0.15):x(x_), y(y_), z(z_), r(r_){}
  template<typename G>
  void operator()(G& g, double) {
    Point center=Point(x,y,z);
    for(auto it=g.node_begin();it!=g.node_end();++it){
      if(norm((*it).position()-center)<r){
        Point temp=(*it).position()-center;
        (*it).position()=r/norm(temp)*temp+center;
        (*it).value().velocity-= dot((*it).value().velocity,temp/norm(temp))*temp/norm(temp);
      }
    }
  }
};

//Box Constraint
struct BoxConstraint{
  BoxConstraint(double h1, double h2, double left,double right, double back, double front): 
               h_lower(h1),h_upper(h2),l(left),r(right), b(back), f(front){}
  template<typename G>
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto node = (*it);
      if (node.position().z < h_lower){
        node.position().elem[2] = h_lower;
        node.value().velocity.elem[2] = 0;
      }
      if (node.position().z > h_upper){
        node.position().elem[2] = h_upper;
        node.value().velocity.elem[2] = 0;
      }
      if (node.position().y < l){
        node.position().elem[1] = l;
        node.value().velocity.elem[1] = 0;
      }
      if (node.position().y > r){
        node.position().elem[1] = r;
        node.value().velocity.elem[1] = 0;
      }
      if (node.position().x < b){
        node.position().elem[0] = b;
        node.value().velocity.elem[0] = 0;
      }
      if (node.position().x > f){
        node.position().elem[0] = f;
        node.value().velocity.elem[0] = 0;
      }
    }
  }
  double h_lower;
  double h_upper;
  double l;
  double r;
  double b;
  double f;
};

struct NullConstraint {

  /** Null Constraint Setter
   * @param[in] g Valid mesh.
   * @param[in] t Valid time.
   * @return Point object that represents the combination of forces of @a f1_ and @a f2_.
   */
  template<typename MESH>
  void operator()(MESH& m, double t) {
    (void) t;     // silence compiler warnings
    (void) m;
  }
};

// Plane Constraint
struct PlainConstraint{
  Point center;
  Point out;
  PlainConstraint(Point c, Point o):center(c), out(o){}
  template<typename MESH>
  void operator()(MESH& m, double) {
    double D = -out.x*center.x-out.y*center.y-out.z*center.z;
    double temp = out.x*out.x+out.y*out.y+out.z*out.z;
    for(auto it=m.node_begin(); it!=m.node_end();++it){
      Point position = (*it).position();
      if(dot(position-center, out)<0){
        double k = -(out.x*position.x+out.y*position.y+out.z*position.z+D)/temp;
        (*it).position().x += k*out.x;
        (*it).position().y += k*out.y;
        (*it).position().z += k*out.z;
      }
    }
  }
};

// Plate Conatraint that push the trampoline
struct PlateConstraint{
  Point center;
  bool is_on;
  PlateConstraint(Point c):center(c), is_on(false){}
  template<typename MESH>
  void operator()(MESH& g, double) {
    if(is_on){
      for(auto it=g.node_begin();it!=g.node_end();++it){
        if(norm((*it).position()-center)<0.5){
          Point temp=(*it).position()-center;
          (*it).position()=0.5/norm(temp)*temp+center;
          (*it).value().velocity-= dot((*it).value().velocity,temp/norm(temp))*temp/norm(temp);
        }
      }
    }
  }
};

//customize my listener to change the wind
struct Listener_Wind: public CS207::SDLViewer::Listener{
  WindForce& wind;
  double pre_level;
  bool is_pause=false;
  Listener_Wind(WindForce& w): wind(w), pre_level(wind.level){}
  void handle(SDL_Event e){
    switch (e.type) {

      case SDL_KEYDOWN: {
        // Keyboard 'arrow up' to increase wind
        if (e.key.keysym.sym == SDLK_UP){
          wind.level+=0.01;
          pre_level=wind.level;
          std::cout<<"Increase wind"<<std::endl;
        }
        // Keyboard 'arrow down' to decrease wind
        if (e.key.keysym.sym == SDLK_DOWN){
          wind.level-=0.01;
          pre_level=wind.level;
          std::cout<<"Decrease wind"<<std::endl;
        }
        // Keyboard 'space' to pause and resume the wind
        if (e.key.keysym.sym == SDLK_SPACE){
          is_pause ^= true;
          if(is_pause){
            wind.level = 0;
            std::cout<<"Pause wind."<<std::endl;
          } 
          else{
            wind.level = pre_level;
            std::cout<<"Resume wind."<<std::endl;
          } 
        }
      } break;

    }
  }
};

// Listener that add a sphere constraint to the trampoline
struct Listener_Trampoline: public CS207::SDLViewer::Listener{
  PlateConstraint& plate;
  Listener_Trampoline(PlateConstraint& p): plate(p){}
  void handle(SDL_Event e){
    switch (e.type) {
      case SDL_KEYDOWN: {
        // Keyboard 'w' to increase z
        if (e.key.keysym.sym == SDLK_w){
          plate.center.z +=0.025;
        }
        // Keyboard 's' to decrease z
        if (e.key.keysym.sym == SDLK_s){
          plate.center.z -=0.025;
        }
        // Keyboard 'd' to increase y
        if (e.key.keysym.sym == SDLK_d){
          plate.center.y +=0.025;
        }
        // Keyboard 'a' to decrease y
        if (e.key.keysym.sym == SDLK_a){
          plate.center.y -=0.025;
        }
        // Keyboard 'f' to increase x
        if (e.key.keysym.sym == SDLK_f){
          plate.center.z +=0.025;
        }
        // Keyboard 'r' to decrease x
        if (e.key.keysym.sym == SDLK_r){
          plate.center.z +=0.025;
        }
        // Keyboard 'e' to activate and inactivate the plate
        if (e.key.keysym.sym == SDLK_e){
          plate.is_on ^= true;
          if(plate.is_on){
            std::cout<<"Get set! Plate reset to (0.5, 0.5, 0.65)"<<std::endl;
            plate.center=Point(0.5, 0.5, 0.65);
          } 
          else std::cout<<"Go!!"<<std::endl;
        }
      } break;
    }
  }
};

// Listener that pause or change the speed of simulation
struct Listener_Pause: public CS207::SDLViewer::Listener{
  double& dt;
  bool &is_pause;
  Listener_Pause(double& t, bool& p): dt(t), is_pause(p){}
  void handle(SDL_Event e){
    switch (e.type) {
      case SDL_KEYDOWN: {
        // Keyboard 'p' to pause and resume
        if (e.key.keysym.sym == SDLK_p){
          is_pause ^= true;
          if(is_pause) std::cout<<"Pause simulation."<<std::endl;
          else std::cout<<"Simulation resume."<<std::endl;
        }
        // Keyboard 'o' to speed up
        if (e.key.keysym.sym == SDLK_o){
          dt +=0.0001;
          std::cout<<"Please don't speed up too much, or they will blow up."<<std::endl;
        }
        // Keyboard 'l' to slow down
        if (e.key.keysym.sym == SDLK_l){
          dt -= 0.0001;
          if(dt<=0){
            std::cout<<"Warning: Minimum dt!!"<<std::endl;
            dt=0;
          }
        }
      } break;
    }
  }
};

/** Change position of two ball when they collide
 *  Simulate a plane between two balls
 * @param[in]     m1, m2   two mesh ball
 * @param[in]     r1, r2   radius of ball @a m1 and @a m2
 * @param[in]     c1, c2   center of ball @a m1 and @a m2
 * @param[in]     t        time of this collision
 *
 */
template<typename MESH>
void collision(MESH& m1, MESH& m2, double r1, double r2, Point c1, Point c2, double t){
  if(norm(c1-c2)<(r1+r2)){
    PlainConstraint((c1+c2)/2, c1-c2)(m1, t);
    PlainConstraint((c1+c2)/2, c2-c1)(m2, t);
    //SphereConstraint(c1.x, c1.y, c1.z, r1)(m2, t);
    //SphereConstraint(c2.x, c2.y, c2.z, r2)(m1, t);
  }
}

/** Change position velocity of ball and trampoline when they touch each other
 *  
 * @param[in]     b   a mesh ball
 * @param[in]     s   a mesh sheet
 * @param[in]     t   time of this collision
 * @param[in]     dt  duration of each time step
 * @param[in]     c   center of @a b
 * @param[in]     r   radius of @a bnstant
 *
 */
template<typename MESH1, typename MESH2>
void interact(MESH1& b, MESH2& s, double t, double dt, Point c, double r=0.2) {
  for (auto it = s.node_begin(); it != s.node_end(); ++it) {
    auto n = *it;
    if (norm(n.position() - c) < r) {
      Point R = n.position() - c;
      assert(norm(R) > 0);
      R /= norm(R);
      // R is a unit vector pointing from the center of the
      // sphere towards our point

      // Mass-spring force that the point in the sheet experiences
      Point f = Point(0,0,0);
      for (auto it2 = n.edge_begin(); it2 != n.edge_end(); ++it2) {
        auto e = *it2;
        auto n2 = e.node2();
        f += -e.value().K * (n.position() - n2.position()) / e.length() * (e.length() - e.value().initial_len);
      }
      // Add gravity force
      //f += Point(0,0,-g*n.value().mass);

      // f2 is the component of f pointing in the radial direction
      Point f2 = R*dot(f, R);

      // use f2 to update the ball velocity
      Point b_velocity = Point(0);
      for (auto it_b = b.node_begin(); it_b != b.node_end(); ++it_b){
        (*it_b).value().velocity += f2 * (dt / (*it_b).value().mass);
        b_velocity += (*it_b).value().velocity;
      }
      b_velocity /= b.num_nodes();

      // set the Node's position to the closest point on the radius of the ball
      n.position() = c + R*r;
      // set the Node's velocity to have radial component equal to the ball
      n.value().velocity -= R*dot(n.value().velocity, R);
      n.value().velocity += R*dot(b_velocity, R);
    }
  }
  (void) t;
}

  
struct makecolor{
  float h, s, v;
  makecolor(float h_, float s_, float v_):h(h_), s(s_), v(v_){}
  template<typename NODE>
  /** Construct a color from hue, saturation, and value.
   * h hue (0 = red, 0.166 = yellow, 0.333 = green, 0.5 = cyan,
   *                   0.666 = blue, 0.833 = magenta)
   * s saturation (0 = white, 1 = saturated)
   * v value (0 = dark,  1 = bright)
   * @pre 0 <= @a h, @a s, @a v <= 1 */
  CS207::Color operator() (NODE&){
    //float f=float(node.value())/max;
    return CS207::Color::make_hsv(h, s, v);
  }
};




int main() {

  std::string MESH1_NODE_PATH = "data/grid1.nodes"; 
  std::string MESH1_EDGE_PATH = "data/grid1.tets";
  std::string MESH2_NODE_PATH = "data/sphere200.nodes"; 
  std::string MESH2_EDGE_PATH = "data/sphere200.tris"; 
  std::string MESH3_NODE_PATH = "data/sphere676.nodes"; 
  std::string MESH3_EDGE_PATH = "data/sphere676.tets"; 
  //std::string MESH5_NODE_PATH = "data/large.nodes";
  //std::string MESH5_EDGE_PATH = "data/large.tets";

  // Construct a mesh to store trampoline
  GraphType mesh_1;

  // Construct a mesh to store bounding box and fixed sphere
  MeshType mesh_2;

  // Construct a mesh to store two bounding balls
  T_MeshType mesh_3, mesh_4;

  // Set radius of balls
  double radius=0.16;

  // Create a g_nodes_file from the first input argument
  std::ifstream g_nodes_file(MESH1_NODE_PATH);
  // Interpret each line of the g_nodes_file as a 3D Point and add to the Graph
  std::vector<G_Node> g_nodes;
  Point p;
  while (CS207::getline_parsed(g_nodes_file, p))
    g_nodes.push_back(mesh_1.add_node(p, NodeData()));

  // Create a tets_file from the second input argument
  std::ifstream g_tets_file(MESH1_EDGE_PATH);
  // Interpret each line of the g_tets_file as four ints which refer to g_nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(g_tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      mesh_1.add_edge(g_nodes[t[0]], g_nodes[t[1]]);
      mesh_1.add_edge(g_nodes[t[0]], g_nodes[t[2]]);

      mesh_1.add_edge(g_nodes[t[0]], g_nodes[t[3]]);
      mesh_1.add_edge(g_nodes[t[1]], g_nodes[t[2]]);

      mesh_1.add_edge(g_nodes[t[1]], g_nodes[t[3]]);
      mesh_1.add_edge(g_nodes[t[2]], g_nodes[t[3]]);
    }
  }
  

  // Create a m_nodes_file from the third input argument
  std::ifstream m_nodes_file(MESH2_NODE_PATH);
  // Interpret each line of the m_nodes_file as a 3D Point and add to the Mesh
  std::vector<Node> m_nodes;
  while (CS207::getline_parsed(m_nodes_file, p))
    m_nodes.push_back(mesh_2.add_node(p*0.15+Point(0.5,0.5,1.0)));

  std::ifstream tris_file(MESH2_EDGE_PATH);
  std::array<int,3> t1;
  while (CS207::getline_parsed(tris_file, t1)) {
    mesh_2.add_triangle(m_nodes[t1[0]], m_nodes[t1[1]], m_nodes[t1[2]]);
  }

  //Declare mesh
  //T_MeshType mesh_3, mesh_4;
  std::vector<typename T_MeshType::node_type> t_node_3;
  std::vector<typename T_MeshType::node_type> t_node_4;

  // Read all Points and add them to the Mesh
  std::ifstream t_nodes_file(MESH3_NODE_PATH);

  while (CS207::getline_parsed(t_nodes_file, p)) {
    //Add nodes
      t_node_3.push_back(mesh_3.add_node(p*radius+Point(0.50, 0.5, 0.65), NodeData()));
      t_node_4.push_back(mesh_4.add_node(p*radius+Point(0.40, 0.3, 0.31), NodeData()));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tet_file(MESH3_EDGE_PATH);
  std::array<int,4> t2;
  while (CS207::getline_parsed(tet_file, t2)) {
  //Initialize each tetrahedral's initial value to be the average of its nodes
    mesh_3.add_tetrahedral(t_node_3[t2[0]], t_node_3[t2[1]], t_node_3[t2[2]], t_node_3[t2[3]]);
    mesh_4.add_tetrahedral(t_node_4[t2[0]], t_node_4[t2[1]], t_node_4[t2[2]], t_node_4[t2[3]]);
  }

  // Print out the stats
  std::cout << mesh_1.num_nodes() << " " << mesh_1.num_edges() << std::endl;
  std::cout << mesh_2.num_nodes() << " " << mesh_2.num_edges() << " " << mesh_2.num_triangles() << std::endl;
  std::cout << mesh_3.num_nodes() << " "
            << mesh_3.num_edges() << " "
            << mesh_3.num_tetrahedral() << std::endl;
  std::cout << mesh_4.num_nodes() << " "
            << mesh_4.num_edges() << " "
            << mesh_4.num_tetrahedral() << std::endl;

  // Initialization of mass and length
  // Set initial conditions for your nodes, if necessary.
  // Construct Forces/Constraints
  
  // Set mass of Ball and Sheet
  double ball_mass=1.0;
  double sheet_mass=1.0;

  // Initialization for sheet
  for(auto it=mesh_1.node_begin();it!=mesh_1.node_end();++it){
    (*it).value().mass=sheet_mass/mesh_1.num_nodes();
    (*it).value().c=1.0/mesh_1.num_nodes();
    (*it).value().velocity=Point(0);
  }
  for(auto it=mesh_1.edge_begin();it!=mesh_1.edge_end();++it){
    (*it).value().initial_len=(*it).length();
    (*it).value().K = 1500;
  }

  // Initialization for fixed objects
  //set the mass and velocity of each Node
  for (auto it = mesh_2.node_begin(); it != mesh_2.node_end(); ++it){
    (*it).value().mass = float(1)/mesh_2.num_nodes();
    (*it).value().c=float(1)/mesh_2.num_nodes();
    (*it).value().velocity = Point(0, 0, 0);
  }
  //set K and L for each edge
  for (auto it = mesh_2.node_begin(); it != mesh_2.node_end(); ++it){
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
      (*j).value().initial_len = (*j).length();
      (*j).value().K = 3000;
    }
  }

  //Initialization for balls
  //Zero initial velocity
  //Set mass
  for (auto it = mesh_3.node_begin(); it != mesh_3.node_end(); ++it) {
    auto n = *it;
    n.value().velocity = Point(0,0,0);
    n.value().mass = ball_mass/mesh_3.num_nodes();
  }

  //To set rest length for all of the Edges to their initial length
  for (auto ei = mesh_3.edge_begin(); ei != mesh_3.edge_end(); ++ei ) {
    (*ei).value().initial_len = (*ei).length();
  }

  //To set initial Volume for all of the Tetrahedral to their initial Volume
  for (auto ti = mesh_3.tetrahedral_begin(); ti != mesh_3.tetrahedral_end(); ++ti ) {
    (*ti).value().initialVolume = (*ti).volume();
  }

  //Zero initial velocity
  //Set mass
  for (auto it = mesh_4.node_begin(); it != mesh_4.node_end(); ++it) {
    auto n = *it;
    n.value().velocity = Point(0,0,0);
    n.value().mass = ball_mass/mesh_4.num_nodes();
  }

  //To set rest length for all of the Edges to their initial length
  for (auto ei = mesh_4.edge_begin(); ei != mesh_4.edge_end(); ++ei ) {
    (*ei).value().initial_len = (*ei).length();
  }

  //To set initial Volume for all of the Tetrahedral to their initial Volume
  for (auto ti = mesh_4.tetrahedral_begin(); ti != mesh_4.tetrahedral_end(); ++ti ) {
    (*ti).value().initialVolume = (*ti).volume();
  }

  // Draw a box
  Point vertex[8];
  vertex[0] = Point(0.0, 0.0, -0.5);
  vertex[1] = Point(1.0, 0.0, -0.5);
  vertex[2] = Point(1.0, 1.0, -0.5);
  vertex[3] = Point(0.0, 1.0, -0.5);
  vertex[4] = Point(0.0, 0.0, 1.5);
  vertex[5] = Point(1.0, 0.0, 1.5);
  vertex[6] = Point(1.0, 1.0, 1.5);
  vertex[7] = Point(0.0, 1.0, 1.5);
  Node Vertex_[8];
  for(int i=0; i<8; ++i){
    Vertex_[i] = mesh_2.add_node(vertex[i]);
  }
  mesh_2.add_edge(Vertex_[0], Vertex_[1]);
  mesh_2.add_edge(Vertex_[1], Vertex_[2]);
  mesh_2.add_edge(Vertex_[2], Vertex_[3]);
  mesh_2.add_edge(Vertex_[3], Vertex_[0]);
  mesh_2.add_edge(Vertex_[4], Vertex_[5]);
  mesh_2.add_edge(Vertex_[5], Vertex_[6]);
  mesh_2.add_edge(Vertex_[6], Vertex_[7]);
  mesh_2.add_edge(Vertex_[7], Vertex_[4]);
  mesh_2.add_edge(Vertex_[0], Vertex_[4]);
  mesh_2.add_edge(Vertex_[1], Vertex_[5]);
  mesh_2.add_edge(Vertex_[2], Vertex_[6]);
  mesh_2.add_edge(Vertex_[3], Vertex_[7]);

  // Customize Forces
  WindForce wind(0.001);
  GravityForce gravity;
  MassSpringForce ms;
  DampingForce damping;
  DashpotForce dashpot(1500,0.001); // 1000 0.001
  VolumePenaltyForce volumePenalty(&mesh_3, 1);

  //Final Force for balls and sheet
  auto f_balls = make_combined_force( gravity, dashpot, volumePenalty, wind);
  auto f_sheet = make_combined_force(gravity, ms, damping, wind);

  //Setting Constraints

  //Final Constraints for mesh_3
  //auto c = make_combined_constraints( ground, NullConstraint());
  
  // Customize Constraint
  SphereConstraint sphere(0.5,0.5,1.0);
  PlateConstraint plate(Point(0.5, 0.5, 0.65));

  //auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));
  //Initialize constriants
  auto box = BoxConstraint(-0.5, 1.5, 0.0, 1.0, 0.0, 1.0);

  //Final constraints for balls and sheet;
  auto c_balls = make_combined_constraint(box, sphere);
  auto c_sheet = make_combined_constraint(box, plate);

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  // Initialize Wind Listener
  std::shared_ptr<CS207::SDLViewer::Listener> ptr_wind(new Listener_Wind(wind));
  std::shared_ptr<CS207::SDLViewer::Listener> ptr_tra(new Listener_Trampoline(plate));
  viewer.add_listener(ptr_wind);
  viewer.add_listener(ptr_tra);

  auto graph_node_map = viewer.empty_node_map(mesh_1);
  auto mesh_node_map = viewer.empty_node_map(mesh_2);
  auto t3_node_map = viewer.empty_node_map(mesh_3);
  auto t4_node_map = viewer.empty_node_map(mesh_3);

  viewer.launch();

  viewer.add_nodes(mesh_1.node_begin(), mesh_1.node_end(), graph_node_map);
  viewer.add_edges(mesh_1.edge_begin(), mesh_1.edge_end(), graph_node_map);

  viewer.add_nodes(mesh_2.node_begin(), mesh_2.node_end(), makecolor(0.833, 0.5, 1), mesh_node_map);
  viewer.add_edges(mesh_2.edge_begin(), mesh_2.edge_end(), mesh_node_map);

  viewer.add_nodes(mesh_3.node_begin(), mesh_3.node_end(), t3_node_map);
  viewer.add_edges(mesh_3.edge_begin(), mesh_3.edge_end(), t3_node_map);
  viewer.add_nodes(mesh_4.node_begin(), mesh_4.node_end(), t4_node_map);
  viewer.add_edges(mesh_4.edge_begin(), mesh_4.edge_end(), t4_node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.0005;
  double t_start = 0.0;
  double t_end   = 5.001;
  bool is_pause=false;

  std::shared_ptr<CS207::SDLViewer::Listener> ptr_p(new Listener_Pause(dt, is_pause));
  viewer.add_listener(ptr_p);

  for (double t = t_start; t < t_end; t += dt) {
    //is_pause=true;
    while(1){
      if(!is_pause) break;
      CS207::sleep(0.01);
    }

    symp_euler_step(mesh_1, t, dt, f_sheet , c_sheet);
    symp_euler_step(mesh_3, t, dt, f_balls, c_balls);
    symp_euler_step(mesh_4, t, dt, f_balls, c_balls);

    Point c1(0,0,0), c2(0,0,0);
    for (auto it = mesh_3.node_begin(); it != mesh_3.node_end(); ++it) {
      c1 += (*it).position();
    }
    for (auto it = mesh_4.node_begin(); it != mesh_4.node_end(); ++it) {
      c2 += (*it).position();
    }
    c1 /= mesh_3.num_nodes();
    c2 /= mesh_4.num_nodes();
    collision(mesh_3, mesh_4, radius, radius, c1, c2, t);
    interact(mesh_3, mesh_1, t, dt, c1, radius);
    interact(mesh_4, mesh_1, t, dt, c2, radius);

    // Update viewer with nodes' new positions
    viewer.add_nodes(mesh_1.node_begin(), mesh_1.node_end(), makecolor(0, 0.5, 1), graph_node_map);
    viewer.add_nodes(mesh_3.node_begin(), mesh_3.node_end(), makecolor(0.3, 0.5, 1), t3_node_map);
    viewer.add_nodes(mesh_4.node_begin(), mesh_4.node_end(), makecolor(0.7, 0.5, 1), t4_node_map);

    viewer.set_label(t);
  }
  return 0;
}


