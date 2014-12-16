//Team Members: Ruitao Du, Yingzhuo Zhang
#pragma once
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */
#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"
/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {
 private:
  struct internal_edge;
  struct internal_node;
  struct Edge_value;
  struct Tri_value;


 public:
 
   /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////
  
  typedef N node_value_type;
  typedef E edge_value_type;
  typedef T tri_value_type;

  /** Type of this mesh. */
  typedef Mesh<N,E,T> mesh_type;
  
  /** Type of indexes and sizes.
   *  Return type of Mesh::Node::index(), Mesh::Node::degree(),
   *  Mesh::num_nodes(), Mesh::num_edges(),
   *  and argument type of Mesh::node(size_type) */
  typedef unsigned size_type;


  /** Predeclarations of Node, Edge, Triangle type. */
  class Node;
  class Edge;
  class Triangle;
  
  /** Synonym for Node, Edge, Triangle (following STL conventions). */
  typedef Node node_type;
  typedef Edge edge_type;
  typedef Triangle tri_type;
  
  /** Type of node iterators, which iterate over all mesh nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all mesh edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of triangle iterators, which iterate over all mesh triangles. */
  class TriIterator;
  /** Synonym for EdgeIterator */
  typedef TriIterator tri_iterator;

  /** Type of edge incident iterators, which iterate incident edges to a node. */
  class EdgeIncidentIterator;
  /** Synonym for EdgeIncidentIterator */
  typedef EdgeIncidentIterator edge_incident_iterator;

  /** Type of triangle incident iterators, which iterate incident edges to a node. */
  class TriIncidentIterator;
  /** Synonym for TriIncidentIterator */
  typedef TriIncidentIterator tri_incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////
  
  /** Construct an empty mesh. */
  Mesh(): points_(), i2u_(), edges_(), e_value_(), num_edge(0), triangles_() {
  }
  /** Default destructor */
  ~Mesh(){
    clear();
  }
  
  /** Remove all nodes, edges and triangles from this mesh.
   * @post num_nodes() == 0 && num_edges() == 0 && num_triangles()==0
   *
   * Invalidates all outstanding Node, Edge and Triangle objects.
   */
  void clear() {
    points_.clear();
    i2u_.clear();
    edges_.clear();
    e_value_.clear();
    num_edge = 0;
    triangles_.clear();
  }


  /////////////////
  // MESH NODES //
  /////////////////

  /** @class Mesh::Node
   * @brief Class representing the mesh's nodes.
   *
   * Node objects are used to access information about the Mesh's nodes.
   */
  class Node: private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Mesh class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Mesh::node_type x;
     * if (...should pick the first node...)
     *   x = mesh.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return find_element().pst;
    }

    /** To modify this node's position. */
    Point& position(){
      return find_element().pst;
    }

    /** Return this node's index, a number in the range [0, mesh_size). */
    size_type index() const {
      return find_element().idx;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same mesh and the same index.
     */
    bool operator==(const Node& x) const {
      if ((mesh_ == x.mesh_) && ( uid_ == x.uid_ )) return true;
      return false;
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
      if (uid_ < x.uid_) return true;
      return false;
    }
    
    /** Return the value associated with this Node which is of type node_value_type
     * Complexity: O(1)
     */
    node_value_type& value(){
      return find_element().value;
    }

    /** To modify the value associated with this Node which is of type node_value_type
     * Complexity: O(1)
     */
    const node_value_type& value() const{
      return find_element().value;
    }

    /** Return the degree of this Node, which is how many incident triangles it has
     * @pre this is a valid Node
     * @return the number of incident triangles of type size_type
     * Complexity: O(1)
     */
    size_type degree() const{
      valid_node();
      return mesh_->edges_[index()].size();
    }

    /** Return an edge_incident_iterator that points to 
     *  the begin position of vector of incident edges 
     *	Complexity: O(1).
     */
    edge_incident_iterator edge_begin() const{
      valid_node();
      return edge_incident_iterator(mesh_, index(), 0);
    }

    /** Return an edge_incident_iterator that points to 
     *  the end position of vector of incident edges 
     *	Complexity: O(1).
     */
    edge_incident_iterator edge_end() const{
      valid_node();
      return edge_incident_iterator(mesh_, index(), mesh_->edges_[index()].size());
    }

    /** Return a tri_incident_iterator that points to 
     *  to the first adjacent triangle of this Node
     *  Complexity: O(1)
     */
    tri_incident_iterator tri_begin() const{
      valid_node();
      return tri_incident_iterator(mesh_, index(), 0);
    }

    /** Return a tri_incident_iterator that points to 
     *  to the end adjacent triangle of this Node
     *  Complexity: O(1)
     */
    tri_incident_iterator tri_end() const{
      valid_node();
      return tri_incident_iterator(mesh_, index(), find_element().adj_triangles.size());
    }

    /** assertion that checks whether this node is valid*/
    void valid_node() const{
      assert(mesh_->i2u_[mesh_->points_[uid_].idx]==uid_);
    }
   private:
    // Allow Mesh to access Node's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type uid_;
    
    /* private constructor */
    Node(const Mesh* mesh1, size_type uid):mesh_(const_cast<Mesh*>(mesh1)), uid_(uid){
    }
    
    /* Internal function that helps finding the content of this Node
     * Return the content associated with this Node of type internal_node
     */
    internal_node& find_element() const{
      valid_node();
      return mesh_->points_[uid_];
    }
  };

  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return i2u_.size();
  }
  

  /** Add a node to the mesh, returning the added node.
   * @param[in] @a position The new node's position
   * @param[in] @a value The value associated with the node
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Push point to the vector
    size_type u = points_.size();
    size_type i = i2u_.size();
    i2u_.push_back(u);
    points_.push_back(internal_node(position, value, i));
    edges_.push_back(std::vector<internal_edge>());     //initialize memory of edge
    // Returns a Node that points to the new node
    return Node(this, u);
  }

  /** Determine if this Node belongs to this Mesh
   * @return True if @a n is currently a Node of this Mesh
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(node(n.index())==n) return true;
    else{
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert( i < num_nodes());
    return Node(this, i2u_[i]);
  }

  /////////////////
  //  MESH EDGES //
  /////////////////

  /** @class Mesh::Edge
   * @brief Class representing the mesh's edges.
   *
   * Edges  are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(mesh_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(mesh_, uid2_);
    }

	/** Return the  */
    size_type uid() const {
      return find_node2().uid;
    }

    double length() const{
        return norm(node1().position()-node2().position());
    }

    edge_value_type& value(){
      return mesh_->e_value_[uid_].edge_info;
    }

    const edge_value_type& value() const{
      return mesh_->e_value_[uid_].edge_info;
    }

    size_type degree(){
      return mesh_->e_value_[uid_].adj_triangles.size();
    }


    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same mesh and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      if((mesh_ == x.mesh_) && ((node1()==x.node1() && node2()==x.node2()) || (node2()==x.node1() && node1()==x.node2()))){
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
      if(mesh_<x.mesh_) return true;
      if (node1()<x.node1())
        return true;
      else if (node1()==x.node1()){
        if(node2()<x.node2()){
          return true;
        }
        else return false;
      }
      else{
        return false;
      }
    }
    const std::vector<size_type> adj_triangles() const{
      return mesh_->e_value_[uid_].adj_triangles;
    }
   private:
    // Allow Mesh to access Edge's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type uid1_;
    size_type uid2_;
    size_type uid_;
    
    /* private constructor */
    Edge(const Mesh* mesh1, size_type u1, size_type u2, size_type u):mesh_(const_cast<Mesh*>(mesh1)), uid1_(u1), uid2_(u2), uid_(u){
    }
    
    /** Internal function that adds triangle indicies to the vector of adjacency triangles
     */
    size_type add_triangle(const size_type i){
      mesh_->e_value_[uid_].adj_triangles.push_back(i);
      return mesh_->e_value_[uid_].adj_triangles.size();
    }
    
    /* Internal function that helps find the content of this Edge
     * Return the content associated with this Edge of type internal_edge
     */
    internal_edge& find_node2() const{
      size_type idx1=mesh_->points_[uid1_].idx;
      for(auto it=mesh_->edges_[idx1].begin();it!=mesh_->edges_[idx1].end();++it){
        if((*it).node2_uid==uid2_) return *it;
      }
      assert(false);
    }
  };

  /** Return the total number of edges in the mesh.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this mesh
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type size = edges_[a.index()].size();
    for(size_type i = 0; i!=size; ++i){
      if((edges_[a.index()][i].node2_uid) == b.uid_){
        return true;
      }
    }
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert( i < num_edges() );
    auto it = edge_begin();
    for (; it!=edge_end(); ++it){
      if ((*it).uid() == i) return *it;
    }
    assert(false);
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Mesh::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }      
                                   
	/** Dereferencing operator
	 *  Returns the node this iterator points to
	 *  Complexity: O(1)
	 */
    Node operator*() const{
      assert(it!=mesh_->num_nodes());
      return Node(mesh_, mesh_->i2u_[it]);
    }
    
    /** Increment operator
     *  @post this iterator points to the a valid node
     *  @return iterator pointing to the next valid node
     * Complexity: O(1)
     */
    NodeIterator& operator++(){
      ++it;
      return *this;
    }
    
    /** Comparison operator
     *  @param[in] NodeIterator @a iterator to be compared with 
     *  @return true if the two iterators are equal
     *  Complexity: O(1)
     */
    bool operator==(const NodeIterator& iterator) const{
      return (it==iterator.it);
    }

   private:
    // Allow Mesh to access Edge's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type it;
    
    /** private constructor */
    NodeIterator(const Mesh* g, size_type i):mesh_(const_cast<Mesh*>(g)), it(i){
    }        
  };

  /** Returns an iterator that points to the first node of the graph */
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }
  
  /** Returns an iterator that points to the position after
   *  the last valid node of the graph
   */
  node_iterator node_end() const{
    return node_iterator(this, i2u_.size());
  }

  /** @class Mesh::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

	/** Dereference operator
	 *  @return the edge that this iterator points to
	 */
    Edge operator*() const{
      return Edge(mesh_, mesh_->i2u_[it1], mesh_->edges_[it1][it2].node2_uid, mesh_->edges_[it1][it2].uid);
    }
    
    /** Increment operator
     *  @post this iterator points to the next valid edge
     *  @return iterator pointing to the next valid edge
     */
    EdgeIterator& operator++(){
      ++it2;                              
      if (it2==mesh_->edges_[it1].size()){
        //check whether it is the end of edges_[it1]
        ++it1;
        it2=0;
      }
      while(it1<mesh_->edges_.size() && mesh_->edges_[it1].size()==0){
         //check whether the size of edges_[it1] is zero. 
 		 ++it1;                                                             
 		 //If yes, go to next node until edges_[it1] is not empty
      }
      while(it1<mesh_->edges_.size() &&  mesh_->points_[mesh_->edges_[it1][it2].node2_uid].idx<it1){ 
        // skip duplicate by only iterating 
        ++it2;                                                                 
        //edges with node1().index() smaller than node2().index()
        //do the checking above again
        if (it2==mesh_->edges_[it1].size()){
          ++it1;
          it2=0;
        }
        while(it1<mesh_->edges_.size() && mesh_->edges_[it1].size()==0){
    	  ++it1;
        }
      }
      return *this;
    }
    
   /** Comparaison operator
    *  @param[in] EdgeIterator @a iterator to be compared with 
    *  @return true if the two iterators are equal
    */
    bool operator==(const EdgeIterator& iterator) const{
      return ((it1==iterator.it1)&&(it2==iterator.it2));
    }

   private:
    // Allow Mesh to access Edge's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type it1;
    size_type it2;
    
    /** private constructor */
    EdgeIterator(const Mesh* m, size_type i1, size_type i2):
    	mesh_(const_cast<Mesh*>(m)), it1(i1), it2(i2){
    } 
  };


  /** Returns an edge_iterator that points to the first edge of this graph*/
  edge_iterator edge_begin() const{      
    size_type i=0;
    while(edges_[i].size()==0){           //find first node that is not isolated
      ++i;
    }
    return edge_iterator(this, i, 0);
  }

  /** Returns an edge_iterator that points to the position after the last edge of this graph*/
  edge_iterator edge_end() const{
    return edge_iterator(this, edges_.size(), 0);
  }

  /** @class Mesh::EdgeIncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class EdgeIncidentIterator:private totally_ordered<EdgeIncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIncidentIterator. */
    EdgeIncidentIterator() {
    }

	/** Dereference operator
	*   @return the edge that this EdgeIncidentIterator points to
	*/
    Edge operator*() const{
      return Edge(mesh_, mesh_->i2u_[it1], mesh_->edges_[it1][it2].node2_uid, mesh_->edges_[it1][it2].uid);
    }
    
    /** Increment operator
     *  @post this iterator points to the next valid incident edge
     *  @return iterator pointing to the next valid incident edge
     */    
    EdgeIncidentIterator& operator++(){
      assert(it2<(mesh_->edges_[it1].size()));
      ++it2;
      return *this;
    }

   /** Comparaison operator
    *  @param[in] IncidentIterator @a x to be compared with 
    *  @return true if the two iterators are equal
    */
    bool operator==(const EdgeIncidentIterator& iterator) const{
      return ((it1==iterator.it1)&&(it2==iterator.it2));
    }

   private:
    // Allow Mesh to access Edge's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type it1, it2;
    
    /** private constructor */
    EdgeIncidentIterator(const Mesh* m, size_type i1, size_type i2):mesh_(const_cast<Mesh*>(m)), it1(i1), it2(i2){
    } 
  };
  
  ////////////////////
  // MESH TRIANGLES //
  ////////////////////

  /** @class Mesh::Triangle
   *  @brief Class representing the mesh's edges.
   *
   *  Each Triangle is composed of three Nodes and three Edges.
   */
  class Triangle:private totally_ordered<Triangle>{
   public:
    /** Construct an invalid Triangle. */
    Triangle():valid_(false){
    }
    
    /** Return a boolean indicating whether or not this triangle is valid 
     * @ return true if the triangle is valid, false otherwise
     */
    bool valid() const{
      return valid_;
    }
    
    /** Accessing value associated with Triangle
     *  @return the reference of the value
     *
     *  Complexity: O(1) 
     */
    tri_value_type& value(){
      return fetch().tri_info;
    }

    /** Return the value associated with Triangle
     *  @return the reference of the value
     *
     *  Complexity: is O(1) 
     */
    const tri_value_type& value() const{
      return fetch().tri_info;
    }

    /** Return the index of this Triangle, a number in the range [0, num_triangles()]. 
     *  @return The unique id of this Triangle, which is the same as the triangle index.
     * Complexity: O(1)
     */
    size_type index(){
      return t_index;
    }

    /** Return the area of this Triangle
     *  Complexity: O(1)
     */
    double area(){
      Point a = mesh_->node(fetch().tri_nodes[0]).position();
      Point b = mesh_->node(fetch().tri_nodes[1]).position();
      Point c = mesh_->node(fetch().tri_nodes[2]).position();
      
      double area = std::abs(a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y));
      return area;
    }
		
    /** Return the Node of this Triangle corresponding to the given index
     *  @param[in] @a i The index of the Node that needs to be retrieved.
     *  @pre 0 <= i < 3
     *  @post result Node.index() == mesh_->triangles_[t_index].tri_nodes[i]
     *  @note triangleNode(i) is the node opposite to the edge triangleEdge(i)
     *
     *  Complexity: O(1)
     */    
     node_type node(size_type i){
      return mesh_->node(fetch().tri_nodes[i]);
     }
		
    /** Return the Edge of this Triangle corresponding to the given index
     *  @param[in] @a i The index of the Edge that needs to be retrieved.
     *  @pre 0 <= i < 3
     *  @note triangleEdge(i) is the edge opposite to the node triangleNode(i)
     *
     *  Complexity: O(1)
     */   
    edge_type edge(int i){
      return fetch().tri_edges[i];
    }


  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
   
    /** Return the i-th adjacent triangle of this Triangle. 
     * @param[in] @a i The index of adjacent triangle.
     * @pre 0 <= i < 3
     *
     * Complexity: O(1)
     */
    tri_type adj_triangles(size_type i){
      std::vector<size_type> adj = edge(i).adj_triangles();
      if(adj.size()==1) return Triangle();
      if(adj[0]==t_index) return Triangle(mesh_, adj[1]);
      else return Triangle(mesh_, adj[0]);
    }

    /** Return the outward normal vector
     * @param[in] an index of an edge @a i
     * @pre has_edge(e)
     * @post result*(e.node1().position()-e.node2().position())=0
     * @post norm(result)=1
     *
     * Complexity: O(1)
     */
    Point out_normal(int i){
      Point a = edge(i).node1().position();
      Point b = edge(i).node2().position();
      Point c = node(i).position();
      
      Point result = Point();
	  
      Point ac(c.x-a.x, c.y-a.y, 0);
      Point n1(b.y-a.y, a.x-b.x, 0);
      Point n2(a.y-b.y, b.x-a.x, 0);
      
      if(dot(ac,n1) < 0) result = n1;
      else result = n2;

      return result;
    }

   private:
    // Allow Mesh to access Triangle's private member data and functions.
    friend class Mesh;
    Mesh* mesh_;
    size_type t_index;	//Index of triangle
    bool valid_;
    
    /* Internal function that helps finding the content of this Triangle
     * Return the content associated with this Triangle of type internal_node
     */    Tri_value& fetch() const{
      return mesh_->triangles_[t_index];
    }
    
    /* private constructor */
    Triangle(const Mesh* m, size_type index): mesh_(const_cast<Mesh*>(m)), t_index(index), valid_(true){
    }
  };

  /** Add a triangle into the Mesh
   *  @param[in] @n1 Three node positions
   *  @pre has_node(Node @a n1), has_node(Node @a n2), has_node(Node @a n3)
   *  @post old num_tri()=new num_tri()-1
   * 
   *  Complexity: is O(1) 
   */
  tri_type add_triangle(node_type n1, node_type n2, node_type n3){
    edge_type e3=add_edge(n1, n2);
    edge_type e1=add_edge(n2, n3);
    edge_type e2=add_edge(n3, n1);
    size_type tri_index = triangles_.size();
    triangles_.push_back(Tri_value());
    Tri_value& curr = triangles_.back();
    curr.tri_nodes.push_back(n1.index());
    curr.tri_nodes.push_back(n2.index());
    curr.tri_nodes.push_back(n3.index());
    curr.tri_edges.push_back(e1);
    curr.tri_edges.push_back(e2);
    curr.tri_edges.push_back(e3);
    n1.find_element().adj_triangles.push_back(tri_index);
    n2.find_element().adj_triangles.push_back(tri_index);
    n3.find_element().adj_triangles.push_back(tri_index);
    e1.add_triangle(tri_index);
    e2.add_triangle(tri_index);
    e3.add_triangle(tri_index);
    return Triangle(this, tri_index);
  }
   
  /** @class Mesh::TriIterator
   * @brief Iterator class for triangles. A forward iterator. */
  class TriIterator:private totally_ordered<TriIterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;
    
    /** Construct an invalid NodeIterator. */
    TriIterator() {
    }
    
	/** Dereferencing operator
	 *  Returns the Triangle this iterator points to
	 *  Complexity: O(1)
	 */
    Triangle operator*() const{
      assert(tri_index<(mesh_->num_triangles()));
      return Triangle(mesh_, tri_index);
    }
    
    /** Increment operator
     * @post the iterator moves to the next triangle
     * @return the new iterator that points to the next triangle
     * Complexity: O(1)
     */
    TriIterator operator++(){
      ++tri_index;
      //assert(tri_index<=mesh_->num_triangles());
      return *this;
    }
    
    /** Comparison operator
     *  @param[in] @a t a triangle iterator
     *  @return true if they are equal
     *  Complexity: O(1)
     */
    bool operator==(const tri_iterator& t) const{
      assert(t.mesh_==mesh_);
      return t.tri_index==tri_index;
    }
   private:
    friend class Mesh;
    Mesh* mesh_;
    size_type tri_index;
    
    /* private constructor */
    TriIterator(const Mesh* m, size_type i):mesh_(const_cast<Mesh*>(m)), tri_index(i){
    }
  };
  
  /** Return the the first adjacent triangle 
   * @return Tri_iterator pointing to the first adjacent triangle of a node
   *
   * Complexity: O(1)
   */
  TriIterator tri_begin(){
    return TriIterator(this, 0);
  }

  /** Return the the lsat adjacent triangle 
   * @return tri_iterator pointing to the last adjacent triangle of a node
   *
   * Complexity O(1)
   */
  TriIterator tri_end(){
    return TriIterator(this, triangles_.size());
  }
	
  /** Return the triangle with index @a i.
   * @pre 0 <= @a i < num_triangle()
   * @post result_triangle.index() == i
   *
   * Complexity: O(1).
   */
  Triangle triangle(size_type i){
    return Triangle(this, i);
  }

  /** Return the number of triangles in the mesh. */
  size_type num_triangles() const {
    return triangles_.size();
  }

  /** @class Mesh::TriIncidentIterator
   * @brief Incident iterator class for triangles. */
  class TriIncidentIterator:private totally_ordered<TriIncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid TriIncidentIterator. */
    TriIncidentIterator() {
    }

	/** Dereferencing operator
	 *  Returns the Triangle this iterator points to
	 *  Complexity: O(1)
	 */
    Triangle operator*() const{
      return Triangle(mesh_, mesh_->points_[mesh_->i2u_[node_index]].adj_triangles[it]);
    }
    
    /** Increment operator
     * @post the incident iterator moves to the next valid incident triangle
     * @return the new iterator that points to the next valid incident triangle
     * Complexity: O(1)
     */    
     TriIncidentIterator& operator++(){
      ++it;
      return *this;
    }
    
   /**Comparation operator
    * @param[in] TriIncidentIterator @a iterator to be compared with 
    * @return true if the two iterators are equal
    */
    bool operator==(const TriIncidentIterator& iterator) const{
      return it==iterator.it;
    }

   private:
    friend class Mesh;
    Mesh* mesh_;
    size_type node_index;
    size_type it;
    
    /** private constructor */
    TriIncidentIterator(const Mesh* m, size_type n, size_type i):mesh_(const_cast<Mesh*>(m)), node_index(n), it(i){
    }  
  };

    /** Add an edge to the mesh, or return the current edge if it already exists.
   * @param[in] @a a, @a b two nodes between which an edge will be added
   * @pre @a a and @a b are distinct valid nodes of this mesh
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */ 
  Edge add_edge(const Node& a, const Node& b) {
    if( has_edge(a,b) ){
      for(size_type i = 0; i!=edges_[a.index()].size(); ++i){
        if((edges_[a.index()][i].node2_uid) == b.uid_){
          return Edge(this, a.uid_, b.uid_, edges_[a.index()][i].uid);
        }
      }
    }
    e_value_.push_back(Edge_value());
    edges_[a.index()].push_back(internal_edge(b.uid_, num_edge));
    edges_[b.index()].push_back(internal_edge(a.uid_, num_edge));
    ++num_edge;
    return Edge(this, a.uid_, b.uid_, num_edge-1);
  }


 private:
  struct internal_edge{
    size_type node2_uid;
    size_type uid;
    std::vector<size_type> adj_triangles;
    internal_edge(){}
    internal_edge(size_type n, size_type u): node2_uid(n), uid(u), adj_triangles(){
    }
  };

  struct internal_node{
    Point pst;
    node_value_type value;
    size_type idx;
    std::vector<size_type> adj_triangles;
    internal_node(const Point& p, const node_value_type& v, const size_type i): pst(p), value(v), idx(i), adj_triangles(){}
  };

  struct Edge_value{
    edge_value_type edge_info;
    std::vector<size_type> adj_triangles;
    Edge_value(): edge_info(), adj_triangles(){}
  };
  struct Tri_value{
    tri_value_type tri_info;
    std::vector<size_type> tri_nodes;
    std::vector<edge_type> tri_edges;
  };



  //vector that stores nodes
  std::vector<internal_node> points_;
  //vector that maps indicies to uids
  std::vector<size_type> i2u_;
  //vector that stores edges
  std::vector< std::vector<internal_edge> > edges_;
  //vector that stores value of edges by uid
  std::vector<Edge_value> e_value_;
  //number of edges
  size_type num_edge;
  //vector that stores triangles
  std::vector<Tri_value> triangles_;
};
