#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template<typename V, typename E>
class Graph {
 private:

  struct internal_edge;
  struct internal_node;
  
 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  typedef V node_value_type;
  typedef E edge_value_type;

  /** Type of this graph. */
  typedef Graph<V, E> graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph(): points_(), i2u_(), edges_(), e_value_(), num_edge(0) {
  }
  /** Default destructor */
  ~Graph(){
    clear();
  }

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points_.clear();
    i2u_.clear();
    edges_.clear();
    e_value_.clear();
    num_edge = 0;
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
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

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return find_element().idx;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      if ((graph_ == x.graph_) && ( uid_ == x.uid_ )) return true;
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
    
    node_value_type& value(){
      return find_element().value;
    }

    const node_value_type& value() const{
      return find_element().value;
    }

    //return node degree
    size_type degree() const{
      valid_node();
      return graph_->edges_[index()].size();
    }

    //Set begin of incidenct iterator
    incident_iterator edge_begin() const{
      valid_node();
      return incident_iterator(graph_, index(), 0);
    }

    //Set end of incident iterator
    incident_iterator edge_end() const{
      valid_node();
      return incident_iterator(graph_, index(), graph_->edges_[index()].size());
    }
    //Check whether node is valid
    void valid_node() const{
      assert(graph_->i2u_[graph_->points_[uid_].idx]==uid_);
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph1, size_type uid):graph_(const_cast<Graph*>(graph1)), uid_(uid){
    }
    //Function to find the position

    internal_node& find_element() const{
      valid_node();
      return graph_->points_[uid_];
    }
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return i2u_.size();
  }

  const Point& node_position(const Node& node) const {
    return node.position();

  }
  

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
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

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
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
    assert( i < size());
    return Node(this, i2u_[i]);
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges 	are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    size_type uid() const {
      return find_node2().uid;
    }

    double length() const{
        return norm(node1().position()-node2().position());
    }

    edge_value_type& value(){
      return graph_->e_value_[uid_];
    }

    const edge_value_type& value() const{
      return graph_->e_value_[uid_];
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      if((graph_ == x.graph_) && ((node1()==x.node1() && node2()==x.node2()) || (node2()==x.node1() && node1()==x.node2()))){
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
      if(graph_<x.graph_) return true;
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

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid1_;
    size_type uid2_;
    size_type uid_;
    Edge(const Graph* graph1, size_type u1, size_type u2, size_type u):graph_(const_cast<Graph*>(graph1)), uid1_(u1), uid2_(u2), uid_(u){
    }
    internal_edge& find_node2() const{
      size_type idx1=graph_->points_[uid1_].idx;
      for(auto it=graph_->edges_[idx1].begin();it!=graph_->edges_[idx1].end();++it){
        if((*it).node2_uid==uid2_) return *it;
      }
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edge;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */ 
  Edge add_edge(const Node& a, const Node& b) {
    if( has_edge(a,b) ){
      for(size_type i = 0; i!=edges_[a.index()].size(); ++i){
        if((edges_[a.index()][i].node2_uid) == b.uid_){
          return Edge(this, a.uid_, b.uid_, edges_[a.index()][i].uid);
        }
      }
    }
    e_value_.push_back(edge_value_type());
    edges_[a.index()].push_back(internal_edge(b.uid_, num_edge));
    edges_[b.index()].push_back(internal_edge(a.uid_, num_edge));
    ++num_edge;
    return Edge(this, a.uid_, b.uid_, num_edge-1);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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

  /** @class Graph::NodeIterator
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
    NodeIterator(const Graph* g, size_type i):graph_(const_cast<Graph*>(g)), it(i){
    }                                     //constructor of valid NodeIterator
    //dereference of NodeIterator
    Node operator*() const{
      assert(it!=graph_->num_nodes());
      return Node(graph_, graph_->i2u_[it]);
    }
    NodeIterator& operator++(){
      ++it;
      return *this;
    }
    bool operator==(const NodeIterator& iterator) const{
      return (it==iterator.it);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type it;    
  };

  //Begin of node iterator
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }
  //End of node iterator
  node_iterator node_end() const{
    return node_iterator(this, i2u_.size());
  }

  /** @class Graph::EdgeIterator
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

    EdgeIterator(const Graph* g, size_type i1, size_type i2):graph_(const_cast<Graph*>(g)), it1(i1), it2(i2){
    }                              //constructor of valid EdgeIterator
        //dereference of EdgeIterator
    Edge operator*() const{
      return Edge(graph_, graph_->i2u_[it1], graph_->edges_[it1][it2].node2_uid, graph_->edges_[it1][it2].uid);
    }
    EdgeIterator& operator++(){
      ++it2;                              
      if (it2==graph_->edges_[it1].size()){              //check whether it is the end of edges_[it1]
        ++it1;
        it2=0;
      }
      while(it1<graph_->edges_.size() && graph_->edges_[it1].size()==0){   //check whether the size of edges_[it1] is zero. 
	++it1;                                                             //If yes, go to next node until edges_[it1] is not empty
      }
      while(it1<graph_->edges_.size() &&  graph_->points_[graph_->edges_[it1][it2].node2_uid].idx<it1){ // skip duplicate by only iterating 
       	++it2;                                                                 //edges with node1().index() smaller than node2().index()
        //do checking above again
        if (it2==graph_->edges_[it1].size()){
          ++it1;
          it2=0;
        }
        while(it1<graph_->edges_.size() && graph_->edges_[it1].size()==0){
	  ++it1;
        }
      }
      return *this;
    }
    bool operator==(const EdgeIterator& iterator) const{
      return ((it1==iterator.it1)&&(it2==iterator.it2));
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type it1;
    size_type it2;
  };


  //set begin of edge iterator
  edge_iterator edge_begin() const{      
    size_type i=0;
    while(edges_[i].size()==0){           //find first node that is not isolated
      ++i;
    }
    return edge_iterator(this, i, 0);
  }

  //set end of edge iterator
  edge_iterator edge_end() const{
    return edge_iterator(this, edges_.size(), 0);
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
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

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    IncidentIterator(const Graph* g, size_type i1, size_type i2):graph_(const_cast<Graph*>(g)), it1(i1), it2(i2){
    }                                           //constructor of valid EdgeIterator

    Edge operator*() const{
      return Edge(graph_, graph_->i2u_[it1], graph_->edges_[it1][it2].node2_uid, graph_->edges_[it1][it2].uid);
    }
    //dereference of IncidentIterator
    IncidentIterator& operator++(){
      assert(it2<graph_->edges_[it1].size());
      ++it2;
      return *this;
    }

    bool operator==(const IncidentIterator& iterator) const{
      return ((it1==iterator.it1)&&(it2==iterator.it2));
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type it1, it2;
  };

//Define remove funtions

/** Remove a node from the graph.
  * @param[in] @a n Node to be removed
  * @return 1 if old has_node(@a n), 0 otherwise
  *
  * @post new size() == old size() - result.
  * @post new num_edges() == old num_edges() - (old @a n.degree())
  *
  * Can invalidate outstanding iterators. 
  * If old has_node(@a n), then @a n becomes invalid, as do any
  * other Node objects equal to @a n. All other Node objects remain valid.
  *
  * Complexity: n.degree()*d, where d=maximum degree().
  */
  size_type remove_node ( const Node & n){
    if(!has_node(n)) return 0;
    size_type n_idx = n.index();
    for(auto it=edges_[n_idx].begin();it!=edges_[n_idx].end();++it){
      size_type n2_idx=points_[(*it).node2_uid].idx;
      for(auto it=edges_[n2_idx].begin();it!=edges_[n2_idx].end();++it){
        if((*it).node2_uid==n.uid_){
	  edges_[n2_idx].erase(it);
	  --num_edge;
	  break;
        }
      }
    }
    edges_.erase(edges_.begin()+n_idx);
    i2u_.erase(i2u_.begin()+n_idx);
    for(auto it=points_.begin();it!=points_.end();++it){
      if((*it).idx>n_idx){
	--(*it).idx;
      }
    }
    return 1;
  }

/** Remove a node from the graph.
  * @param[in] @a n_it NodeIterator of node to be removed
  * @return @a n_it if old @a n_it==node_end()
  *         result=++(old @a n_it) otherwise
  *
  * @post new size() == old size() - old remove_node(*n_it).
  * @post new num_edges() == old num_edges() - (old @a *(n_it).degree())
  *
  * Can invalidate outstanding iterators including old @a n_it. 
  * If old has_node(@a *n_it), then any other Node objects 
  * equal to @a *n_it becomes invalid. All other Node objects remain valid.
  *
  * Complexity: n.degree()*d, where d=maximum degree().
  */
  node_iterator remove_node ( node_iterator n_it ){
    if(n_it!=node_end()){ 
      remove_node(*n_it);
    }
    return n_it;
  }

/** Remove an edge from the graph.
  * @param[in] @a n1 and @a n2 Nodes of an edge to remove
  * @return 1 if old has_edge(@a n1,@a n2), 0 otherwise
  *
  * @post new num_edges() == old num_edges() - 2*result
  *
  * Can invalidate outstanding iterators. 
  * If old has_edge(@a n1,@a n2), then any Edge objects uid_
  * equals to this edge uid_. All other Edge objects remain valid.
  *
  * Complexity: O(num_edges()).
  */
  size_type remove_edge ( const Node & n1, const Node & n2){
    if(!has_edge(n1,n2)) return 0;
    for(auto it=edges_[n1.index()].begin();it!=edges_[n1.index()].end();++it){
      if((*it).node2_uid==n2.uid_){
	edges_[n1.index()].erase(it);
	break;
      }
    }
    for(auto it=edges_[n2.index()].begin();it!=edges_[n2.index()].end();++it){
      if((*it).node2_uid==n1.uid_){
	edges_[n2.index()].erase(it);
	break;
      }
    }
    --num_edge;
    return 0;
  }
/** Remove an edge from the graph.
  * @param[in] @a e an edge to remove
  * @return 1 if old has_edge(@a e.node1(),@a e.node2()), 0 otherwise
  *
  * @post new num_edges() == old num_edges() - 2*old has_edge(@a e.node1(),@a e.node2())
  *
  * Can invalidate outstanding iterators. 
  * If old has_edge(@a e.node1(),@a e.node2()), then any Edge objects with uid_
  * equals to this @a e.uid_ is invalid. All other Edge objects remain valid.
  *
  * Complexity: n.degree()*d, where d=maximum degree().
  */
  size_type remove_edge ( const Edge &e){
    return remove_edge(e.node1(), e.node2());
  }
/** Remove an edge from the graph.
  * @param[in] @a e_it an EdgeIterator edge to remove
  * @return @a e_it if old @a e_it==edge_end()
  *         result==++(old @a e_it) otherwise
  *
  * @post new num_edges() == old num_edges() - 2*old has_edge(@a *e_it.node1(),@a *e_it.node2())
  *
  * Can invalidate outstanding iterators including old @a e_it. 
  * If old has_edge(@a *e_it.node1(),@a *e_it.node2()), then any other Edge objects 
  * equal to @a *e_it becomes invalid. All other Edge objects remain valid.
  *
  * Complexity: n.degree()*d, where d=maximum degree().
  */
  edge_iterator remove_edge ( edge_iterator e_it ){
    if(e_it==edge_end()){
      remove_edge(*e_it);
    }
    return e_it;
  }




 private:

  struct internal_edge{
    size_type node2_uid;
    size_type uid;
    internal_edge(){}
    internal_edge(size_type n, size_type u): node2_uid(n), uid(u){
      
    }
  };

  struct internal_node{
    Point pst;
    node_value_type value;
    size_type idx;
    internal_node(const Point& p, const node_value_type& v, const size_type i): pst(p), value(v), idx(i){}
  };
  //vector stores nodes
  std::vector<internal_node> points_;
  //vector maps index to uid
  std::vector<size_type> i2u_;
  //vector stores edges
  std::vector< std::vector<internal_edge> > edges_;
  //vector stores value of edges by uid
  std::vector<edge_value_type> e_value_;
  //number of edges
  size_type num_edge;
  

};

#endif
