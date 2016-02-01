%module RCSP

%{
  #include "../include/CommonTypes.h"
%}
%include "../include/CommonTypes.h"

// Rename comparison operators so they can be used seamlessly from python
%rename(__ne__) operator!=;
%rename(__eq__) operator==;
%rename(__lt__) operator<;
%rename(__le__) operator<=;
%rename(__gt__) operator>;
%rename(__ge__) operator>=;


// Rename functions so they match old api
%rename(shortest_path_between_rois) r2r_shortest_path;

// Make it possible to use len(Property)
%rename(__len__) dim; 

%include "exception.i"
%catches(std::out_of_range) __getitem__;
%catches(std::out_of_range) __setitem__;

%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_string.i"

// Get the typedefs
%{
  #include "CommonTypes.h"
%}

using namespace BrainGraph;
%template(IdVector) std::vector<id_type>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;

// Property
%{
  #include "../include/Property.h"
%}
%include "../include/Property.h"

%extend BrainGraph::Property {
  double __getitem__(size_t i) {
    return $self->operator[](i);
  }
  void __setitem__(size_t i, double value) {
    $self->operator[](i) = value;
  }
}


// Node 
namespace std {
  %template(StringPropertyMap) map<string, BrainGraph::Property>;
};
%{
  #include "../include/Node.h"
%}
%include "../include/Node.h"


// NodeList
namespace std {
  %template(NodeVector) vector<BrainGraph::Node>;
};


// Edge
%{
  #include "../include/Edge.h"
%}
%include "../include/Edge.h"


// EdgeList
%template(EdgeVector) std::vector<BrainGraph::Edge>;
%template(EdgeListVector) std::vector< std::vector<BrainGraph::Edge> >;


// Graph		
%{
   #include "../include/Graph.h"
%}
%template(LongLongVector) std::vector<long long>;
%template(GraphVector) std::vector<Graph>;
%include "../include/Graph.h"


// ROItoROi
%ignore BrainGraph::ROItoROI::begin ;
%ignore BrainGraph::ROItoROI::end ;
%{
   #include "../include/ROItoROI.h"
%}
%template(DijkstraPath) std::vector<BrainGraph::ShortestPathNode>;
%include "../include/ROItoROI.h"


// Dijkstra
%{
   #include "../include/Dijkstra.h"
%}
%template(ULongVector) std::vector<unsigned long>;
%template(IntListMap) std::map<int, std::vector<int> >;
%include "../include/Dijkstra.h"
%catches(std::out_of_range) single_source_shortest_path;
%catches(std::out_of_range) dijkstra;