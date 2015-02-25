#ifndef LAGRANGE_INTERFACE_HPP
#define LAGRANGE_INTERFACE_HPP

#include <Python.h>
#include "bempp/space/py_space_variants.hpp"
#include "bempp/common/shared_ptr.hpp"
#include "bempp/grid/grid.hpp"


namespace Bempp {

   static inline PyObject* _py_p1_vertex_map(
           const SpaceVariants& spaceVariant){

       const Space<double>& space = *(_py_get_space_ptr<double>(spaceVariant));

       auto grid = space.grid();
       auto view = grid->leafView();
       const IndexSet& indexSet = view->indexSet();

       npy_intp globalDofCount = space.globalDofCount();


       PyObject* vertexMap;
       vertexMap = PyArray_SimpleNew(1,&globalDofCount,NPY_INT);
       for (npy_intp i = 0; i < globalDofCount; ++i)
           *((int*)PyArray_GETPTR1(vertexMap,i)) = -1;
           
       std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
       while(!it->finished()){
           const Entity<0> &element = it->entity();
           std::vector<GlobalDofIndex> dofs;
           space.getGlobalDofs(element,dofs);
           for (int i = 0; i < 3;++i){
               int index = indexSet.subEntityIndex(element,i,2);
               *((int*)PyArray_GETPTR1(vertexMap,index)) = dofs[i];
           }
           it->next();
       }
       return vertexMap;


       } 
}
#endif
