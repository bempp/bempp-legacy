// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_lib_grid_3d_geometry_decl_hpp
#define bempp_lib_grid_3d_geometry_decl_hpp

#include <dune/common/fvector.hh>

namespace Bempp {

/** Abstract wrapper of a geometry

 \todo Implement remaining methods (currently only center() is implemented,
 for the sake of an example)
 */
class Geometry {
public:
	typedef double ctype;
	enum {
		cdim = 3
	};

	/** Destructor */
	virtual ~Geometry() {
	}

	virtual Dune::FieldVector<ctype, cdim> center() const = 0;
};

/** Wrapper of a Dune geometry of type DuneGeometry */

template<typename DuneGeometry>
class ConcreteGeometry: public Geometry {
private:
	const DuneGeometry* m_dune_geometry;

	void setDuneGeometry(const DuneGeometry* dune_geometry) {
		m_dune_geometry = dune_geometry;
	}

	template<int codim, typename DuneEntity> friend class ConcreteEntity;

public:
	/** Default constructor */
	ConcreteGeometry() :
			m_dune_geometry(0) {
	}

	/** Constructor from a pointer to DuneGeometry */
	explicit ConcreteGeometry(const DuneGeometry* dune_geometry) :
			m_dune_geometry(dune_geometry) {
	}

	const DuneGeometry& duneGeometry() const {
		return *m_dune_geometry;
	}

	virtual Dune::FieldVector<ctype, cdim> center() const {
		return m_dune_geometry->center();
	}
};

} // namespace Bempp

#endif
