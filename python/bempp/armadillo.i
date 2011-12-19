// Generate a minimal wrapper for arma::Mat<double>
// (used by the wrapper of Bempp::VtkWriter, not intended to be accessed by the user)
namespace arma {
    template<typename eT>
    class Mat {
    public:
        ~Mat();
    };

    %template(__Mat_double) Mat<double>;
}

// Define typemaps
%include "armadillo2numpy.i"





