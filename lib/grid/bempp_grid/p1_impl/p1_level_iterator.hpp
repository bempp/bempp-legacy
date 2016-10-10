#ifndef bempp_p1_level_iterator_hpp
#define bempp_p1_level_iterator_hpp

namespace BemppGrid {

class P1DataContainer;    

template <int codim, typename pitype, class>
class P1LevelIteratorImp {

    public:

        P1LevelIteratorImp(const shared_ptr<const P1DataContainer>& data, int level) :
            m_data(data), m_level(level), m_index(0) {}

    private:

        shared_ptr<const P1DataContainer> m_data;
        int m_level;
        int m_index;




};

}




#endif
