template <typename CoordinateType>
class QuadratureDescriptorSelectorFactory
{
public:
    virtual ~QuadratureDescriptorSelectorFactory() {}

    virtual shared_ptr<QuadratureDescriptorSelectorForIntegralOperators<CoordinateType> >
    makeQuadratureDescriptorSelectorForIntegralOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<int> >& testShapesetOrders,
        const shared_ptr<const std::vector<int> >& trialShapesetOrders) const = 0;

    virtual shared_ptr<QuadratureDescriptorSelectorForGridFunctions<CoordinateType> >
    makeQuadratureDescriptorSelectorForGridFunctions(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<int> >& testShapesetOrders) const = 0;

    virtual shared_ptr<QuadratureDescriptorSelectorForLocalOperators<CoordinateType> >
    makeQuadratureDescriptorSelectorForLocalOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<int> >& testShapesetOrders,
        const shared_ptr<const std::vector<int> >& trialShapesetOrders) const = 0;
};

template <typename CoordinateType>
class DefaultQuadratureDescriptorSelectorFactory
{
};

template <typename CoordinateType>
shared_ptr<QuadratureDescriptorSelectorForIntegralOperators<CoordinateType> >
DefaultQuadratureDescriptorSelectorFactory<CoordinateType>::
makeQuadratureDescriptorSelectorForIntegralOperators(
    const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
    const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
    const shared_ptr<const std::vector<int> >& testShapesetOrders,
    const shared_ptr<const std::vector<int> >& trialShapesetOrders) const
{
    typedef DefaultQuadratureDescriptorSelectorForIntegralOperators<CoordinateType> Selector;
    return boost::make_shared<Selector>(testRawGeometry, trialRawGeometry,
                                        testShapesetOrders, trialShapesetOrders);
}

template <typename CoordinateType>
class QuadratureDescriptorSelectorForIntegralOperators
{
public:
    virtual ~QuadratureDescriptorSelectorForIntegralOperators() {}

    virtual DoubleQuadratureDescriptor quadratureDescriptor(
        int testElementIndex, int trialElementIndex,
        CoordinateType nominalDistance) const = 0;
};

template <typename CoordinateType>
class DefaultQuadratureDescriptorSelectorForIntegralOperators :
    public QuadratureDescriptorSelectorForIntegralOperators<CoordinateType>
{
public:
    DefaultQuadratureDescriptorSelectorForIntegralOperators(
        const AccuracyOptionsEx& accuracyOptions,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& trialShapesets) :
        m_testRawGeometry(testRawGeometry),
        m_trialRawGeometry(trialRawGeometry),
        m_testShapesets(testShapesets),
        m_trialShapesets(trialShapesets)
    {
        Utilities::checkConsistencyOfGeometryAndShapesets(
            *testRawGeometry, *testShapesets);
        Utilities::checkConsistencyOfGeometryAndShapesets(
            *trialRawGeometry, *trialShapesets);

        precalculateElementSizesAndCenters();
    }

    virtual DoubleQuadratureDescriptor quadratureDescriptor(
        int testElementIndex, int trialElementIndex,
        CoordinateType nominalDistance) const = 0;
};

template <typename CoordinateType>
class QuadratureDescriptorSelectorForLocalOperators
{
public:
    virtual ~QuadratureDescriptorSelectorForLocalOperators() {}

    virtual SingleQuadratureDescriptor quadratureDescriptor(
        int elementIndex) const = 0;
};

template <typename CoordinateType>
class QuadratureDescriptorSelectorForGridFunctions
{
public:
    virtual ~QuadratureDescriptorSelectorForGridFunctions() {}

    virtual SingleQuadratureDescriptor quadratureDescriptor(
        int elementIndex) const = 0;
};

template <typename CoordinateType>
class DoubleQuadratureRuleFamily
{
public:
    virtual ~DoubleQuadratureRuleFamily() {}

    virtual void fillQuadraturePointsAndWeights(
        const DoubleQuadratureDescriptor& desc,
        arma::Mat<CoordinateType>& testPoints,
        arma::Mat<CoordinateType>& trialPoints,
        std::vector<CoordinateType>& testWeights,
        std::vector<CoordinateType>& trialWeights,
        bool& isTensor) const = 0;
};

template <typename CoordinateType>
class SingleQuadratureRuleFamily
{
public:
    virtual ~SingleQuadratureRuleFamily() {}

    virtual void fillQuadraturePointsAndWeights(
        const SingleQuadratureDescriptor& desc,
        arma::Mat<CoordinateType>& points,
        std::vector<CoordinateType>& weights) const = 0;
};

