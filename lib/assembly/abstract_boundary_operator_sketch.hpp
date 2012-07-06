#ifndef bempp_abstract_boundary_operator_hpp
#define bempp_abstract_boundary_operator_hpp

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class Context
{
public:
    Context(const AssemblyOptions& assemblyOptions,
            const shared_ptr<LocalAssemblerFactory<BasisFunctionType, ResultType> >&
            localAssemblerFactory);

    static shared_ptr<const Context<BasisFunctionType, ResultType> >& defaultContext();
    static void setDefaultContext(
            const shared_ptr<const Context<BasisFunctionType, ResultType> >&
            newContext);

    shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
    getWeakForm(const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const {
        return m_cache.getWeakForm(*this, op);
    }

private:
    static shared_ptr<const Context<BasisFunctionType, ResultType> > m_defaultContext;

    AssemblyOptions m_assemblyOptions;
    shared_ptr<LocalAssemblerFactory<BasisFunctionType, ResultType> > m_localAssemblerFactory;

    mutable WeakFormCache<BasisFunctionType, ResultType> m_cache;
};

// TODO: deal with template params (this should be instantiated for specific
// types, like double-double, double-complex<double>...
shared_ptr<Context<BasisFunctionType, ResultType> > Context::m_defaultContext(
        make_shared<DefaultLocalAssemblerFactoryForSurfaceOperators<BasisFunctionType, ResultType> >());

template <typename BasisFunctionType, typename ResultType>
class WeakFormCache
{
public:
    shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
    getWeakForm(const Context<BasisFunctionType, ResultType>& context,
                const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const {
        shared_ptr<const Key> key = op.key();
        if (key) // operator is cacheable
            ; // look into cache, if not there, call op.assembleWeakForm(context), try inserting result into map
        else
            return op.assembleWeakForm(context);
    }
};

class Key
{
public:
    virtual ~Key() {}
    virtual bool less(const Key& other) const = 0;
};

template <typename BasisFunctionType, typename ResultType>
class AbstractBoundaryOperator
{
public:
    // most methods from current BoundaryOperator: domain() ...
    /** \brief Domain.
     *
     *  Return a reference to the function space being the domain of the operator. */
    const Space<BasisFunctionType>& domain() const;

    /** \brief Range.
     *
     *  Return a reference to the function space being the range of the operator. */
    const Space<BasisFunctionType>& range() const;

    /** \brief Dual to range.
     *
     *  Return a reference to the function space dual to the range of the operator. */
    const Space<BasisFunctionType>& dualToRange() const;

    /** @}
     *  @name Label
     *  @{ */

    /** \brief Return the label of the operator. */
    std::string label() const;

    /** @}
     *  @name Assembly
     *  @{ */

    /** \brief Return true if this operator supports the representation \p repr. */
    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const = 0;

    virtual shared_ptr<DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
    assembleWeakForm(const Context& context) = 0;

    // default implementation returns an invalid key. Override if you
    // want the weak form to be cached.
    virtual shared_ptr<const Key> key() { return shared_ptr<const Key>(); }

    bool isCacheable() const { return key().get() != 0; }
};

//template <typename BasisFunctionType, typename ResultType>
//class CacheableBoundaryOperator // replaces ElementaryBO
//{
//public:
//    virtual shared_ptr<DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
//    weakForm(const Context& context);

//    virtual std::auto_ptr<DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
//    assembleWeakForm(const Context& context) = 0;

//};

//template <typename BasisFunctionType, typename ResultType>
//shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
//CacheableBoundaryOperator::weakForm(const Context& context) const
//{
//    return context->getWeakForm(*this); // this may call assembleWeakForm(context)
//}

//template <typename BasisFunctionType, typename ResultType>
//shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
//ScaledBoundaryOperator<BasisFunctionType, ResultType>::weakForm(const Context& context) const
//{
//    return assembleWeakForm(context); // don't try to look into cache, it's a lightweight operator
//}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
ScaledBoundaryOperator<BasisFunctionType, ResultType>::assembleWeakForm(const Context& context) const
{
    shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> > wf =
            m_op.weakForm(context);
    // create a ScaledDiscreteBO wrapping wf and return it
}

template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator
{
public:
    BoundaryOperator(const shared_ptr<const AbstractBoundaryOperator<
                     BasisFunctionType, ResultType> >& abstractOp,
                     const shared_ptr<const Context<
                     BasisFunctionType, ResultType> >& context =
                     Context<BasisFunctionType, ResultType>::defaultContext());

    shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> >
    abstractOperator() const;
    shared_ptr<const Context> context() const;

    shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
    weakForm() const;

private:
    shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> >
    m_abstractOp;
    shared_ptr<const Context> m_context;
};

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<BasisFunctionType, ResultType> >
BoundaryOperator::weakForm() const
{
    return m_context->getWeakForm(m_abstractOp);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSingleLayerBoundaryOperator(
        const shared_ptr<const Space>& domain,
        const shared_ptr<const Space>& range,
        const shared_ptr<const Space>& dualToRange,
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context =
        Context<BasisFunctionType, ResultType>::defaultContext())
{
    return BoundaryOperator<BasisFunctionType, ResultType>(
                make_shared<Laplace3dSingleLayerBoundaryOperator<
                BasisFunctionType, ResultType> >(domain, range, dualToRange),
                context);
}

template <typename BasisFunctionType>
class Laplace3dKey : public Key
{
public:
    Laplace3dKey(const std::type_info* operatorType,
                 const Space<BasisFunctionType>* domain,
                 const Space<BasisFunctionType>* range,
                 const Space<BasisFunctionType>* dualToRange);

    virtual bool less(const Key& other) const {
        const Laplace3dKey* otherLaplaceKey =
                dynamic_cast<const Laplace3dKey*>(&other);
        if (otherLaplaceKey == 0)
            return typeid(*this).before(typeid(other));
        else
            return (m_operatorType->before(*otherLaplaceKey->m_operatorType) ||
                    (*m_operatorType == *otherLaplaceKey->m_operatorType &&
                     (m_domain < otherLaplaceKey->m_domain ||
                      (m_domain == otherLaplaceKey->m_domain &&
                       (m_range < otherLaplaceKey->m_range ||
                        (m_range == otherLaplaceKey->m_range &&
                         m_dualToRange < otherLaplaceKey->m_dualToRange))));
    }

private:
    const std::type_info* m_operatorType;
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

} // namespace Bempp

#endif // bempp_abstract_boundary_operator_hpp
