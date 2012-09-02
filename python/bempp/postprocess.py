import fileinput, re

# boost::shared_ptr< Bempp::AbstractBoundaryOperator< double,double > const >
# boost::shared_ptr< Bempp::DiscreteBoundaryOperator< std::complex< double > > const >
# boost::shared_ptr< Fiber::QuadratureStrategy< std::complex< double >,std::complex< double >,Bempp::GeometryFactory,void > >
# boost::shared_ptr< Bempp::Space< Bempp::AbstractBoundaryOperator< double,double >::BasisFunctionType > const >
# boost::shared_ptr< Bempp::Space< double > const >
# boost::shared_ptr< Bempp::Context< double,double > const >
# boost::shared_ptr< Bempp::AbstractBoundaryOperatorId const >

reSingleType = re.compile(
    r"-> boost::shared_ptr< (?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " > (?:const )?>")
reDoubleType = re.compile(
    r"-> boost::shared_ptr< (?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    "(?:,Bempp::GeometryFactory,void)?"
    " > (?:const )?>")
reComplicatedSpace = re.compile(
    r"-> boost::shared_ptr< Bempp::(Space)< Bempp::AbstractBoundaryOperator< "
    "(float|double|std::complex< float >|std::complex< double >)"
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::BasisFunctionType"
    " > (?:const )?>")
reString = re.compile(
    r"-> std::string")
reDouble = re.compile(
    r"-> double")
reUnsignedInt = re.compile(
    r"-> (?:unsigned int|size_t)")
reConstRef = re.compile(
    r"-> (.*) const &\"\"\"")
# reQuadratureStrategy = re.compile(
#     r"-> boost::shared_ptr< Bempp::(\w+)< "
#     "(float|double|std::complex< float >|std::complex< double >)"
#     ","
#     "(float|double|std::complex< float >|std::complex< double >)"
#     ",Bempp::GeometryFactory,void"
#     " > (?:const )?>")
pythonType = {
    "float": "float32",
    "double": "float64",
    "std::complex< float >": "complex64",
    "std::complex< double >": "complex128"
    }

for l in fileinput.input(("/home/wojtek/Projects/BEM/bempp-work/bempp/build-debug/python/bempp/core.py", "/home/wojtek/Projects/BEM/bempp-work/bempp/build-debug/python/bempp/lib.py"),inplace=1):
    # print l
    # if "-> boost::shared_ptr" in l:
    #     print l,
    orig_l = l
    m = reSingleType.search(l)
    if m:
        type_ = m.group(2)
        replacement = r"-> \1_%s" % pythonType[type_]
        l = reSingleType.sub(replacement, l)
    m = reDoubleType.search(l)
    if m:
        type1 = m.group(2)
        type2 = m.group(3)
        replacement = r"-> \1_%s_%s" % (pythonType[type1], pythonType[type1])
        l = reDoubleType.sub(replacement, l)
    m = reComplicatedSpace.search(l)
    if m:
        type1 = m.group(2)
        type2 = m.group(3)
        replacement = r"-> \1_%s_%s" % (pythonType[type1], pythonType[type1])
        l = reComplicatedSpace.sub(replacement, l)
    m = reString.search(l)
    if m:
        replacement = r"-> string"
        l = reString.sub(replacement, l)
    m = reDouble.search(l)
    if m:
        replacement = r"-> float"
        l = reDouble.sub(replacement, l)
    m = reUnsignedInt.search(l)
    if m:
        replacement = r"-> int"
        l = reUnsignedInt.sub(replacement, l)
    m = reConstRef.search(l)
    if m:
        replacement = "-> \\1\"\"\""
        l = reConstRef.sub(replacement, l)
    print l,
