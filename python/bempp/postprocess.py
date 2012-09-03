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
reBempp = re.compile(
    r"-> Bempp::")
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

core_fname = "/home/wojtek/Projects/BEM/bempp-work/bempp/build-debug-enthought-mkl-7-recursive/python/bempp/core.py"
lib_fname = "/home/wojtek/Projects/BEM/bempp-work/bempp/build-debug-enthought-mkl-7-recursive/python/bempp/lib.py"
for l in fileinput.input((core_fname, lib_fname), inplace=1):
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
        replacement = r"-> \1_%s_%s" % (pythonType[type1], pythonType[type2])
        l = reDoubleType.sub(replacement, l)
    m = reComplicatedSpace.search(l)
    if m:
        type1 = m.group(2)
        type2 = m.group(3)
        replacement = r"-> \1_%s_%s" % (pythonType[type1], pythonType[type2])
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
    l = l.replace("-> Bempp::", "-> ")
    l = l.replace("-> IndexSet::IndexType", "-> int")
    print l,

# Remove the useless 'self' parameter from the Parameters list
f = open(core_fname, "r")
text = f.read()
f.close()

def replaceSingleType(m):
    return r"%s_%s" % (m.group(1), pythonType[m.group(2)])

def replaceDoubleType(m):
    return r"%s_%s_%s" % (m.group(1),
                          pythonType[m.group(2)],
                          pythonType[m.group(3)])

def replaceParameters(m):
    text = m.group(0)
    lines = text.splitlines(True)
    assert(len(lines) >= 2)
    lines[0] = lines[0].replace("Parameters:", ":Parameters:")
    for n in range(1, len(lines)):
        line = reParameter.sub(r"\1- \2 (\3)", lines[n])
        line = reSingleType.sub(replaceSingleType, line)
        line = reDoubleType.sub(replaceDoubleType, line)
        line = reAutoPtr.sub(r"\1", line)
        line = reMatrix.sub("2D array", line)
        line = reColumn.sub("1D array", line)
        line = reMagnitudeType.sub("float", line)
        line = line.replace("char const *", "string")
        line = line.replace(" const", "")
        line = line.replace(" &", "")
        line = line.replace("size_t", "int")
        line = line.replace("unsigned int ", "int")
        line = line.replace("std::complex< double >", "complex")
        line = line.replace("std::complex< float >", "complex")
        line = line.replace("double", "float")
        line = line.replace("std::string", "string")
        line = line.replace("Bempp::", "")
        line = line.replace("Fiber::", "")
        lines[n] = line
    return "".join(lines)

reMatrix = re.compile(
    r"arma::Mat< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >")
reColumn = re.compile(
    r"arma::Col< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >")
reSingleType = re.compile(
    r"(?:boost::shared_ptr< )?(?:Bempp|Fiber)::"
    "(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >(?: (?:const )?>)?")
reDoubleType = re.compile(
    r"(?:boost::shared_ptr< )?(?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    " >(?: (?:const )?>)?")
reMagnitudeType = re.compile(
    r"(?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::MagnitudeType")
reAutoPtr = re.compile(
    r"std:::auto_ptr< (?:Bempp|Fiber)::(\w+) >")
reConstRef = re.compile(
    r" const &")
reParameter = re.compile(
    r"^( +)(\w+): (.*)")
reParameters = re.compile(
    r"^( *)Parameters:\n(\1 +(.*)\n)+", re.MULTILINE)
reParametersWithSelfOnly = re.compile(
    r"^ *Parameters:\n *self: .*\n\n", re.MULTILINE)
reParametersBeginningWithSelf = re.compile(
    r"^( *Parameters:\n) *self: .*\n", re.MULTILINE)
text = reParametersWithSelfOnly.sub("", text)
text = reParametersBeginningWithSelf.sub("\1", text)
text = reParameters.sub(replaceParameters, text)

f = open(core_fname, "w")
f.write(text)
f.close()

