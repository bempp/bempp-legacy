import fileinput, re, sys

# should be the path to core.py
core_fname = sys.argv[1]

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
reMagnitudeOrCoordinateType = re.compile(
    r"-> (?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    "("
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    ")?"
    " >::(?:MagnitudeType|CoordinateType)")
reString = re.compile(
    r"-> std::string")
reDouble = re.compile(
    r"-> double")
reUnsignedInt = re.compile(
    r"-> (?:unsigned int|size_t)")
reConstRef = re.compile(
    r"-> (.*) const &\"\"\"")
reAutoPtr = re.compile(
    r"-> std::auto_ptr< (?:Bempp|Fiber)::(\w+) >")
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

magnitudeType = {
    "float": "float",
    "double": "double",
    "std::complex< float >": "float",
    "std::complex< double >": "double"
    }

for l in fileinput.input(core_fname, inplace=1):
    orig_l = l
    l = reMagnitudeOrCoordinateType.sub("-> float", l)
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
        replacement = r"-> \1_%s" % (pythonType[type1])
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
    m = reAutoPtr.search(l)
    if m:
        replacement = r"-> \1"
        l = reAutoPtr.sub(replacement, l)
    l = l.replace("-> Bempp::", "-> ")
    l = l.replace("-> IndexSet::IndexType", "-> int")

# Remove the useless 'self' parameter from the Parameters list
f = open(core_fname, "r")
text = f.read()
f.close()

def replaceMagnitudeType(m):
    return magnitudeType[m.group(2)]

def replaceSingleType(m):
    return r"%s_%s" % (m.group(1), pythonType[m.group(2)])

def replaceDoubleType(m):
    return r"%s_%s_%s" % (m.group(1),
                          pythonType[m.group(2)],
                          pythonType[m.group(3)])

def replaceGridFunctionVector(m):
    return r"list of GridFunction_%s_%s objects" % (
        pythonType[m.group(1)],
        pythonType[m.group(2)])

def replacePreconditionerDiscreteBoundaryOperatorPtr(m):
    return r"DiscreteBoundaryOperator_%s" % (
        pythonType[m.group(1)])

def replaceParameters(m):
    text = m.group(0)
    lines = text.splitlines(True)
    assert(len(lines) >= 2)
    lines[0] = lines[0].replace("Parameters:", "*Parameters:*")
    for n in range(1, len(lines)):
        line = reParameter.sub(r"\1- \2 (\3)", lines[n])
        line = reContextQuadratureStrategy.sub(
            r"Fiber::QuadratureStrategy< \1,\2,Bempp::GeometryFactory,void >",
            line)
        line = reQuadratureStrategy.sub(r"\1 >", line)
        line = reMagnitudeType.sub(replaceMagnitudeType, line)
        line = reCoordinateType.sub(replaceMagnitudeType, line)
        line = reComplexType.sub(r"\2", line)
        line = reValueType.sub(r"\2", line)
        line = reVector.sub(r"list of \1 objects", line)
        line = rePreconditionerDiscreteBoundaryOperatorPtr.sub(
            replacePreconditionerDiscreteBoundaryOperatorPtr, line)
        line = reSingleType.sub(replaceSingleType, line)
        line = reDoubleType.sub(replaceDoubleType, line)
        line = reThyraSolveStatus.sub(
            r"Thyra::SolveStatus< \1PROTECT >", line)
        line = reAutoPtr.sub(r"\1", line)
        line = reMatrix.sub("2D array", line)
        line = reColumn.sub("1D array", line)
        line = line.replace("char const *", "string")
        line = line.replace(" const", "")
        line = line.replace(" &", "")
        line = line.replace("size_t", "int")
        line = line.replace("unsigned int ", "int")
        line = reComplexDouble.sub("complex", line)
        line = reComplexFloat.sub("complex", line)
        line = reDouble.sub("float", line)
        line = line.replace("PROTECT", "")
        line = line.replace("std::string", "string")
        line = line.replace("Bempp::", "")
        line = line.replace("Fiber::", "")
        line = line.replace("enum VtkWriter::DataType",
                            "'cell_data' or 'vertex_data'")
        line = line.replace("enum VtkWriter::OutputType",
                            "'ascii', 'base64', 'appendedraw' "
                            "or 'appendedbase64'")
        line = line.replace("enum Bempp::GridParameters::Topology",
                            "string")
        line = line.replace("enum TranspositionMode",
                            "'n', 't', 'c' or 'h'")
        lines[n] = line
    # print "".join(lines), "\n"
    return "".join(lines)

reContextQuadratureStrategy = re.compile(
    r"Bempp::Context< "
    "(float|double|std::complex< float >|std::complex< double >)"
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::QuadratureStrategy")
reQuadratureStrategy = re.compile(
    "(QuadratureStrategy< "
    "float|double|std::complex< float >|std::complex< double >"
    ","
    "float|double|std::complex< float >|std::complex< double >)"
    ",Bempp::GeometryFactory,void >")
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
    "("
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    ")?"
    " >::MagnitudeType")
reCoordinateType = re.compile(
    r"(?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    "("
    ","
    "(float|double|std::complex< float >|std::complex< double >)"
    ")?"
    " >::CoordinateType")
reThyraSolveStatus = re.compile(
    r"Thyra::SolveStatus< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >")
reComplexType = re.compile(
    r"(?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::ComplexType")
reValueType = re.compile(
    r"(?:Bempp|Fiber)::(\w+)< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::ValueType")
reVector = re.compile(
    r"std::vector< (.+),std::allocator< \1 > >")
rePreconditionerDiscreteBoundaryOperatorPtr = re.compile(
    r"Bempp::Preconditioner< "
    "(float|double|std::complex< float >|std::complex< double >)"
    " >::DiscreteBoundaryOperatorPtr")
reAutoPtr = re.compile(
    r"std::auto_ptr< (?:Bempp|Fiber)::(\w+) >")
reConstRef = re.compile(
    r" const &")
reComplexDouble = re.compile(
    r"std::complex< double >(?!PROTECT)")
reComplexFloat = re.compile(
    r"std::complex< float >(?!PROTECT)")
reDouble = re.compile(
    r"double(?!PROTECT)(?! >)")
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

