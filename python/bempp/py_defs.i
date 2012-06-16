// Useful Python tools

%pythoncode %{
    
class Template1(object):
    def __init__(self,name,dtype1,*args,**kwargs):
        keys=globals()
        dtype1=checkType(dtype1)
        fullName=name+'_'+dtype1
        if not fullName in keys: 
            raise Exception('Class '+name+' with template type '+dtype1+' does not exist.')
        self.impl=keys[fullName](*args,**kwargs)
        self.templateType=[dtype1]

    def __getattribute__(self,attr):
        if attr=='templateType':
            return object.__getattribute__(self,'templateType')
        f=getattr(object.__getattribute__(self,'impl'),attr)
        return f



class Template2(object):
    def __init__(self,name,dtype1,dtype2,*args,**kwargs):
        keys=globals()
	dtype1=checkType(dtype1)
	dtype2=checkType(dtype2)
	fullName=name+'_'+dtype1+'_'+dtype2
        if not fullName in keys: 
            raise Exception('Class '+name+' with template types '+(dtype1,dtype2)+' does not exist.')
        self.impl=keys[fullName](*args,**kwargs)
	self.templateType=[dtype1,dtype2]
        
    def __getattribute__(self,attr):
        if attr=='templateType':
            return object.__getattribute__(self,'templateType')
        f=getattr(object.__getattribute__(self,'impl'),attr)
        return f


class Template3(object):
    def __init__(self,name,dtype1,dtype2,dtype3,*args,**kwargs):
        keys=globals()
	dtype1=checkType(dtype1)
	dtype2=checkType(dtype2)
	dtype3=checkType(dtype3)
	fullName=name+'_'+dtype1+'_'+dtype2+'_'+dtype3
        if not fullName in keys: 
            raise Exception('Class '+name+' with template type '+(dtype1,dtype2,dtype3)+' does not exist.')
        self.impl=globals()[fullName](*args,**kwargs)
	self.templateType=[dtype1,dtype2,dtype3]

    def __getattribute__(self,attr):
        if attr=='templateType':
            return object.__getattribute__(self,'templateType')
        f=getattr(object.__getattribute__(self,'impl'),attr)
        return f


def checkType(dtype):
    dtypes={'float':'float64',
            'float32':'float32',
            'float64':'float64',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

def promoteTypeToComplex(dtype):
    dtypes={'float':'complex128',
            'float32':'complex64',
            'float64':'complex128',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

  %}
