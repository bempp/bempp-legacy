#Copyright (C) 2011 by the BEM++ Authors
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

import os,urllib,shutil,subprocess,sys
from py_modules.tools import extract_file, to_bool


def configureAhmed(root,config):

    if not config.has_section('AHMED'):
        config.add_section('AHMED')
        config.set('AHMED','enable_ahmed','OFF')
    else:
        if not config.has_option('AHMED','enable_ahmed'): 
            config.set('AHMED','enable_ahmed','ON')
    if to_bool(config.get('AHMED','enable_ahmed')):
        config.set('AHMED','enable_ahmed','ON') # Ensure that the option has the right format
        if not (config.has_option('AHMED','lib')
                and config.has_option('AHMED','metis_lib')
                and config.has_option('AHMED','include_dir')): raise Exception('AHMED libraries or include not defined')
    else:
        config.set('AHMED','enable_ahmed','OFF')
                 

def buildAhmed(root,config):
    pass

        
    
        
        
        
            

        
