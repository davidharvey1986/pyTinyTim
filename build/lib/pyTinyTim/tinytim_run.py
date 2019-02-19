import os as os
import subprocess as process
import ipdb as pdb

def run(x, y, **args):

    '''
    Wrapper for the IDL function tiny tim as created by
    Richard

    INPUTS : X,Y : A VECTOR CONTAINING THE PIXEL POSITIONS
             OF THE STARS WE WANT

    ARGS : ALL THE ARGS WE WANT TO GO INTO TINYTIM_CREATE
           IF IT IS A STRING IT MUST HAVE QUOTES AROUND IT

    EXAMPLE :
    tinytim.run( x, y, focus_range=focus[iImage], \
                            pixel_scale=0.03, \
                            output_dir="'"+outputDir+"'", \
                            filter=filter_string,\
                            raw=1, fitsname="'"+images[iImage]+"_TT'",\
                            exact_position=1)

    REQUIRES : TINYTIM
               TINYTIM_CREATE.PRO (IDL)

    '''
    #Set up the the string
    args_string = [ str(iArg)+'='+str(args[iArg]) for iArg in args ]
    code_dir='/Users/DavidHarvey/Library/Code/IDL/tinytim/IDL'
    #write coordfile
    TT_file = open(code_dir+"/coordfile.TT","wb")
    try:
        len(x)
    except:
        x=[x]
        y=[y]
    for i in xrange(len(x)):
        TT_file.write( str(x[i])+"   "+str(y[i])+"\n")
    TT_file.close()
    args_string.append("coordfile='"+code_dir+"/coordfile.TT'")
    args_list = ','.join(args_string)
    command_str ='idl -e "tinytim_create,'+args_list+'"'
    os.system(command_str)
    print 'DONE'
