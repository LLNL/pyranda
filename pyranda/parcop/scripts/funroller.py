import copy
import sys,os
import re

try:
    file_name = sys.argv[1]
except:
    print("Error: Please specify filename")
    exit()

try:
    fsuff = sys.argv[2]
except:
    fsuff = ".f"
    #print("Using %s as fortran file suffix" % fsuff)

try:
    profile = bool(int(sys.argv[3]))
except:
    profile = True
    #print("Adding profiling tags")



# Asynch mode currently does not work.  A known bug has been reported to IBM
#  that causes asynchronous code to fail when spanning subroutine calls and map statements
#  When bug is fixed, we can nominally switch both these flags to True.
asyncMode = False  # Only switches FEXL statments that have been flagged
asyncAll  = False  # Switches all statements to async 

sync_var = "sync_var"

verbose = False #True
replace = True #False


# Some markup keywords
start_unroll = "!$FEXL"
end_unroll = "!$END FEXL"
declare_unroll = "!$DEF-FEXL"

# fexl-async macro
fexl_async = "!$fexl-async"

# Parameters
MAXLINES = 1e6

# Wild card symbol.. trailing only
wcs = "@"

# Indices used in new loops
INDICES = ['ifunr','jfunr','kfunr','nfunr']
INDrank = [1      ,2      ,3      ,4      ]

# OMP pragmas
ompStart = """
!$omp target teams distribute parallel do collapse(%s)"""

ompEnd = """!$omp end target teams distribute parallel do

"""
# Map - target enter|exit data
ompMap = "!$omp target %s "

## Async option
ompAsync = " nowait depend(inout:%s) "

newline = "\n"


# Caliper profiling
caliperStart = """CALL exosim_annotation_begin("%s-L.%s")\n"""
caliperEnd   = """CALL exosim_annotation_end("%s-L.%s")\n"""


# Open file and get lines
pid = open(file_name,'r')
lines = pid.readlines()
pid.close()

has_any = False
for line in lines:
    if start_unroll in line:
        has_any = True
        break
    
# Also check for async macro
for line in lines:
    if ("!$omp " in line) and (fexl_async in line):
        has_any = True
    
    
if not(has_any):
    #print('funroller - No loops found in %s' % file_name )
    exit()

# For replace... make a backup
if replace:
    os.system('mv %s %s' % (file_name,file_name.replace(fsuff,fsuff + ".fexl") ) )


def matchPar(text):
    istart = []  # stack of indices of opening parentheses
    d = []
    for i, c in enumerate(text):
        if c == '(':
             istart.append(i)
        if c == ')':
            try:
                #d[istart.pop()] = i
                d.append( [istart.pop(), i ] )
            except IndexError:
                print('Too many closing parentheses')
    if istart:  # check if stack is empty afterwards
        print('Too many opening parentheses')

    return(d)


class iFunc():
    def __init__(self,start,end,kind,filename=""):
        self.start = start      # Star line number
        self.end   = end        # End line number
        self.header = None    # Header line number
        self.kind  = kind
        self.filename = filename
        self.funcBody = []
        self.funcName = ""
        self.loops = None
        self.returns = []
        self.defcode = []
        self.endcode = []

        self.kind = ""
        if 'sub' in kind:
            self.kind = "subroutine"
        else:
            self.kind = "function"


    def getFuncBody(self,lines):
        """
        Given a loop, extract the body
        """
        self.funcBody = lines[self.start:self.end+1]


    def getHeader(self):
        """
        Find the line number of the end of the DEF-FEXL area
        """
        for ii in range(len(self.funcBody)):
            fline = self.funcBody[ii].strip().lower()
            if declare_unroll.lower() in fline:
                self.header = ii + self.start
        
        return
        

    def getFuncName(self):
        
        self.funcName = self.loopBody[0].lower().strip().split(self.kind)[0]


    def getDEFcode(self):        
        defcode = []

        # Only add def indices if there are loops
        if ( len(self.loops) > 0 and self.header):
            defcode += ['integer :: %s,%s,%s,%s \n' % tuple(INDICES)]

        # Add loop add ons
        for loop in self.loops:
            if loop.AddToDEFFEXL:
                defcode += loop.AddToDEFFEXL

        self.defcode = defcode
    
    def sortLoops(self):
        """ Sort loops by start line number """
        # In-place sort
        self.loops.sort(key=lambda x: x.start, reverse=False)


    def getENDcode(self):
        endcode = []
        for loop in self.loops:
            if loop.AddToEndFunction:
                endcode += loop.AddToEndFunction
        
        self.endcode = endcode

    def makeFrontMatter(self):
        """ Make the code and add it to the loops for function front matter
        """
        self.getDEFcode()
        if self.defcode:
            dloop = iLoop(self.header,self.header )  # Make it a loop object
            dloop.newBody = self.defcode             # Add defcode as body
            parms = self.loops[0].parms              # Check parms
            self.loops.append( dloop )               # Add it to list
            
    def makeMiddleMatter(self):


        # For functions with memory maps, free before "RETURN" statements
        mem_parms = ['memory','vars_']
        has_mem = False
        for ll in self.loops:
            if ll.parms:
                if any( ll.parms.has_key(mm) for mm in mem_parms):
                    has_mem = True
        
        if has_mem:
            for ii in range(len(self.funcBody)):
                fline = " " + self.funcBody[ii].strip().lower() + " "
                if " return " in fline:
                    dloop = iLoop(ii+self.start,ii+self.start)
                    dloop.newBody = self.endcode + ["\n"] + ["RETURN \n"]
                    self.loops.append(dloop)


    def makeEndMatter(self):
        self.getENDcode()
        if self.endcode:
            dloop = iLoop(self.end,self.end)
            #endf = lines[self.end+offset]
            endf = lines[self.end]
            dloop.newBody = self.endcode + ["\n"] + [endf]
            self.loops.append(dloop)



class iLoop():

    def __init__(self,start,end,filename=''):
        self.start  = start   # Starting line of loop
        self.end    = end     # Ending line of the loop
        self.loopBody = []
        self.parent = None

        self.parms = None
        self.newBody = []
        self.indices = INDICES
        self.ranks = INDrank
        self.dim = 'dim'
        self.var = 'var'
        self.fix = 'fix'
        self.sumvar = 'sumvar'
        self.bounds = 'bounds'
        self.replace = 'replace'
        self.private = 'private'
        self.shared  = 'shared'
        self.async   = 'async'
        self.vars_to = "vars_to"
        self.vars_alloc = "vars_alloc"
        self.vars_alloc_end = "vars_alloc_end"
        self.vars_from = "vars_from"
        self.vars_delete = "vars_delete"
        self.vars_update_to = "vars_update_to"
        self.vars_update_from = "vars_update_from"
        self.pool = "pool"
        self.return_early = "return_early"
        self.memory = "memory"
        self.filename = filename
        self.AddToEndFunction = []
        self.AddToDEFFEXL = []
        
    def getLoopBody(self,lines):
        """
        Given a loop, extract the body
        """
        self.loopBody = lines[self.start:self.end+1]

    def showConvert(self):

        print("Original loop body")
        for code in self.loopBody:
            print(code.replace('\n',''))
        
        print("Converted loop body")
        for code in self.newBody:
            print(code.replace('\n',''))

            
    def convert(self,lines,offset):

        n1 = len(lines)
        # Delete old code
        Tlines = self.end - self.start
        istart = self.start + offset
        for ii in range(Tlines+1):
            del lines[istart]
            
        # Add new code
        for code in self.newBody[-1::-1]:
            lines.insert( istart,code)

        n2 = len(lines)
        return n2-n1
            

            
            
    def convertLoop(self):

        new_code = []
        
        # Get unroll args
        dim    = self.dim
        var    = self.var
        sumvar = self.sumvar
        fix    = self.fix
        bounds = self.bounds
        replace= self.replace
        private= self.private
        shared = self.shared
        async  = self.async
        vars_to     = self.vars_to
        vars_alloc  = self.vars_alloc
        vars_alloc_end  = self.vars_alloc_end
        vars_from   = self.vars_from
        vars_delete = self.vars_delete 
        vars_update_to = self.vars_update_to
        vars_update_from = self.vars_update_from
        pool = self.pool
        return_early = self.return_early
        memory = self.memory

        
        try:
            unroll = self.loopBody[0]
            #parms = eval(unroll.strip().replace(start_unroll,'').strip())

            # Get FEXL argument list (can be multi-lined)
            unroll = ''
            for ii in range(len(self.loopBody)):
                if self.loopBody[ii].strip().startswith( start_unroll ):
                    unroll += self.loopBody[ii].strip().replace(start_unroll,'').strip()
                else:
                    break
            loopStart = ii
            
            parms = eval(unroll.strip().replace(start_unroll,'').strip())                    
            
            comment = "" #unroll.split('!$')[0]
            indent  = ""

        except:
            print("Error in loop (l.%s) : %s" % (self.start+1,unroll) )
            import pdb
            pdb.set_trace()

        # Check parms and case out for empty
        self.parms = parms
        if not(self.parms):
            return


        # Make data/target asynchonous
        extra = ""
        if ( (async in parms) and asyncMode):
            extra = ompAsync % parms[async] # Data

        # For asyncAll, override all omp
        if asyncAll:
            extra = ompAsync % sync_var
            
            
        omp_map   = ompMap   + extra  #+ newline # Data
        # nuclear option
        #omp_map   ="""
        #!$omp barrier
        #!$omp taskwait
        #%s
        #!$omp taskwait
        #!$omp barrier""" % omp_map        
        
        omp_start = ompStart + extra #+ newline # Compute
            

        #####################################            
        # FEXL-DATA + OMP
        # Look for any data kind
        if ( (vars_to     in parms) or
             (vars_from   in parms) or
             (vars_alloc  in parms) or
             (vars_alloc_end  in parms) or
             (vars_delete in parms) or
             (vars_update_from in parms) or
             (vars_update_to in parms)  or 
             (memory in parms) ):
             
            # Check for enters/exit/update
            dataEnter  = (vars_to in parms)    or (vars_alloc in parms)
            dataExit   = (vars_from in parms)  or (vars_delete in parms)
            dataUpdate = (vars_update_to in parms) or (vars_update_from in parms)

            # Check for pool request
            makePoolAlloc   = (pool in parms) and (vars_alloc     in parms)
            makePoolAllocEnd   = (pool in parms) and (vars_alloc_end  in parms)
            makePoolRelease = (pool in parms) and (vars_delete    in parms)
            endAssociate    = True

            # Check for inline defined temporarly memory, dynamic
            dynamicMemory = (memory in parms)


            if dynamicMemory:
                varFlag = "!$VAR-FEXL"
                dimMatch = "dimension("
                alloc_type = parms[memory]
                for iline in self.loopBody[1:-1]:
                    if varFlag in iline:
                        vlist = iline.split("::")[1].split("!")[0].strip().split(',')  # Variables in this line
                        # FEXL POOL VERSION
                        if alloc_type == 'pool':
                            cmd = ""
                            # Get dimension info
                            dimstr = iline[iline.lower().find(dimMatch):]
                            pars = matchPar(dimstr.split("!")[0])
                            dim_args = dimstr[ :pars[-1][1] + 1]
                            size_args = dimstr[ pars[-1][0]+1:pars[-1][1] ]
                            size_ = re.sub( "\(.*?\)", " ", size_args )
                            dim = len( size_.split(',') )
                            code = iline.replace(dim_args,dim_args + " ,POINTER ")
                            code = code.replace(size_args,":" + ",:"*(dim-1))
                            code = code.replace(varFlag,"!")
                            new_code.append( code )
                            myPool = iline.split(varFlag)[1].strip()
                            mkpool = []
                            for ivar in vlist:
                                mkpool.append("%s=>fexlPool%s_get() \n" % (ivar,myPool) )
                            endPool = "pool_error = fexlPool%s_free(%s)\n" % (myPool,len(vlist)) 
                            self.AddToDEFFEXL +=  mkpool            # Order of appearance
                            self.AddToEndFunction.insert(0,endPool) # Reverse order
                            #import pdb
                            #pdb.set_trace()
                        # OMP VERSION
                        elif alloc_type == 'omp':
                            cmd = iline.replace(varFlag,'!')
                            new_code.append( cmd )
                            enterData = ompMap % "enter data map(alloc: %s )\n"  % ",".join(vlist)
                            exitData  = ompMap % "exit  data map(release: %s)\n" % ",".join(vlist)
                            self.AddToDEFFEXL.append(enterData)      # Order of appearance
                            self.AddToEndFunction.insert(0,exitData) # Reverse order
                            
                    else:
                        new_code.append( iline )
                


            if return_early in parms:
                endAssociate = not parms[return_early]
                                    
            # If pools, then turn off dataEnter/Exit
            if makePoolAlloc:
                dataEnter = False
            if makePoolRelease:
                dataExit = False

            if makePoolRelease:
                myPool  = parms[pool]
                omp_cmd = ""
                if endAssociate:
                    omp_cmd += "END ASSOCIATE" + newline
                vlist = parms[vars_delete].split(",")
                omp_cmd += "fexlpool_error = fexlPool%s_free(%s)" % (myPool,len(vlist))
                omp_cmd +=  newline
                new_code.append( omp_cmd )

            if makePoolAllocEnd:
                myPool  = parms[pool]
                omp_cmd = ""
                if endAssociate:
                    omp_cmd += "END ASSOCIATE" + newline
                new_code.append( omp_cmd )

            if makePoolAlloc:
                myPool  = parms[pool]
                omp_cmd = "ASSOCIATE("
                vlist = parms[vars_alloc].split(",")
                omp_cmd += "%s=>fexlPool%s_get()" % (vlist[0],myPool)
                for ivar in vlist[1:]:
                    omp_cmd += ",%s=>fexlPool%s_get()" % (ivar,myPool)
                omp_cmd += ")" + newline
                new_code.append( omp_cmd )
            
            
            
            # Supports all 3, but goes in order.... exit, enter, update
            if dataExit:
                omp_cmd = omp_map % "exit data %s" 

                mapp = ""
                # From 
                if ( vars_from in parms ):
                    mapp += " map(from:%s) " % parms[vars_from]
                
                # Delete->release
                if (vars_delete in parms):
                    mapp += " map(release:%s) " % parms[vars_delete]

                omp_cmd = omp_cmd % mapp
                omp_cmd += newline

                # Drop nowait and add barrier after the exit
                if asyncMode:
                    omp_cmd = omp_cmd.replace('nowait','')
                    omp_cmd += "!$omp barrier" + newline
                
                new_code.append( omp_cmd )

            if dataEnter:
                omp_cmd = omp_map % "enter data %s" 

                mapp = ""
                # To
                if ( vars_to in parms ):
                    mapp += " map(to:%s) " % parms[vars_to]
                
                # Alloc
                if (vars_alloc in parms):
                    mapp += " map(alloc:%s) " % parms[vars_alloc]

                omp_cmd = omp_cmd % mapp
                omp_cmd += newline
                new_code.append( omp_cmd )

            if dataUpdate:

                # update to
                omp_cmd = omp_map % "update %s"
                if ( vars_update_to in parms ):
                    omp_cmd = omp_cmd % ( " to(%s) " % parms[vars_update_to] )
                    omp_cmd += newline
                    new_code.append( omp_cmd )            
                
                # update from
                omp_cmd = omp_map % "update %"
                if (vars_update_from in parms):
                    omp_cmd = omp_cmd % (" from(%s) " % parms[vars_update_from] )
                    omp_cmd += newline
                    new_code.append( omp_cmd )
                    
            self.newBody = new_code
            return None
        
        #####################################
        # FEXL-LOOPS + OMP
        skipdo = False

        if not parms.has_key(self.dim):
            parms[self.dim] = 0

        collapse = parms[self.dim]
        if ( fix in parms ):
            collapse -= 1
        if collapse <= 0:
            skipdo = True

        # Caliper profiling
        if profile and not skipdo:
            new_code.append(caliperStart % (self.filename,self.start))

        # OMP loops start and other options
        if collapse > 0:
            new_code.append(omp_start % collapse)

        # OMP private variable list
        if ( private in parms):            
            s_priv = parms[private][0]
            for pp in parms[private][1:]:
                s_priv += ", %s" % pp
            new_code.append(" & " + newline)
            new_code.append("!$omp   private(%s)" % s_priv )

        # OMP shared variable list
        if ( shared in parms):            
            s_shar = parms[shared][0]
            for pp in parms[shared][1:]:
                s_shar += ", %s" % pp
            new_code.append(" & " + newline)
            new_code.append("!$omp   shared(%s)" % s_shar )

        # OMP newline
        new_code.append( newline )
                

        # Check for var = [] cases... omp only
        if (not var in parms ):
            skipdo = True
            parms.update({'var':[]})

        
        # Skip the collapse/Do loops for certain cases
        if not skipdo:

            # Unrolled header
            my_indices = self.indices[:parms[self.dim]]
            my_ranks   = self.ranks[  :parms[self.dim]]

            # Check for fixed indices (fortran index, 1 is first)
            if ( fix in parms ):
                ifix = parms[fix] - 1  # Convert fortran to python starting index
                my_indices.pop(ifix)
                my_ranks.pop(ifix)

            my_indices.reverse()
            my_ranks.reverse()

            # Check for explcit loop bounds (all or nothing)
            expBounds = False
            if ( bounds in parms ):
                ibnd = parms[ bounds ]  # bounds: "imin:imax,jmin:jmax,kmin:kmax"
                expBounds = True


            for ind,rnk in zip(my_indices,my_ranks):

                # Automatic bnds here (default)
                size_var = ""
                for vv in parms[self.var]:
                    if wcs not in vv:
                        size_var =  vv  # Use this var for size, first one that doesnt have a wild card
                        break
                if not size_var:
                    print("Error: File: %s, line: %s -- Must give a non-wild card variable in vars list" %(self.filename,self.start) )
                    exit()

                bndStr = "do %s=%s,size(%s,%s)\n" %  (ind,1, size_var , rnk )

                # Explicit range option here
                if expBounds:
                    min1 = ibnd.split(',')[rnk-1].split(':')[0]
                    max1 = ibnd.split(',')[rnk-1].split(':')[1]
                    bndStr = "do %s=%s,%s \n" % (ind,min1,max1)
                                
                code = (comment + indent + bndStr)
                            
                new_code.append( code )
                indent += '  '

            
            # Form variable "index" which is thing appended to the variable "var(i,j,k)" -> index = '(i,j,k)'
            if ( not fix in parms ):
                my_indices.reverse()
                index = '(%s' % my_indices[0]
                for dd in range(1,parms[self.dim]):
                    index += ',' + my_indices[dd]
                index += ')'


                # Colon syntax
                colon = "(:"
                for dd in range(1,parms[self.dim]):
                    colon += ",:"

            else:

                my_indices.reverse()
                my_ranks.reverse()            

                colonL = ""
                indexL = ""
                colonR = ""
                indexR = ""
                cnt = 1
                for col,ind in zip(my_ranks,my_indices):
                    if col == 1:
                        colonL = "(:"
                        indexL = '(' + ind
                    elif col == parms[self.dim]:
                        colonR += ",:)"
                        indexR += ',' + ind + ")"
                    elif col == cnt:
                        colonL +=",:"
                        indexL += ',' + ind
                    else :
                        colonR +=",:"
                        indexR += ',' + ind
                    cnt += 1

            # Protect
            protect = {}
            protect["PRESENT("] = "#present#"
            protect["SIZE("]    = "#size#"

            leftColon  = "#LEFTCOLON#"
            rightColon = "#RIGHTCOLON#"
            sColon = "#COLON#"
        
        for code in self.loopBody[loopStart:-1]:

            esc = ['\n','+','-','*','/',';','=','MAX','MIN','max','min']
            esc += ['(',')',',']


            # Do global find and replace
            if ( replace in parms):
                replaces = parms[replace]
                for rep in replaces:
                    code = code.replace(rep,replaces[rep])

            
            if not skipdo:
                 

                # Find colon syntax

                # Left and right colon for fixed index cases
                if ( fix in parms ):
                    if colonL:
                        code = code.replace(colonL, '%s ' % leftColon )
                    if colonR:
                        code = code.replace(colonR,'%s ' % rightColon )
                else:
                    # Find colon syntax
                    code = code.replace(colon,'%s ' % sColon )

                # SUM reductions case
                code = code.replace("DIM=4","DIM=1")

                # Where statements
                if ( "WHERE" in code ):
                    code = code.replace("ELSEWHERE","ELSE")
                    code = code.replace("WHERE"," IF")
                    if ("END" not in code) and ("ELSE" not in code):
                        #import pdb
                        #pdb.set_trace()
                        code = code.replace("\n"," THEN \n")

                
                # Protected strings
                for p in protect:
                    if ( p in code):
                        code = code.replace( p , protect[p] )

                # White space separate esc string for correct parsing
                for e in esc:
                    code = code.replace(e," %s "%e)   
                        
                # Add white space around all chars/operators 
                mycode =  [e+" " for e in code.split(" ") if e]
                
                
                icode = []
                # Loop over items in this line of code
                for myvar in mycode:

                    if ( not fix in parms ):
                        
                        # Check for 1:1 match
                        tvar = myvar.strip()

                        # For a var with a wildcard (wcs)
                        has_wildcard = [1 if wcs in vv else 0 for vv in parms[self.var]]
                        wc_matches = False
                        if any(has_wildcard):
                            dd = []
                            for vv in parms[self.var]:
                                vstr = vv.strip()
                                if wcs in vstr[0]:  # First char wild.. match 'ends with'
                                    idd = tvar.endswith( vstr.replace(wcs,'') )
                                elif wcs in vstr[-1]:  # Last char wild... match 'startswith'
                                    idd = tvar.startswith( vstr.replace(wcs,'') )
                                else:
                                    idd = False
                                dd.append( idd )         
                            wc_matches = any( dd )

                        # Check for colon operators
                        if sColon in myvar.strip():             
                            code = myvar.strip().replace(sColon,index[:-1])
                            icode.append(code)
                            
                        # Simple check for direct match or wc_match
                        elif (tvar in parms[self.var]) or wc_matches:
                            code =  myvar.strip() + index
                            icode.append(code)


                        # Special case for sum on 4th index
                        elif 'sumvar' in parms:                        
                            if myvar.strip() in parms[self.sumvar]:
                                code = myvar.strip() + index.replace(')',',:)')
                                icode.append(code)
                            else:
                                icode.append(myvar)

                        # Everything else
                        else:
                            icode.append(myvar)

                    else:
                        
                        # Check for colon operators for fixed index cases
                        if leftColon in myvar.strip():
                            
                            code = myvar.strip().replace(leftColon,indexL)
                            icode.append(code)

                        elif rightColon in myvar.strip():             
                            code = myvar.strip().replace(rightColon,indexR)
                            icode.append(code)                                                        

                        # Everything else
                        else:
                            icode.append(myvar)
                            
                ## End "mycode" loop
                code = indent + comment + ''.join(icode) 
                for e in esc:
                    code = code.replace(" %s "%e,e)

                # Protected strings
                for p in protect:
                    code = code.replace( protect[p] , p )
                    
                new_code.append( code )
                
            # End skipdo if statement
            else:
                new_code.append( code )
            
        # End code in loopBody for loop
        if not skipdo:

            # Add in the loop end do
            for dd in range(len(my_indices)):
                indent = indent[0:-2]
                code = (comment + indent
                        + "end do\n" ) 
                new_code.append( code )

            # OMP end
            new_code.append( ompEnd )

            # Caliper profiling (END)
            if profile:
                new_code.append(caliperEnd % (self.filename,self.start))

            
            
        self.newBody = new_code
            
        
        
# DO THIS FOR EACH SUBROUTINE/FUNCTION


# Interate over lines
def_line_number = 0

unroll_defs = [ ]
unroll_list = [ ]  
unroll_item = [ ]

def findFunctions(ilines):
    # Step 0: parse subroutines
    def_line_number = 0
    prog = ["program","end program"]
    subr = ["subroutine","end subroutine"]
    func = ["function","end function"]
    modu = ["module","contains"]
    blks = [prog,subr,func,modu]
    start = MAXLINES
    myFuncs = []
    blk_i = -1
    for ii in range(len(ilines)):

        fline = ilines[ii].strip().lower()
    
        if blk_i == -1:
            for bb in range(len(blks)):
                blk = blks[bb]
                if ( fline.startswith(blk[0] ) ):
                    start = min(ii,start)
                    blk_i = bb
                    break

        if blk_i >= 0:
            if ( fline.startswith( blks[blk_i][1] ) ):
                end = ii
                kind = fline[0:3]
                ifun = iFunc(start,end,kind,os.path.basename(file_name) )
                ifun.getFuncBody(ilines)
                ifun.getHeader()
                myFuncs.append( ifun )
                start = MAXLINES
                blk_i = -1


    return myFuncs


def findLoops(ilines):
    # Step 1: parse the code
    def_line_number = 0
    start = MAXLINES
    iloops = []
    for ii in range(len(ilines)):

        fline = ilines[ii]
    
        if start_unroll in fline:
            start = min(ii,start)

        if end_unroll in fline:
            end = ii
            iloops.append(  iLoop(start,end,os.path.basename(file_name) ) )
            start = MAXLINES

    return iloops


myFuncs = findFunctions(lines)
myLoops = []

# Step 1: parse the code (loop over functions then fexl statements)
for f in myFuncs:
    start_function = f.start
    end_function = f.end
    floops = findLoops( lines[start_function:end_function] )    
    for loop in floops:
        loop.start += start_function
        loop.end += start_function
        loop.parent = f
    f.loops = floops       # Function loops
    myLoops.extend(floops) # Global loops


# Step 2: Get the contents of the FEXL loops
cnt = 1
for func in myFuncs:
    for loops in func.loops:
        loops.getLoopBody( lines )
        loops.convertLoop()
        if verbose:
            print("++++++++++ Loop # %s ++++++++++++" % cnt)
            loops.showConvert()        
        cnt += 1

# Step 3: Make function level code modifications
for func in myFuncs:    
    # Function front/end matter
    func.makeFrontMatter()
    func.makeEndMatter()
    func.makeMiddleMatter()
    func.sortLoops()


# Step 4: Actually change the code
offset = 0
for func in myFuncs:
    
    # Sort loops and loop over them
    for loops in func.loops:
        offset += loops.convert(lines,offset)
    

# Simple global find and replaces here
for ii in range(len(lines)):
    # Asynchronous macro switch
    # !$fexl-async.sync_var  -> nowait depend(sync_var)
    if ( (asyncMode or asyncAll) and (fexl_async in lines[ii])):

        try:

            if asyncAll:
                rep   = fexl_async 
                front = lines[ii].split(rep)[0]
                var   = sync_var
                back  = " \n"
            else:
                rep   = fexl_async + "."
                front = lines[ii].split(rep)[0]
                var   = lines[ii].split(rep)[1].split(".")[0]
                back  = lines[ii].split(rep)[1].split(".")[1]
                    
            wth   = ompAsync % var        
            lines[ii] = "%s %s %s" % ( front, wth, back )
            #lines[ii].replace(rep,wth) + ")"
        except:
            print("Warning: Line %s of %s is invalid syntax" % (ii,file_name) )
            print("   %s" %  lines[ii] )
    


if not replace:
    pid = open(file_name.replace(fsuff,'.FEXL-TESTING'+fsuff) ,'w')
else:
    pid = open(file_name ,'w')

pid.writelines( lines )
pid.close()

print('FEXL - Converted %s loops in %s' % (len(myLoops),file_name))

           
