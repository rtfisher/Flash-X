## This is the template for the stub Makefile generated by setup script
## All lines starting with ## are comments
## See Readme.template for more details regarding the syntax of this file
##
## The variables which may be used in this template
##
## unitname -> string
##
##
#\tMakefile generated by setup for unit %(unitname)s.  This Makefile is
#\tempty; if make complains about unresolved references at link time, you
#\tneed to create stub routines for this unit and create a Makefile
#\twhich sets the value of the macro $(%(unitname)s) equal to the list
#\tof object files containing your stub routines.  There should be a stub
#\troutine for every routine which is called outside the unit, to be
#\tused in case the unit is not included (or in case the included
#\tsub-units of this unit do not define/override some routines).

#\tFor example, if your stub routines reside in files named stub1.F,
#\tstub2.F, and stub3.F in the %(unitname)s directory, your Makefile should
#\tcontain the line
#\t\t%(unitname)s = stub1.o stub2.o stub3.o

#\tFor sub-units which redefine these routines (and give the files the
#\tsame names), no additional Makefile is required.  If a sub-unit
#\tnamed \"urp\" adds another file named myfile.F, then in %(unitname)s/urp
#\tcreate a Makefile containing

#\t\t%(unitname)s += myfile.o

# %(unitname)s =
