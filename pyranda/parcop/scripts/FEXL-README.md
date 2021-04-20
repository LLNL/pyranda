# FEXL README File
### Simple instructions for adding FEXL to array syntax code



## Simple setup
The following files need to be customized for FEXL to work with your build:
- make-funroller.py: The driver for the pre-processor and tells which files should be parsed and converted [Needs to be modified]
- ../fexl: Short cut for the forward conversion [Maybe modify]
- ../undo_fexl: Short cut for the backward conversion [Maybe modify]
- funroller.py: Actual script that will do the conversion itself.  [Maybe modify]


## Simple examples of FEXL syntax

#### Always add DEF-FEXL to include
You'll need to tell FEXL where it is appropriate to add new vaiable definition lines to the code.  If you have a FEXL loop inside
of a subroutine, you'll need to add
```fortran
SUBROUTINE foobar(arg1,arg2)
USE globals
IMPLICIT NONE
REAL*8, DIMENSION(:,:,:) :: arg1
INTEGER :: arg2
!$DEF-FEXL
...
```

#### Trivial case - Array syntax
```fortran90
!$FEXL {dim=3,var=['arg1','tmp','vari']}
tmp = arg1 * vari**2 / 3.2
vari = tmp * SQRT( arg1 )
!$FEXL
```

#### Multi-line input syntax
```fortran90
!$FEXL { dim=3,
!$FEXL   var=['arg1','tmp','vari']}
tmp = arg1 * vari**2 / 3.2
vari = tmp * SQRT( arg1 )
!$FEXL
```

#### Specify loop bounds
```fortran90
!$FEXL {dim=3,var=['arg1','tmp','vari'],
!$FEXL  bounds="x1:xn,y1:yn,z1:zn"}
tmp = arg1 * vari**2 / 3.2
vari = tmp * SQRT( arg1 )
!$FEXL
```

#### Blended array syntax
```fortran90
!$FEXL {dim=3,var=['arg1','tmp','vari']}
tmp = arg1 * vari**2 / 3.2
vari = tmp4d[:,:,:,n] * SQRT( arg1 )
!$FEXL
```

#### Find and replace
```fortran90
!$FEXL {dim=3,var=['arg1','tmp','vari'],
        replace={'tmp3d':'tmp'} }
tmp = arg1 * vari**2 / 3.2
vari = tmp3d * SQRT( arg1 )
!$FEXL
```

#### More examples to come


